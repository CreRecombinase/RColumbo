#include "RcppArmadillo.h"
#include <vector>
#include <memory>
#include <H5Cpp.h>
#include <H5Library.h>
#include "h5func.hpp"
#include "blosc_filter.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <hdf5.h>




#if defined(__GNUC__)
#define PUSH_ERR(func, minor, str, ...) H5Epush(H5E_DEFAULT, __FILE__, func, __LINE__, H5E_ERR_CLS, H5E_PLINE, minor, str, ##__VA_ARGS__)
#elif defined(_MSC_VER)
#define PUSH_ERR(func, minor, str, ...) H5Epush(H5E_DEFAULT, __FILE__, func, __LINE__, H5E_ERR_CLS, H5E_PLINE, minor, str, __VA_ARGS__)
#else
/* This version is portable but it's better to use compiler-supported
approaches for handling the trailing comma issue when possible. */
#define PUSH_ERR(func, minor, ...) H5Epush(H5E_DEFAULT, __FILE__, func, __LINE__, H5E_ERR_CLS, H5E_PLINE, minor, __VA_ARGS__)
#endif	/* defined(__GNUC__) */

#define GET_FILTER(a,b,c,d,e,f,g) H5Pget_filter_by_id(a,b,c,d,e,f,g,NULL)


size_t blosc_filter(unsigned flags, size_t cd_nelmts,
                    const unsigned cd_values[], size_t nbytes,
                    size_t *buf_size, void **buf);

herr_t blosc_set_local(hid_t dcpl, hid_t type, hid_t space);


/* Register the filter, passing on the HDF5 return value */
int register_blosc(char **version, char **date){

  int retval;

  H5Z_class_t filter_class = {
    H5Z_CLASS_T_VERS,
    (H5Z_filter_t)(FILTER_BLOSC),
    1, 1,
    "blosc",
    NULL,
    (H5Z_set_local_func_t)(blosc_set_local),
    (H5Z_func_t)(blosc_filter)
  };

  retval = H5Zregister(&filter_class);
  if(retval<0){
    PUSH_ERR("register_blosc", H5E_CANTREGISTER, "Can't register Blosc filter");
  }
  *version = strdup(BLOSC_VERSION_STRING);
  *date = strdup(BLOSC_VERSION_DATE);
  return 1; /* lib is available */
}

/*  Filter setup.  Records the following inside the DCPL:

1. If version information is not present, set slots 0 and 1 to the filter
revision and Blosc version, respectively.

2. Compute the type size in bytes and store it in slot 2.

3. Compute the chunk size in bytes and store it in slot 3.
*/
herr_t blosc_set_local(hid_t dcpl, hid_t type, hid_t space){

  int ndims;
  int i;
  herr_t r;

  unsigned int typesize, basetypesize;
  unsigned int bufsize;
  hsize_t chunkdims[32];
  unsigned int flags;
  size_t nelements = 8;
  unsigned int values[] = {0,0,0,0,0,0,0,0};
  hid_t super_type;
  H5T_class_t classt;

  r = GET_FILTER(dcpl, FILTER_BLOSC, &flags, &nelements, values, 0, NULL);
  if(r<0) return -1;

  if(nelements < 4) nelements = 4;  /* First 4 slots reserved. */

/* Set Blosc info in first two slots */
values[0] = FILTER_BLOSC_VERSION;
values[1] = BLOSC_VERSION_FORMAT;

ndims = H5Pget_chunk(dcpl, 32, chunkdims);
if(ndims<0) return -1;
if(ndims>32){
  PUSH_ERR("blosc_set_local", H5E_CALLBACK, "Chunk rank exceeds limit");
  return -1;
}

typesize = H5Tget_size(type);
if (typesize==0) return -1;
/* Get the size of the base type, even for ARRAY types */
classt = H5Tget_class(type);
if (classt == H5T_ARRAY) {
  /* Get the array base component */
  super_type = H5Tget_super(type);
  basetypesize = H5Tget_size(super_type);
  /* Release resources */
  H5Tclose(super_type);
}
else {
  basetypesize = typesize;
}

/* Limit large typesizes (they are pretty inneficient to shuffle
and, in addition, Blosc does not handle typesizes larger than
256 bytes). */
if (basetypesize > BLOSC_MAX_TYPESIZE) basetypesize = 1;
values[2] = basetypesize;

/* Get the size of the chunk */
bufsize = typesize;
for (i=0; i<ndims; i++) {
  bufsize *= chunkdims[i];
}
values[3] = bufsize;

#ifdef BLOSC_DEBUG
fprintf(stderr, "Blosc: Computed buffer size %d\n", bufsize);
#endif

r = H5Pmodify_filter(dcpl, FILTER_BLOSC, flags, nelements, values);
if(r<0) return -1;

return 1;
}


/* The filter function */
size_t blosc_filter(unsigned flags, size_t cd_nelmts,
                    const unsigned cd_values[], size_t nbytes,
                    size_t *buf_size, void **buf){

  void* outbuf = NULL;
  int status = 0;                /* Return code from Blosc routines */
size_t typesize;
size_t outbuf_size;
int clevel = 5;                /* Compression level default */
int doshuffle = 1;             /* Shuffle default */
int compcode;                  /* Blosc compressor */
int code;
char *compname = "blosclz";    /* The compressor by default */
char *complist;
char errmsg[256];

/* Filter params that are always set */
typesize = cd_values[2];      /* The datatype size */
outbuf_size = cd_values[3];   /* Precomputed buffer guess */
/* Optional params */
if (cd_nelmts >= 5) {
  clevel = cd_values[4];        /* The compression level */
}
if (cd_nelmts >= 6) {
  doshuffle = cd_values[5];  /* BLOSC_SHUFFLE, BLOSC_BITSHUFFLE */
/* bitshuffle is only meant for production in >= 1.8.0 */
#if ( (BLOSC_VERSION_MAJOR <= 1) && (BLOSC_VERSION_MINOR < 8) )
if (doshuffle == BLOSC_BITSHUFFLE) {
  PUSH_ERR("blosc_filter", H5E_CALLBACK,
           "this Blosc library version is not supported.  Please update to >= 1.8");
  goto failed;
}
#endif
}
if (cd_nelmts >= 7) {
  compcode = cd_values[6];     /* The Blosc compressor used */
/* Check that we actually have support for the compressor code */
complist = blosc_list_compressors();
code = blosc_compcode_to_compname(compcode, &compname);
if (code == -1) {
  PUSH_ERR("blosc_filter", H5E_CALLBACK,
           "this Blosc library does not have support for "
           "the '%s' compressor, but only for: %s",
           compname, complist);
  goto failed;
}
}

/* We're compressing */
if(!(flags & H5Z_FLAG_REVERSE)){

  /* Allocate an output buffer exactly as long as the input data; if
  the result is larger, we simply return 0.  The filter is flagged
  as optional, so HDF5 marks the chunk as uncompressed and
  proceeds.
  */

  outbuf_size = (*buf_size);

#ifdef BLOSC_DEBUG
  fprintf(stderr, "Blosc: Compress %zd chunk w/buffer %zd\n",
          nbytes, outbuf_size);
#endif

  outbuf = malloc(outbuf_size);

  if (outbuf == NULL){
    PUSH_ERR("blosc_filter", H5E_CALLBACK,
             "Can't allocate compression buffer");
    goto failed;
  }

  blosc_set_compressor(compname);
  status = blosc_compress(clevel, doshuffle, typesize, nbytes,
                          *buf, outbuf, nbytes);
  if (status < 0) {
    PUSH_ERR("blosc_filter", H5E_CALLBACK, "Blosc compression error");
    goto failed;
  }

  /* We're decompressing */
} else {
  /* declare dummy variables */
  size_t cbytes, blocksize;

  free(outbuf);

  /* Extract the exact outbuf_size from the buffer header.
  *
  * NOTE: the guess value got from "cd_values" corresponds to the
  * uncompressed chunk size but it should not be used in a general
  * cases since other filters in the pipeline can modify the buffere
  *  size.
  */
  blosc_cbuffer_sizes(*buf, &outbuf_size, &cbytes, &blocksize);

#ifdef BLOSC_DEBUG
  fprintf(stderr, "Blosc: Decompress %zd chunk w/buffer %zd\n", nbytes, outbuf_size);
#endif

  outbuf = malloc(outbuf_size);

  if(outbuf == NULL){
    PUSH_ERR("blosc_filter", H5E_CALLBACK, "Can't allocate decompression buffer");
    goto failed;
  }

  status = blosc_decompress(*buf, outbuf, outbuf_size);
  if(status <= 0){    /* decompression failed */
  PUSH_ERR("blosc_filter", H5E_CALLBACK, "Blosc decompression error");
    goto failed;
  } /* if !status */

} /* compressing vs decompressing */

  if(status != 0){
    free(*buf);
    *buf = outbuf;
    *buf_size = outbuf_size;
    return status;  /* Size of compressed/decompressed data */
  }

  failed:
    free(outbuf);
  return 0;

} /* End filter function */



using namespace H5;

bool file_exists (const std::string& name) {
  if (FILE *file = fopen(name.c_str(), "r")) {
    fclose(file);
    return true;
  } else {
    return false;
  }
}


H5FilePtr create_or_open_file(const std::string &fname)
{
  H5::H5File* file;
  if(!file_exists(fname)){
    try{
      file = new H5::H5File(fname.c_str(), H5F_ACC_EXCL);
    }catch(FileIException error){
      error.printError();
      Rcpp::stop("Error creating file!");
    }
  }
  else{
    try{
      file = new H5::H5File(fname.c_str(), H5F_ACC_RDWR);
    }catch(FileIException error){
      error.printError();
      Rcpp::stop("Error opening file!");
    }
  }
  return H5FilePtr(file);
}

H5FilePtr open_file(const std::string &fname)
{
  H5::H5File* file;
  if(!file_exists(fname)){
    Rcpp::stop("File does not exist!");
  }
  else{
    try{
      file = new H5::H5File(fname.c_str(), H5F_ACC_RDONLY);
    }catch(FileIException error){
      error.printError();
      Rcpp::stop("Error opening file!");
    }
  }
  return H5FilePtr(file);
}


H5GroupPtr create_or_open_group(H5FilePtr &file,const std::string &groupname)
{
  Group* group;
  Group *rg = new Group(file->openGroup("/"));
  if(groupname==""||groupname=="/"){
    return H5GroupPtr(rg);
  }
  hsize_t objc= rg->getNumObjs();
  bool fgroup=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=rg->getObjnameByIdx(i);
      if(tst==groupname){
        fgroup = true;
      }
    }
  }
  if(fgroup){
    try{
      group = new Group(file->openGroup(groupname));
    }catch(GroupIException error){
      error.printError();
      Rcpp::stop("Error opening group");
    }
  }else{
    try{
      group = new Group(file->createGroup("/"+groupname));
    }catch(GroupIException error){
      error.printError();
      Rcpp::stop("Error creating group");
    }
  }
  return H5GroupPtr(group);
}

H5GroupPtr open_group(H5FilePtr &file,const std::string &groupname)
{
  Group* group;
  Group *rg;
  try{
    rg= new Group(file->openGroup("/"));
  }catch(GroupIException error){
    error.printError();
    Rcpp::stop("Error opening group");
  }
  if(groupname==""||groupname=="/"){
    return H5GroupPtr(rg);
  }
  hsize_t objc= rg->getNumObjs();
  bool fgroup=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=rg->getObjnameByIdx(i);
      if(tst==groupname){
        fgroup = true;
      }
    }
  }
  if(fgroup){
    try{
      group = new Group(file->openGroup(groupname));
    }catch(GroupIException error){
      error.printError();
      Rcpp::stop("Error opening group");
    }
  }else{
    Rcpp::stop("Group does not exist!");
  }
  return H5GroupPtr(group);
}






// H5DataSetPtr create_or_open_dataset(H5GroupPtr &group,const std::string &dataname, const DataType &data_type,std::vector<hsize_t> &cdatadim,std::vector<hsize_t> &mdatadim,std::vector<hsize_t> &chunkdim,const unsigned int deflate_level)
// {
//   H5::Exception::dontPrint();
//   DataSet* dataset;
//   DataSpace* fdataspace;
//   hsize_t objc= group->getNumObjs();
//   bool fdat=false;
//   if(objc!=0){
//     for(hsize_t i=0; i<objc;i++){
//       std::string tst=group->getObjnameByIdx(i);
//       if(tst==dataname){
//         fdat = true;
//       }
//     }
//   }
//   if(fdat){
//     try{
//       dataset = new DataSet(group->openDataSet(dataname));
//     }  catch( DataSetIException error )
//     {
//       error.printError();
//       Rcpp::stop("Error opening dataset");
//     }
//   }else{
//     hsize_t *cumdima = new hsize_t[cdatadim.size()];
//     hsize_t *chunkdima = new hsize_t[chunkdim.size()];
//     hsize_t*mdima = new hsize_t[mdatadim.size()];
//
//     for(int i=0; i<cdatadim.size();i++){
//       cumdima[i]=cdatadim[i];
//       mdima[i]=mdatadim[i];
//       chunkdima[i]=chunkdim[i];
//     }
//     try{
//       fdataspace= new DataSpace(cdatadim.size(),cumdima,mdima); //Create dataspace for dataset (on disk)
//     }catch(DataSpaceIException error)
//     {
//       error.printError();
//       Rcpp::stop("Error creating file dataspace");
//     }
//     DSetCreatPropList cparms; //Create chunksize file parameters
//     cparms.setChunk(chunkdim.size(),chunkdima); //Set chunksize
//     cparms.setDeflate(deflate_level);
//     try{
//       dataset = new DataSet(group->createDataSet(dataname,data_type,*fdataspace,cparms));
//     }  catch( DataSetIException error )
//     {
//       error.printError();
//       Rcpp::stop("Error creating dataset");
//     }
//   }
//   return H5DataSetPtr(dataset);
// }



H5DataSetPtr create_or_open_dataset(H5GroupPtr &group,const std::string &dataname, const DataType &data_type,std::vector<hsize_t> &cdatadim,std::vector<hsize_t> &mdatadim,std::vector<hsize_t> &chunkdim,const int deflate_level)
{
  H5::Exception::dontPrint();
  DataSet* dataset;
  DataSpace* fdataspace;
  hsize_t objc= group->getNumObjs();



  char* version;
  char* date;
  int r=0;
  unsigned int cd_values[7];
  if(deflate_level>0){
    r = register_blosc(&version,&date);

//    printf("Blosc version info: %s (%s)\n", version, date);

    cd_values[4]=deflate_level;
    cd_values[5]=1;
    cd_values[6]=BLOSC_BLOSCLZ;
  }
  bool fdat=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=group->getObjnameByIdx(i);
      if(tst==dataname){
        fdat = true;
      }
    }
  }
  if(fdat){
    try{
      dataset = new DataSet(group->openDataSet(dataname));
    }  catch( DataSetIException error )
    {
      error.printError();
      Rcpp::stop("Error opening dataset");
    }
  }else{
    hsize_t *cumdima = new hsize_t[cdatadim.size()];
    hsize_t *chunkdima = new hsize_t[chunkdim.size()];
    hsize_t *chunksizea = new hsize_t[chunkdim.size()];
    hsize_t*mdima = new hsize_t[mdatadim.size()];

    for(int i=0; i<cdatadim.size();i++){
      cumdima[i]=cdatadim[i];
      mdima[i]=mdatadim[i];
      chunkdima[i]=chunkdim[i];
      hsize_t tchunk=std::min((int)chunkdima[i],1000);
      chunksizea[i]=tchunk;
    }
    try{
      fdataspace= new DataSpace(cdatadim.size(),cumdima,mdima); //Create dataspace for dataset (on disk)
    }catch(DataSpaceIException error)
    {
      error.printError();
      Rcpp::stop("Error creating file dataspace");
    }
    DSetCreatPropList cparms; //Create chunksize file parameters
    cparms.setChunk(chunkdim.size(),chunkdima); //Set chunksize
  if(deflate_level>0){
    cparms.setFilter(FILTER_BLOSC,H5Z_FLAG_OPTIONAL,7,cd_values);
  }else{
//    std::cout<<"Using DEFLATE for compression"<<std::endl;
    cparms.setDeflate(-deflate_level);
//    cparms.setFilter(FILTER_BLOSC,H5Z_FLAG_OPTIONAL,7,cd_values);
  }
    try{
      dataset = new DataSet(group->createDataSet(dataname,data_type,*fdataspace,cparms));
    }  catch( DataSetIException error )
    {
      error.printError();
      Rcpp::stop("Error creating dataset");
    }
  }
  return H5DataSetPtr(dataset);
}



H5DataSetPtr create_or_open_ref_dataset(H5GroupPtr &group,const std::string &dataname,std::vector<hsize_t> &cdatadim,std::vector<hsize_t> &mdatadim)
{
  H5::Exception::dontPrint();
  DataSet* dataset;
  DataSpace* fdataspace;
  hsize_t objc= group->getNumObjs();
  std::string refdataname=dataname+"_ref";

  bool fdat=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=group->getObjnameByIdx(i);
      if(tst==refdataname){
        fdat = true;
      }
    }
  }
  if(fdat){
    try{
      dataset = new DataSet(group->openDataSet(refdataname));
    }  catch( DataSetIException error )
    {
      error.printError();
      Rcpp::stop("Error opening dataset");
    }
  }else{
    //    std::cout<<"Creating dataset"<<std::endl;
    hsize_t *cumdima = new hsize_t[cdatadim.size()];
    hsize_t*mdima = new hsize_t[mdatadim.size()];

    for(int i=0; i<cdatadim.size();i++){
      cumdima[i]=cdatadim[i];
      mdima[i]=mdatadim[i];
    }
    try{
//      std::cout<<"Creating dataspace of dimensions"<<cumdima[0]<<" up to "<<mdima[0]<<std::endl;
      fdataspace= new DataSpace(cdatadim.size(),cumdima,mdima); //Create dataspace for dataset (on disk)
    }catch(DataSpaceIException error)
    {
      error.printError();
      Rcpp::stop("Error creating file dataspace");
    }
    try{
//      std::cout<<"Creating datatype: "<<refdataname<<std::endl;
      DataType dt(PredType::STD_REF_DSETREG);
//      std::cout<<"Creating dataset: "<<refdataname<<std::endl;
      dataset = new DataSet(group->createDataSet(refdataname,dt,*fdataspace));
    }  catch( DataSetIException error )
      {
	error.printError();
	Rcpp::stop("Error creating dataset");
      }
//      std::cout<<"Dataset created: "<<refdataname<<std::endl;
  }
  return H5DataSetPtr(dataset);
}


H5DataSetPtr open_ref_dataset(H5GroupPtr &group,const std::string &dataname)
{
  H5::Exception::dontPrint();
  DataSet* dataset;
  DataSpace* fdataspace;
  hsize_t objc= group->getNumObjs();
  std::string refdataname=dataname+"_ref";

  bool fdat=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=group->getObjnameByIdx(i);
      if(tst==refdataname){
        fdat = true;
      }
    }
  }
  if(fdat){
    try{
      dataset = new DataSet(group->openDataSet(refdataname));
    }  catch( DataSetIException error )
    {
      error.printError();
      Rcpp::stop("Error opening dataset");
    }
  }else{
    Rcpp::stop("Error opening dataset reference dataset");
  }
  return H5DataSetPtr(dataset);
}



H5DataSetPtr open_blosc_dataset(H5GroupPtr &group,const std::string &dataname)
{
  H5::Exception::dontPrint();
  DataSet* dataset;
  DataSpace* fdataspace;
  hsize_t objc= group->getNumObjs();
  char* version;
  char* date;
  int r=0;
//  std::cout<<"registering blosc"<<std::endl;
  r = register_blosc(&version,&date);
  unsigned int cd_values[7];
//  printf("Blosc version info: %s (%s)\n", version, date);
  // cd_values[4]=deflate_level;
  // cd_values[5]=1;
  // cd_values[6]=BLOSC_BLOSCLZ;
  bool fdat=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=group->getObjnameByIdx(i);
      if(tst==dataname){
        fdat = true;
      }
    }
  }
  if(fdat){
    try{
      dataset = new DataSet(group->openDataSet(dataname));
    }  catch( DataSetIException error )
    {
      error.printError();
      Rcpp::stop("Error opening dataset");
    }
  }else{
    Rcpp::stop("Dataset not found!");
  }
  return H5DataSetPtr(dataset);
}










H5DataSetPtr get_dataset(const std::string h5filename, const std::string groupname, const std::string dataname){
  H5File* file;
  DataSet* dataset;
  try{
    file= new H5File( h5filename, H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
    Rcpp::stop("FileIexception!");
  }
  Group haplogroup = file->openGroup(groupname);
  try{
    dataset = new DataSet(haplogroup.openDataSet(dataname));
  }
  catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("DataSetIException!");
  }
  return H5DataSetPtr(dataset);
}



ArrayType read_arraytype(H5DataSetPtr dataset, const PredType pt){
  hsize_t adim[1];
  DataType dt = dataset->getDataType();
  DataType dtss = dt.getSuper();
  size_t dts=dt.getSize();
  size_t bcs =dtss.getSize();
  hsize_t dims= (dts/bcs);
  adim[0]=dims;
  ArrayType mem_arraytype(pt,1,adim);
  return(mem_arraytype);
}




hsize_t get_arraysize(ArrayType &atype){
  hsize_t asize[1];
  atype.getArrayDims(asize);
  return(asize[0]);
}

//[[Rcpp::export]]
size_t data_bytes(const std::string h5file,const std::string groupname, const std::string dataname){
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group = open_group(file,groupname);
  H5DataSetPtr dataset = open_blosc_dataset(group,dataname);
  DataType dt=dataset->getDataType();
  size_t bytes=dt.getSize();
  dt.close();
  dataset->close();
  group->close();
  file->close();
  return(bytes);
}

//[[Rcpp::export]]
size_t get_rownum_h5(const std::string h5file,const std::string groupname, const std::string dataname){
  H5FilePtr file = open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_blosc_dataset(group,dataname);
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  fspace.close();
  dataset->close();
  file->close();
  return(datadim[0]);
}



//[[Rcpp::export]]
size_t get_colnum_h5(const std::string h5file,const std::string groupname, const std::string dataname){
  H5FilePtr file = open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_blosc_dataset(group,dataname);
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[]={0,0};;
  fspace.getSimpleExtentDims(datadim,NULL);
  fspace.close();
  dataset->close();
  file->close();
  return(datadim[1]);
}


bool isArrayType(H5DataSetPtr dataset){
//  Rcpp::Rcout<<"Getting array type"<<std::endl;
  hid_t did=dataset->getId();
  hid_t dtype = H5Dget_type(did);
  if(H5T_ARRAY==H5Tget_class(dtype)){
    return(true);
  }else{
    return(false);
  }
}

bool isStrType(H5DataSetPtr dataset){
  hid_t did=dataset->getId();
  hid_t dtype = H5Dget_type(did);
  if(H5T_STRING==H5Tget_class(dtype)){
    return(true);
  }else{
    return(false);
  }
}

bool isNumeric(H5DataSetPtr dataset){
  hid_t did=dataset->getId();
  hid_t dtype = H5Dget_type(did);
  if(H5T_FLOAT==H5Tget_class(dtype)){
    return(true);
  }else{
    return(false);
  }
}

bool isInteger(H5DataSetPtr dataset){
  hid_t did=dataset->getId();
  hid_t dtype = H5Dget_type(did);
  if(H5T_INTEGER==H5Tget_class(dtype)){
    return(true);
  }else{
    return(false);
  }
}
bool isFloatArray(H5DataSetPtr dataset){
  return(dataset->getArrayType().getSuper().getClass()==H5T_FLOAT);
}

bool isIntArray(H5DataSetPtr dataset){
  return(dataset->getArrayType().getSuper().getClass()==H5T_INTEGER);
}



size_t copy_subset(const std::string oh5file,const std::string nh5file,const std::string groupname, const std::string dataname,const arma::uvec index){

  H5DataSetPtr dataset=get_dataset(oh5file,groupname,dataname);
  if(isArrayType(dataset)){
    if(isFloatArray(dataset)){
//      Rcpp::Rcout<<"Reading /"<<groupname<<"/"<<dataname<<" ..."<<std::endl;
      arma::fmat subdata= read_fmat_chunk_ind(oh5file,groupname,dataname,index);
      write_mat_h5(nh5file,groupname,dataname,subdata.n_cols,subdata.n_rows,subdata,2);
    }else{
      Rcpp::stop("copying non-float matrices not yet supported");
      // if(isIntArray(dataset)){
      //   Rcpp::Rcout<<"Reading /"<<groupname<<"/"<<dataname<<" ..."<<std::endl;
      //   arma::fmat subdata= read_fmat_chunk_ind(oh5file,groupname,dataname,index);
      //   write_mat_h5(nh5file,groupname,dataname,subdata.n_cols,subdata.n_rows,subdata,2);
      // }
    }
  }else{
    if(isInteger(dataset)){
      arma::uvec tdata =arma::conv_to<arma::uvec>::from(read_uint_h5(oh5file,groupname,dataname));
      tdata=tdata.elem(index-1);
      write_uint_h5(nh5file,groupname,dataname,tdata.n_elem,tdata,2);
    }else{
      if(isNumeric(dataset)){
        arma::fvec tdata = arma::conv_to<arma::fvec>::from(read_float_h5(oh5file,groupname,dataname));
        tdata = tdata.elem(index-1);
        write_float_h5(nh5file,groupname,dataname,tdata.n_elem,tdata,2);
      }
    }
  }
return(1);
}

//[[Rcpp::export]]
std::vector<std::string> getGroups(const std::string h5file){

  H5File* file;

  try{
    file= new H5File( h5file.c_str(), H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
  }
  // file->getNumObjs()
  // Group group = file->openGroup(groupname);
  hsize_t objc= file->getNumObjs();
  std::vector<std::string> retvec(objc);
  bool fdat=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=file->getObjnameByIdx(i);
      retvec[i]=tst;
    }
  }
  file->close();
  return(retvec);
}


//[[Rcpp::export]]
std::vector<std::string> getObjects(const std::string h5file, const std::string groupname){

  H5File* file;

  try{
    file= new H5File( h5file.c_str(), H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
  }
  Group group = file->openGroup(groupname);
  hsize_t objc= group.getNumObjs();
  std::vector<std::string> retvec(objc);
  bool fdat=false;
  if(objc!=0){
    for(hsize_t i=0; i<objc;i++){
      std::string tst=group.getObjnameByIdx(i);
      retvec[i]=tst;
    }
  }
  group.close();
  file->close();
  return(retvec);
}


//[[Rcpp::export]]
size_t get_arraydim(const std::string h5file, const std::string groupname,const std::string dataname){

  H5File* file;
  DataSet* dataset;
  try{
    file= new H5File( h5file.c_str(), H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
  }
  Group haplogroup = file->openGroup(groupname);
  try{
    dataset = new DataSet(haplogroup.openDataSet(dataname));
  }
  catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("DataSetIException!");
  }
  DataSpace fspace =dataset->getSpace();
  DataType dt = dataset->getDataType();
  hsize_t adim[1];
  DataType dtss = dt.getSuper();
  size_t dts=dt.getSize();
  size_t bcs =dtss.getSize();
  hsize_t dims= (dts/bcs);
  adim[0]=dims;

  dataset->close();
  dt.close();
  dtss.close();
  dataset->close();
  haplogroup.close();
  fspace.close();

  return(dims);
}

H5::DataType read_data_type(const std::string h5file, const std::string groupname, const std::string dataname){
  H5File* file;
  DataSet* dataset;
  try{
    file= new H5File( h5file.c_str(), H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
    Rcpp::stop("FileIexception!");
  }
  Group haplogroup = file->openGroup(groupname);
  try{
    dataset = new DataSet(haplogroup.openDataSet(dataname));
  }
  catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("DataSetIException!");
  }
  DataSpace fspace =dataset->getSpace();
  DataType dt = dataset->getDataType();

  dataset->close();
  dataset->close();
  haplogroup.close();
  fspace.close();
  //  Rcpp::Rcout<<dt.getSize()<<std::endl;
  return(dt);
}


//[[Rcpp::export]]
arma::fmat read_fmat_h5(const std::string hap_h5file,const std::string groupname, const std::string dataname, size_t offset, size_t chunksize){

  H5FilePtr file=open_file(hap_h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset=open_blosc_dataset(group,dataname);
  // std::cout<<"opened dataset"<<std::endl;
  // std::cout<<"Getting array type"<<std::endl;
  ArrayType mem_arraytype=read_arraytype(dataset,PredType::NATIVE_FLOAT);


  // std::cout<<"Getting Array dimensions"<<std::endl;
  // std::cout<<"Getting Dataspace"<<std::endl;
  DataSpace fspace =dataset->getSpace();

//  std::cout<<"Getting Data dimensions"<<std::endl;
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  size_t matrix_dims[2];
  if(offset+chunksize>datadim[0]){
    chunksize=datadim[0]-offset;
  }
  datadim[0]=chunksize;
  matrix_dims[0]=get_arraysize(mem_arraytype);
  matrix_dims[1]=chunksize;


  //   std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadim,offseta);
  //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  arma::fmat retmat(matrix_dims[0],matrix_dims[1]);
  DataSpace memspace(1,datadim);
  //    std::cout<<"Reading data"<<std::endl;
  dataset->read(retmat.memptr(),mem_arraytype,memspace,fspace);
  mem_arraytype.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  return(retmat);
}

//[[Rcpp::export]]
arma::mat read_dmat_h5(const std::string hap_h5file,const std::string groupname, const std::string dataname, const size_t offset, const size_t chunksize){


  H5FilePtr file=open_file(hap_h5file);
  H5GroupPtr group=open_group(file,groupname);
  H5DataSetPtr dataset=open_blosc_dataset(group,dataname);


  ArrayType mem_arraytype=read_arraytype(dataset,PredType::NATIVE_DOUBLE);
  //  std::cout<<"Getting Array dimensions"<<std::endl;
  //  std::cout<<"Getting Dataspace"<<std::endl;
  DataSpace fspace =dataset->getSpace();

  //  std::cout<<"Getting Data dimensions"<<std::endl;
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  size_t matrix_dims[2];
  size_t nchunksize=chunksize;
  if(offset+chunksize>datadim[0]){
    nchunksize=datadim[0]-offset;
  }
  datadim[0]=nchunksize;
  matrix_dims[0]=get_arraysize(mem_arraytype);
  matrix_dims[1]=nchunksize;


  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadim,offseta);
//  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  arma::mat retmat(matrix_dims[0],matrix_dims[1]);
  DataSpace memspace(1,datadim);
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read(retmat.memptr(),mem_arraytype,memspace,fspace);
  mem_arraytype.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  return(retmat);
}


//[[Rcpp::export]]
arma::mat read_dmat_chunk_ind(const std::string h5file,const std::string groupname, const std::string dataname, const arma::uvec indvec){
  arma::uword minind=arma::min(indvec);
  arma::uword maxind=arma::max(indvec);
  arma::uword chunksize=maxind-minind+1;
  arma::uvec newind = indvec-minind;
  arma::mat retmat=read_dmat_h5(h5file,groupname,dataname,minind-1,chunksize);
  return(retmat.cols(newind));
}

//[[Rcpp::export]]
arma::fmat read_fmat_chunk_ind(const std::string h5file,const std::string groupname, const std::string dataname, const arma::uvec indvec){
  arma::uword minind=arma::min(indvec);
  arma::uword maxind=arma::max(indvec);
  arma::uword chunksize=maxind-minind+1;
  arma::uvec newind = indvec-minind;
  arma::fmat retmat=read_fmat_h5(h5file,groupname,dataname,minind-1,chunksize);
  return(retmat.cols(newind));
}




std::vector<int> read_int_h5(const std::string h5file, const std::string groupname, const std::string dataname){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group=open_group(file,groupname);
  H5DataSetPtr dataset=open_blosc_dataset(group,dataname);


  DataType dt = dataset->getDataType();
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  hsize_t readrows=datadim[0];
  DataSpace memspace(1,datadim);
  std::vector<int> retdat(datadim[0]);

//  std::cout<<"Reading in int dataset dataname"<<std::endl;

  dataset->read(&retdat[0],dt,memspace,fspace);
  fspace.close();
  memspace.close();
  dataset->close();
  file->close();

//  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
  return(retdat);

}


//
//
// template <int RTYPE> Vector<RTYPE> read_h5_col(const std::string h5file, const std::string groupname, const std::string dataname, const std::string dclass, const arma::uvec col_index){
//
//   using namespace Rcpp;
//   const size_t indnum=col_index.n_elem;
//
//   H5FilePtr file=open_file(h5file);
//   H5GroupPtr group=open_group(file,groupname);
//   H5DataSetPtr dataset =open_blosc_dataset(group,dataname);
//   DataType dt = dataset->getDataType();
//   DataSpace fspace =dataset->getSpace();
//
//   hsize_t datadim[]={0,0};
//   fspace.getSimpleExtentDims(datadim,NULL);
//
//   hsize_t readrows=datadim[0];
//   if(indnum==0){
//     std::cout<<"Reading entire dataset"<<std::endl;
//     if(datadim[0]==1&&datadim[1]!=0){
//       fspace.close();
//       dataset->close();
//       file->close();
//       DataSpace memspace(1,datadim);
//       if(dclass=="FLOAT"){
//         if(dt.getSize()==4){
//           arma::fvec retdat=arma::vectorise(read_fmat_h5(h5file,groupname,dataname,0,datadim[1]));
//           memspace.close();
//           return(Rcpp::DataFrame::create(_[dataname]= Rcpp::NumericVector(retdat.begin(),retdat.end())));
//         }else{
//           arma::fvec retdat =arma::conv_to<arma::fvec>::from(arma::vectorise(read_dmat_h5(h5file,groupname,dataname,0,datadim[1])));
//           memspace.close();
//           return(Rcpp::DataFrame::create(_[dataname]= Rcpp::NumericVector(retdat.begin(),retdat.end())));
//         }
//       }else{
//         Rcpp::stop("reading matrix class of given data type is not yet supported");
//       }
//     }
//     if(dclass=="FLOAT"){
//       std::vector<float> retdat(datadim[0]);
//       //std::cout<<"Reading in float dataset "<<dataname<<std::endl;
//       DataSpace memspace(1,datadim);
//       dataset->read(&retdat[0],PredType::NATIVE_FLOAT,memspace,fspace);
//       fspace.close();
//       memspace.close();
//       dataset->close();
//       file->close();
//       //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
//       return(Rcpp::DataFrame::create(_[dataname]= Rcpp::NumericVector(retdat.begin(),retdat.end())));
//     }else{
//       if(dclass=="INTEGER"){
//         std::vector<int> retdat(datadim[0]);
//         //      std::cout<<"Reading in float dataset "<<dataname<<std::endl;
//         DataSpace memspace(1,datadim);
//         dataset->read(&retdat[0],PredType::NATIVE_INT,memspace,fspace);
//         fspace.close();
//         memspace.close();
//         dataset->close();
//         file->close();
//         //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
//         return(Rcpp::DataFrame::create(_[dataname]= Rcpp::IntegerVector(retdat.begin(),retdat.end())));
//       }
//     }
//   }else{
//
//     arma::sword min_ind= arma::min(col_index);
//     arma::sword max_ind= arma::max(col_index);
//
//     arma::uword chunksize=max_ind-min_ind+1;
//     std::cout<<"Reading subset of size: "<<chunksize<<" from "<<min_ind<<" to "<<max_ind<<std::endl;
//     arma::uvec scol_index = col_index-min_ind;
//     hsize_t offseta[1];
//     offseta[0] =min_ind-1;
//     hsize_t memdim[1];
//     memdim[0]=chunksize;
//     fspace.selectHyperslab(H5S_SELECT_SET,memdim,offseta);
//
//     //   std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
//     //    fspace.selectElements(H5S_SELECT_SET,indnum,scol_index.memptr());
//     //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
//     DataSpace memspace(1,memdim);
//
//     if(dclass=="FLOAT"){
//       arma::fvec retdat(chunksize);
//       //      std::cout<<"Reading float data"<<std::endl;
//       dataset->read(retdat.memptr(),PredType::NATIVE_FLOAT,memspace,fspace);
//       retdat = retdat.elem(scol_index);
//       memspace.close();
//       fspace.close();
//       group->close();
//       dataset->close();
//       file->close();
//       return(Rcpp::DataFrame::create(_[dataname]= Rcpp::NumericVector(retdat.begin(),retdat.end())));
//     }else{
//       if(dclass=="INTEGER"){
//         arma::Col<int> retdat(chunksize,arma::fill::zeros);
//         std::cout<<"Reading in integer dataset "<<dataname<<std::endl;
//         dataset->read(retdat.memptr(),PredType::NATIVE_INT,memspace,fspace);
//         std::cout<<"Subsetting integer dataset "<<dataname<<std::endl;
//         retdat=retdat.elem(scol_index);
//         fspace.close();
//         memspace.close();
//         dataset->close();
//         file->close();
//         return(Rcpp::DataFrame::create(_[dataname]= Rcpp::IntegerVector(retdat.begin(),retdat.end())));
//       }else{
//         Rcpp::stop("reading matrix class of given data type is not yet supported");
//       }
//     }
//   }
// }



// }

//[[Rcpp::export]]
Rcpp::DataFrame read_h5_df_col(const std::string h5file, const std::string groupname, const std::string dataname,const std::string dclass,const arma::uvec col_index){

  using namespace Rcpp;
  const size_t indnum=col_index.n_elem;

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group=open_group(file,groupname);
  H5DataSetPtr dataset =open_blosc_dataset(group,dataname);
  DataType dt = dataset->getDataType();
  DataSpace fspace =dataset->getSpace();

  hsize_t datadim[]={0,0};
  fspace.getSimpleExtentDims(datadim,NULL);

  hsize_t readrows=datadim[0];
  if(indnum==0){
    std::cout<<"Reading entire dataset"<<std::endl;
    if(datadim[0]==1&&datadim[1]!=0){
      fspace.close();
      dataset->close();
      file->close();
      DataSpace memspace(1,datadim);
      if(dclass=="FLOAT"){
        if(dt.getSize()==4){
          arma::fvec retdat=arma::vectorise(read_fmat_h5(h5file,groupname,dataname,0,datadim[1]));
          memspace.close();
          return(Rcpp::DataFrame::create(_[dataname]= Rcpp::NumericVector(retdat.begin(),retdat.end())));
        }else{
          arma::fvec retdat =arma::conv_to<arma::fvec>::from(arma::vectorise(read_dmat_h5(h5file,groupname,dataname,0,datadim[1])));
          memspace.close();
          return(Rcpp::DataFrame::create(_[dataname]= Rcpp::NumericVector(retdat.begin(),retdat.end())));
        }
      }else{
        Rcpp::stop("reading matrix class of given data type is not yet supported");
      }
    }
    if(dclass=="FLOAT"){

      std::vector<float> retdat(datadim[0]);
      //std::cout<<"Reading in float dataset "<<dataname<<std::endl;
      DataSpace memspace(1,datadim);
      dataset->read(&retdat[0],PredType::NATIVE_FLOAT,memspace,fspace);
      fspace.close();
      memspace.close();
      dataset->close();
      file->close();
      //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
      return(Rcpp::DataFrame::create(_[dataname]= Rcpp::NumericVector(retdat.begin(),retdat.end())));
    }else{
      if(dclass=="INTEGER"){
        std::vector<int> retdat(datadim[0]);
        //      std::cout<<"Reading in float dataset "<<dataname<<std::endl;
        DataSpace memspace(1,datadim);
        dataset->read(&retdat[0],PredType::NATIVE_INT,memspace,fspace);
        fspace.close();
        memspace.close();
        dataset->close();
        file->close();
        //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
        return(Rcpp::DataFrame::create(_[dataname]= Rcpp::IntegerVector(retdat.begin(),retdat.end())));
      }
    }
  }else{

    arma::sword min_ind= arma::min(col_index);
    arma::sword max_ind= arma::max(col_index);
    // std:vector<float> retdat
    arma::uword chunksize=max_ind-min_ind+1;
    std::cout<<"Reading subset of size: "<<chunksize<<" from "<<min_ind<<" to "<<max_ind<<std::endl;
    arma::uvec scol_index = col_index-min_ind;
    hsize_t offseta[1];
    offseta[0] =min_ind-1;
    hsize_t memdim[1];
    memdim[0]=chunksize;
    fspace.selectHyperslab(H5S_SELECT_SET,memdim,offseta);

    //   std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
    //    fspace.selectElements(H5S_SELECT_SET,indnum,scol_index.memptr());
    //  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
    DataSpace memspace(1,memdim);

    if(dclass=="FLOAT"){
      arma::fvec retdat(chunksize);
      //      std::cout<<"Reading float data"<<std::endl;
      dataset->read(retdat.memptr(),PredType::NATIVE_FLOAT,memspace,fspace);
      retdat = retdat.elem(scol_index);
      memspace.close();
      fspace.close();
      group->close();
      dataset->close();
      file->close();
      return(Rcpp::DataFrame::create(_[dataname]= Rcpp::NumericVector(retdat.begin(),retdat.end())));
    }else{
      if(dclass=="INTEGER"){
        arma::Col<int> retdat(chunksize,arma::fill::zeros);
        std::cout<<"Reading in integer dataset "<<dataname<<std::endl;
        dataset->read(retdat.memptr(),PredType::NATIVE_INT,memspace,fspace);
        std::cout<<"Subsetting integer dataset "<<dataname<<std::endl;
        retdat=retdat.elem(scol_index);
        fspace.close();
        memspace.close();
        dataset->close();
        file->close();
        return(Rcpp::DataFrame::create(_[dataname]= Rcpp::IntegerVector(retdat.begin(),retdat.end())));
      }else{
        Rcpp::stop("reading matrix class of given data type is not yet supported");
      }
    }
  }
}




//[[Rcpp::export]]
Rcpp::IntegerVector read_Rint_h5(const std::string h5file, const std::string groupname, const std::string dataname){
  std::vector<int> retvec=read_int_h5(h5file,groupname,dataname);
  Rcpp::IntegerVector tretvec(retvec.begin(),retvec.end());
  return(tretvec);
}


//[[Rcpp::export]]
std::vector<float> read_float_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t offset,const size_t chunksize){
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group=open_group(file,groupname);
  H5DataSetPtr dataset=open_blosc_dataset(group,dataname);
  DataType dt = dataset->getDataType();
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  if(chunksize+offset>datadim[0]){
    Rcpp::Rcerr<<"Attempting to access element past extent of data ("<<chunksize+offset<<">"<<datadim[0]<<")"<<std::endl;
    Rcpp::stop("Error reading floatvec");
  }
  hsize_t readrows=chunksize;
  datadim[0]=chunksize;
  hsize_t offseta[1];
  offseta[0]=offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadim,offseta);
  DataSpace memspace(1,datadim);
  std::vector<float> retdat(chunksize);

//  std::cout<<"Reading in float dataset dataname"<<std::endl;

  dataset->read(&retdat[0],dt,memspace,fspace);
  fspace.close();
  memspace.close();
  dataset->close();
  file->close();
  //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
  return(retdat);
}

//[[Rcpp::export]]
std::vector<double> read_double_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t offset,const size_t chunksize){
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group=open_group(file,groupname);
  H5DataSetPtr dataset=open_blosc_dataset(group,dataname);
  DataType dt = dataset->getDataType();
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  if(chunksize+offset>datadim[0]){
    Rcpp::Rcerr<<"Attempting to access element past extent of data ("<<chunksize+offset<<">"<<datadim[0]<<")"<<std::endl;
    Rcpp::stop("Error reading floatvec");
  }
  hsize_t readrows=chunksize;
  datadim[0]=chunksize;
  hsize_t offseta[1];
  offseta[0]=offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadim,offseta);
  DataSpace memspace(1,datadim);
  std::vector<double> retdat(chunksize);

//  std::cout<<"Reading in float dataset dataname"<<std::endl;

  dataset->read(&retdat[0],dt,memspace,fspace);
  fspace.close();
  memspace.close();
  dataset->close();
  file->close();
  //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
  return(retdat);
}





std::vector<float> read_float_h5(const std::string h5file, const std::string groupname, const std::string dataname){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group=open_group(file,groupname);
  H5DataSetPtr dataset=open_blosc_dataset(group,dataname);
  DataType dt = dataset->getDataType();
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  hsize_t readrows=datadim[0];
  DataSpace memspace(1,datadim);
  std::vector<float> retdat(datadim[0]);
  //  std::cout<<"Reading in int dataset dataname"<<std::endl;
  dataset->read(&retdat[0],dt,memspace,fspace);
  fspace.close();
  memspace.close();
  dataset->close();
  file->close();
  //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
  return(retdat);
}

//[[Rcpp::export]]
Rcpp::NumericVector read_Rfloat_h5(const std::string h5file, const std::string groupname, const std::string dataname){
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group=open_group(file,groupname);
  H5DataSetPtr dataset=open_blosc_dataset(group,dataname);
  DataType dt = dataset->getDataType();
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[]={0,0};
  fspace.getSimpleExtentDims(datadim,NULL);
  hsize_t readrows=datadim[0];
  std::cout<<readrows<<std::endl;
  DataSpace memspace(1,datadim);
  if(dt.getSize()==4){
    if(datadim[0]==1&&datadim[1]!=0){
      fspace.close();
      memspace.close();
      dataset->close();
      file->close();
      arma::fvec retdat =arma::vectorise(read_fmat_h5(h5file,groupname,dataname,0,datadim[1]));
      Rcpp::NumericVector retvec(retdat.begin(),retdat.end());
      //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
      return(retvec);
    }
    std::vector<float> retdat(datadim[0]);
    //      std::cout<<"Reading in float dataset "<<dataname<<std::endl;
    dataset->read(&retdat[0],dt,memspace,fspace);
    fspace.close();
    memspace.close();
    dataset->close();
    file->close();
    Rcpp::NumericVector retvec(retdat.begin(),retdat.end());
    //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
    return(retvec);
  }else{
    if(datadim[0]==1&&datadim[1]!=0){
      fspace.close();
      memspace.close();
      dataset->close();
      file->close();
      arma::vec retdat =arma::vectorise(read_dmat_h5(h5file,groupname,dataname,0,datadim[1]));
      Rcpp::NumericVector retvec(retdat.begin(),retdat.end());
      //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
      return(retvec);
    }
    std::vector<double> retdat(datadim[0]);
    //    std::cout<<"Reading in double dataset"<<dataname<<std::endl;
    dataset->read(&retdat[0],dt,memspace,fspace);
    fspace.close();
    memspace.close();
    dataset->close();
    file->close();
    Rcpp::NumericVector retvec(retdat.begin(),retdat.end());
    //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
    return(retvec);
  }
}



std::vector<unsigned int> read_uint_h5(const std::string h5file, const std::string groupname, const std::string dataname){

  H5FilePtr file=open_file(h5file);
  H5GroupPtr group=open_group(file,groupname);
  H5DataSetPtr dataset=open_blosc_dataset(group,dataname);

  DataType dt = dataset->getDataType();
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  hsize_t readrows=datadim[0];
  DataSpace memspace(1,datadim);
  std::vector<unsigned int> retdat(datadim[0]);

  //  std::cout<<"Reading in int dataset dataname"<<std::endl;
  dataset->read(&retdat[0],dt,memspace,fspace);
  fspace.close();
  memspace.close();
  dataset->close();
  file->close();
  //  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
  return(retdat);
}


//[[Rcpp::export]]
arma::uvec intersect_col(const std::string h5file1, const std::string h5groupname1, const std::string h5colname1, const std::string h5file2, const std::string h5groupname2, const std::string h5colname2){
  using namespace Rcpp;
//  Rcout<<"Reading col1"<<std::endl;
  std::vector<unsigned int> col1=read_uint_h5(h5file1,h5groupname1,h5colname1);
//  Rcout<<"Sorting col1"<<std::endl;
  std::sort(col1.begin(),col1.end());
//  Rcout<<"Reading col2"<<std::endl;
  std::vector<unsigned int> col2=read_uint_h5(h5file2,h5groupname2,h5colname2);
//  Rcout<<"Sorting col2"<<std::endl;
  std::sort(col2.begin(),col2.end());
//  Rcout<<"Preallocating intersection"<<std::endl;
  std::vector<unsigned int> intersection(std::max(col1.size(),col2.size()));
  std::vector<unsigned int>::iterator it;
//  Rcout<<"Computing intersection"<<std::endl;
  it =std::set_intersection(col1.begin(),col1.end(),col2.begin(),col2.end(),intersection.begin());
//  Rcout<<"Resizing results of size"<<it-intersection.begin()<<std::endl;
  intersection.resize(it-intersection.begin());
  return(arma::conv_to<arma::uvec>::from(intersection));
}




//[[Rcpp::export]]
arma::fmat read_2dfmat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t row_offset,const size_t col_offset,const size_t row_chunksize,const size_t col_chunksize){
  //Try breaking up reads in to chunks
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_blosc_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();
  int nmdc;
  size_t cache_elem=0;
  size_t cache_bytes=0;
  //  unsigned int majv,minv,relv;
  //  H5Library::getLibVersion(majv,minv,relv);
  //  std::cout<<"HDF5 Version:"<<majv<<"."<<minv<<"."<<relv<<std::endl;
  //  double policy;
  // FileAccPropList plist= file->getAccessPlist();
  // plist.getCache(nmdc,cache_elem,cache_bytes,policy);
  // std::cout<<"Metadata cache has  "<<nmdc<<std::endl;
  // std::cout<<"Elements in raw data chunk cache  is"<<cache_elem<<std::endl;
  // std::cout<<"Bytes in raw data chunk cache  is"<<cache_bytes<<std::endl;
  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0,0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
  hsize_t matrix_dims[2];
  size_t rchunksize=row_chunksize;
  size_t cchunksize=col_chunksize;
  if(row_offset+rchunksize>datadims[1]){
    rchunksize=datadims[1]-row_offset;
  }
  if(col_offset+cchunksize>datadims[0]){
    cchunksize=datadims[0]-col_offset;
  }
  hsize_t filechunksize[2];
  cparms.getChunk(2,filechunksize);
  size_t colchunknum=std::ceil((double)cchunksize/(double)filechunksize[0]);
  size_t rowchunknum=std::ceil((double)rchunksize/(double)filechunksize[1]);
  // std::cout<<"read consists of"<<colchunknum*rowchunknum<<"chunks"<<std::endl;
  // size_t slab_bytes= rchunksize*cchunksize*sizeof(float);
  // std::cout<<"Slab is "<<slab_bytes<<"bytes"<<std::endl;
  // plist.setCache(nmdc,cache_elem,slab_bytes,1);


  matrix_dims[0]=rchunksize;
  matrix_dims[1]=cchunksize;
  datadims[0]=cchunksize;
  datadims[1]=rchunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[2];
  offseta[0]=col_offset;
  offseta[1]=row_offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
//  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
//  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  arma::fmat retmat(matrix_dims[0],matrix_dims[1]);
  DataSpace memspace(2,matrix_dims);
//  std::cout<<"Reading data"<<std::endl;
  dataset->read(retmat.memptr(),dt,memspace,fspace);
//  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  return(retmat);
}


//[[Rcpp::export]]
arma::mat read_2ddmat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const size_t row_offset,const size_t col_offset,const size_t row_chunksize,const size_t col_chunksize){
  //Try breaking up reads in to chunks
  H5FilePtr file=open_file(h5file);
  H5GroupPtr group= open_group(file,groupname);
  H5DataSetPtr dataset = open_blosc_dataset(group,dataname);
  DSetCreatPropList cparms= dataset->getCreatePlist();
  int nmdc;
  size_t cache_elem=0;
  size_t cache_bytes=0;
  //  unsigned int majv,minv,relv;
  //  H5Library::getLibVersion(majv,minv,relv);
  //  std::cout<<"HDF5 Version:"<<majv<<"."<<minv<<"."<<relv<<std::endl;
  //  double policy;
  // FileAccPropList plist= file->getAccessPlist();
  // plist.getCache(nmdc,cache_elem,cache_bytes,policy);
  // std::cout<<"Metadata cache has  "<<nmdc<<std::endl;
  // std::cout<<"Elements in raw data chunk cache  is"<<cache_elem<<std::endl;
  // std::cout<<"Bytes in raw data chunk cache  is"<<cache_bytes<<std::endl;
  DataType dt= dataset->getDataType();
  hsize_t datadims[]={0,0};
  DataSpace fspace=dataset->getSpace();
  fspace.getSimpleExtentDims(datadims,NULL);
//  std::cout<<"Full data is of dimensions"<<datadims[0]<<"x"<<datadims[1]<<std::endl;
  hsize_t matrix_dims[2];
  size_t rchunksize=row_chunksize;
  size_t cchunksize=col_chunksize;
  if(row_offset+rchunksize>datadims[1]){
    rchunksize=datadims[1]-row_offset;
  }
  if(col_offset+cchunksize>datadims[0]){
    cchunksize=datadims[0]-col_offset;
  }
  hsize_t filechunksize[2];
  cparms.getChunk(2,filechunksize);
  size_t colchunknum=std::ceil((double)cchunksize/(double)filechunksize[0]);
  size_t rowchunknum=std::ceil((double)rchunksize/(double)filechunksize[1]);
  // std::cout<<"read consists of"<<colchunknum*rowchunknum<<"chunks"<<std::endl;
  // size_t slab_bytes= rchunksize*cchunksize*sizeof(float);
  // std::cout<<"Slab is "<<slab_bytes<<"bytes"<<std::endl;
  // plist.setCache(nmdc,cache_elem,slab_bytes,1);


  matrix_dims[0]=rchunksize;
  matrix_dims[1]=cchunksize;
  datadims[0]=cchunksize;
  datadims[1]=rchunksize;

  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[2];
  offseta[0]=col_offset;
  offseta[1]=row_offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadims,offseta);
//  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
//  std::cout<<"Matrix starts at"<<row_offset<<"x"<<col_offset<<std::endl;
  arma::mat retmat(matrix_dims[0],matrix_dims[1]);
  DataSpace memspace(2,matrix_dims);
//  std::cout<<"Reading data"<<std::endl;
  dataset->read(retmat.memptr(),dt,memspace,fspace);
//  std::cout<<"Read complete!"<<std::endl;
  dt.close();
  memspace.close();
  fspace.close();
  group->close();
  dataset->close();
  file->close();
  return(retmat);
}

//[[Rcpp::export]]
arma::mat dmeansd_h5(const std::string h5file,const std::string groupname,const std::string dataname,int dim=0){
  size_t rownum= get_rownum_h5(h5file,groupname,dataname);
  size_t colnum= get_colnum_h5(h5file,groupname,dataname);
  arma::mat data= read_2ddmat_h5(h5file,groupname,dataname,0,0,rownum,colnum);
  arma::vec meanvec= arma::conv_to<arma::vec>::from(arma::mean(data,dim));
  arma::vec sdvec = arma::conv_to<arma::vec>::from(arma::stddev(data,dim));
  return(arma::join_horiz(meanvec,sdvec));
}






size_t write_covmat_h5(const std::string h5file, const std::string groupname, const std::string dataname, const size_t dimension, arma::fmat &data,const arma::uword rowoffset,const arma::uword coloffset,const size_t rchunksize, const size_t cchunksize){

//  Rcpp::Rcout<<"Initializing write for file:"<<h5file<<std::endl;
//  Rcpp::Rcout<<"Data is of dimensions:"<<data.n_rows<<"x"<<data.n_cols<<std::endl;
  hsize_t chunkstart=0;
  //  size_t rchunksize = data.n_cols/5;
  //  size_t cchunksize = data.n_rows/5;

  //  std::vector<float> tdat= arma::conv_to<std::vector<float>>::from(arma::vectorise(data));
  // if(rchunksize!=data.n_cols){
  //   Rcpp::Rcerr<<"chunksize("<<rchunksize<<") not equal to column number("<<data.n_cols<<")!"<<std::endl;
  //   Rcpp::stop("error in write_covmat_h5");
  // }
  size_t urchunksize=rchunksize;
  if(rchunksize==dimension){
    urchunksize=rchunksize/2;
  }
  size_t ucchunksize=cchunksize;
  if(cchunksize==dimension){
    ucchunksize=cchunksize/2;
  }

  hsize_t mchunknum= std::ceil((double)dimension/(double)rchunksize);
  mchunknum=mchunknum*10;
  std::vector<hsize_t> cumdim{dimension,dimension};
  std::vector<hsize_t> maxdim{dimension,dimension};
  std::vector<hsize_t> chunkdim{ucchunksize,urchunksize};
  std::vector<hsize_t> refdim{1};
  std::vector<hsize_t> mrefdim{mchunknum};
  hsize_t adim[2];//dimension of each array element (ncols)
  FloatType ftypew(PredType::NATIVE_FLOAT);
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  // // // std::cout<<"Creating/opening Reference dataset"<<std::endl;
  // H5DataSetPtr refdataset = create_or_open_ref_dataset(group,dataname,refdim,mrefdim);
//    std::cout<<"Creating/opening Real dataset"<<std::endl;
    size_t deflate_level=4;
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,ftypew,cumdim,maxdim,chunkdim,deflate_level);

  DataSpace* fdataspace;
  try{
    fdataspace= new DataSpace(dataset->getSpace());
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }
  // DataSpace* reffdataspace;
  // try{
  //   reffdataspace= new DataSpace(refdataset->getSpace());
  // }catch(DataSpaceIException error){
  //   error.printError();
  //   Rcpp::stop("Error creating reference dataspace ");
  // }

  // hsize_t refdatadim[]={0};
  // reffdataspace->getSimpleExtentDims(refdatadim,NULL);
  // refdatadim[0]=refdatadim[0]+1;
  // refdataset->extend(refdatadim);
  // hsize_t refoffset[1];
  // refoffset[0]=refdatadim[0]-1;

  hsize_t datadim[2];
  fdataspace->getSimpleExtentDims(datadim,NULL);
//  std::cout<<"old data dim"<<datadim[0]<<"x"<<datadim[1]<<std::endl;
//  std::cout<<"data chunk size is "<<urchunksize<<"x"<<ucchunksize<<std::endl;
  // chunkstart = datadim[0];
  // size_t chunka[]={chunksize,chunksize};
  hsize_t memdim[]={data.n_cols,data.n_rows};
  // fdataspace->close();
  // delete fdataspace;
  // fdataspace = new DataSpace(dataset->getSpace());
  // fdataspace->getSimpleExtentDims(datadim,NULL);
  // std::cout<<"new data dim"<<datadim[0]<<std::endl;
  DataSpace *mspace;
  try{
    mspace= new DataSpace(2,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }
  hsize_t odim[]={coloffset,rowoffset};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1,1};
  hsize_t blocka[]={1,1};


  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim,odim,stridea,blocka);
  hsize_t block_start[]={0,0};
  hsize_t block_stop[]={0,0};


//  std::cout<<"Starting to write data"<<std::endl;
  try{
    dataset->write(data.memptr(),PredType::NATIVE_FLOAT,*mspace,*fdataspace);
  }  catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("Error writing file");
  }
  std::cout<<"Data sucessfully written"<<std::endl;
  try{
    file->flush(H5F_SCOPE_GLOBAL);
  }catch(FileIException error)
  {
    error.printError();
    Rcpp::stop("Error flushing file");
  }
  std::cout<<"File flushed"<<std::endl;
  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  fdataspace->close();
  mspace->close();
  return(rchunksize);
}


size_t write_blosc_covmat_h5(const std::string h5file, const std::string groupname, const std::string dataname, const size_t dimension, arma::fmat &data,const arma::uword rowoffset,const arma::uword coloffset,const int deflate_level){

  Rcpp::Rcout<<"Initializing write for file:"<<h5file<<std::endl;
  Rcpp::Rcout<<"Data is of dimensions:"<<data.n_rows<<"x"<<data.n_cols<<std::endl;

  hsize_t chunkstart=0;
  size_t rchunksize = data.n_cols/2;
  size_t cchunksize = data.n_rows;
  std::vector<float> tdat= arma::conv_to<std::vector<float>>::from(arma::vectorise(data));
  // if(rchunksize!=data.n_cols){
  //   Rcpp::Rcerr<<"chunksize("<<rchunksize<<") not equal to column number("<<data.n_cols<<")!"<<std::endl;
  //   Rcpp::stop("error in write_covmat_h5");
  // }
  size_t urchunksize=rchunksize;
  if(rchunksize==dimension){
    urchunksize=rchunksize;
  }
  size_t ucchunksize=cchunksize;
  if(cchunksize==dimension){
    ucchunksize=cchunksize;
  }
  std::vector<hsize_t> cumdim{dimension,dimension};
  std::vector<hsize_t> maxdim{dimension,dimension};
  std::vector<hsize_t> chunkdim{ucchunksize,urchunksize};

  FloatType ftypew(PredType::NATIVE_FLOAT);
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,ftypew,cumdim,maxdim,chunkdim,deflate_level);



  DataSpace* fdataspace;
  try{
    fdataspace= new DataSpace(dataset->getSpace());
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }
  hsize_t datadim[2];
  fdataspace->getSimpleExtentDims(datadim,NULL);
  std::cout<<"old data dim"<<datadim[0]<<"x"<<datadim[1]<<std::endl;
  std::cout<<"data chunk size is "<<urchunksize<<"x"<<ucchunksize<<std::endl;
  // chunkstart = datadim[0];
  // size_t chunka[]={chunksize,chunksize};
  // datadim[0]=datadim[0]+chunksize[0];
  // datadim[1]=datadim[1]+chunksize[1];
  // dataset->extend(datadim);
  hsize_t memdim[]={data.n_cols,data.n_rows};
  // fdataspace->close();
  // delete fdataspace;
  // fdataspace = new DataSpace(dataset->getSpace());
  // fdataspace->getSimpleExtentDims(datadim,NULL);
  // std::cout<<"new data dim"<<datadim[0]<<std::endl;
  DataSpace *mspace;
  try{
    mspace= new DataSpace(2,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  }catch(DataSpaceIException error){
    error.printError();
    Rcpp::stop("Error creating memory dataspace ");
  }
  hsize_t odim[]={coloffset,rowoffset};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1,1};
  hsize_t blocka[]={1,1};
  std::cout<<"Getting (disk) Data dimensions"<<std::endl;
  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim, odim,stridea,blocka);
  hsize_t block_start[]={0,0};
  hsize_t block_stop[]={0,0};

  fdataspace->getSelectBounds(block_start,block_stop);
  std::cout<<"Writing will be from row"<<block_start[0]<<"to row"<<block_stop[0]<<"size="<<block_stop[0]-block_start[0]+1<<std::endl;
  std::cout<<"Writing will be from row"<<block_start[1]<<"to row"<<block_stop[1]<<"size="<<block_stop[1]-block_start[1]+1<<std::endl;
  std::cout<<"Starting to write data"<<std::endl;
  try{
    dataset->write(&tdat[0],PredType::NATIVE_FLOAT,*mspace,*fdataspace);
  }  catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("Error writing file");
  }
  std::cout<<"Data sucessfully written"<<std::endl;
  try{
    file->flush(H5F_SCOPE_GLOBAL);
  }catch(FileIException error)
  {
    error.printError();
    Rcpp::stop("Error flushing file");
  }
  std::cout<<"File flushed"<<std::endl;
  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  delete fdataspace;
  delete mspace;
  return(rchunksize);
}



size_t write_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, const hsize_t Nind, arma::mat &data,const  int deflate_level){

  Rcpp::Rcout<<"Initializing write for file:"<<h5file<<std::endl;
  Rcpp::Rcout<<"Initiating write_mat on:"<<dataname<<std::endl;
  hsize_t chunkstart=0;
  size_t chunksize = data.size()/Nind;
  if(chunksize!=data.n_cols){
    Rcpp::Rcerr<<"chunksize not equal to column number!"<<std::endl;
    Rcpp::stop("error in write_mat_h5");
  }
  size_t uchunksize=chunksize;
  if(chunksize==Nsnps){
    uchunksize=chunksize/2;
  }
  std::vector<hsize_t> cumdim{0};
  std::vector<hsize_t> maxdim{Nsnps};
  std::vector<hsize_t> chunkdim{uchunksize};
  hsize_t adim[1];//dimension of each array element (ncols)
  adim[0]=Nind; //Assign array size dimension
  ArrayType ahtypew(PredType::NATIVE_DOUBLE,1,adim); //Create array size type
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,ahtypew,cumdim,maxdim,chunkdim,deflate_level);
  DataSpace* fdataspace = new DataSpace(dataset->getSpace());
  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);
  std::cout<<"old data dim"<<datadim[0]<<std::endl;
  chunkstart = datadim[0];
  datadim[0]=datadim[0]+chunksize;
  dataset->extend(datadim);
  hsize_t memdim[]={chunksize};
  fdataspace->close();
  delete fdataspace;
  fdataspace = new DataSpace(dataset->getSpace());
  fdataspace->getSimpleExtentDims(datadim,NULL);
  std::cout<<"new data dim"<<datadim[0]<<std::endl;
  DataSpace *mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  hsize_t odim[]={chunkstart};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};
  //  std::cout<<"Getting Data dimensions"<<std::endl;
  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim, odim,stridea,blocka);
  // std::cout<<"Starting to write data"<<std::endl;
  dataset->write(data.memptr(),ahtypew,*mspace,*fdataspace);
  file->flush(H5F_SCOPE_GLOBAL);
  // std::cout<<"Data sucessfully written"<<std::endl;
  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  delete fdataspace;
  delete mspace;
  return(chunksize);
}


size_t write_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, const hsize_t Nind, arma::fmat &data,const  int deflate_level){

  Rcpp::Rcout<<"Initiating write_mat on:"<<dataname<<std::endl;
  hsize_t chunkstart=0;
  size_t chunksize = data.size()/Nind;
  size_t uchunksize=chunksize;
  if(chunksize==Nsnps){
    uchunksize=100;
  }
  std::vector<hsize_t> cumdim{0};
  std::vector<hsize_t> maxdim{Nsnps};
  std::vector<hsize_t> chunkdim{uchunksize};
  hsize_t adim[1];//dimension of each array element (ncols)
  adim[0]=Nind; //Assign array size dimension
  ArrayType ahtypew(PredType::NATIVE_FLOAT,1,adim); //Create array size type
  Rcpp::Rcout<<"Opening file for writing: "<<h5file<<std::endl;
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,ahtypew,cumdim,maxdim,chunkdim,deflate_level);
  DataSpace* fdataspace = new DataSpace(dataset->getSpace());
  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);
  std::cout<<"old data dim"<<datadim[0]<<std::endl;
  chunkstart = datadim[0];
  datadim[0]=datadim[0]+chunksize;
  try{
  dataset->extend(datadim);
  }    catch( DataSetIException error )
  {
    error.printError();
    Rcpp::stop("Error extending dataset");
  }
  hsize_t memdim[]={chunksize};
  fdataspace->close();
  delete fdataspace;
  fdataspace = new DataSpace(dataset->getSpace());
  fdataspace->getSimpleExtentDims(datadim,NULL);
  std::cout<<"new data dim"<<datadim[0]<<std::endl;
  DataSpace *mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  hsize_t odim[]={chunkstart};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};
  //  std::cout<<"Getting Data dimensions"<<std::endl;
  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim, odim,stridea,blocka);
  std::cout<<"Starting to write data"<<std::endl;
  dataset->write(data.memptr(),ahtypew,*mspace,*fdataspace);
  file->flush(H5F_SCOPE_GLOBAL);
  std::cout<<"Data sucessfully written"<<std::endl;
  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  delete fdataspace;
  delete mspace;
  return(chunksize);
}




size_t write_int_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, arma::uvec &data,const  int deflate_level){


  typedef std::vector<int> istd;
  std::vector<int>tdata = arma::conv_to<istd>::from(data);
  hsize_t chunkstart=0;
  size_t chunksize = tdata.size();
  std::vector<hsize_t> cumdim{0};
  std::vector<hsize_t> maxdim{Nsnps};
  std::vector<hsize_t> chunkdim{chunksize};
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,PredType::NATIVE_INT,cumdim,maxdim,chunkdim,deflate_level);

  DataSpace* fdataspace = new DataSpace(dataset->getSpace());
  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);
  chunkstart = datadim[0];
  datadim[0]=datadim[0]+chunksize;
  dataset->extend(datadim);
  hsize_t memdim[]={chunksize};
  fdataspace->close();
  delete fdataspace;
  fdataspace = new DataSpace(dataset->getSpace());
  fdataspace->getSimpleExtentDims(datadim,NULL);

  DataSpace *mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  hsize_t odim[]={chunkstart};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};
  //  std::cout<<"Getting Data dimensions"<<std::endl;
  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim, odim,stridea,blocka);
  std::cout<<"Starting to write data"<<std::endl;
  dataset->write(&tdata[0],PredType::NATIVE_INT,*mspace,*fdataspace);
  file->flush(H5F_SCOPE_GLOBAL);
  dataset->close();
  mspace->close();
  fdataspace->close();
  group->close();
  file->close();
  return(0);
}



size_t write_float_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, arma::fvec &data,const  int deflate_level){

  typedef std::vector<float> istd;
  std::vector<float>tdata = arma::conv_to<istd>::from(data);
  hsize_t chunkstart=0;
  size_t chunksize = tdata.size();
  std::vector<hsize_t> cumdim{0};
  std::vector<hsize_t> maxdim{Nsnps};
  std::vector<hsize_t> chunkdim{chunksize};
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,PredType::NATIVE_FLOAT,cumdim,maxdim,chunkdim,deflate_level);

  DataSpace* fdataspace = new DataSpace(dataset->getSpace());
  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);
  chunkstart = datadim[0];
  datadim[0]=datadim[0]+chunksize;
  dataset->extend(datadim);
  hsize_t memdim[]={chunksize};
  fdataspace->close();
  delete fdataspace;
  fdataspace = new DataSpace(dataset->getSpace());
  fdataspace->getSimpleExtentDims(datadim,NULL);

  DataSpace *mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  hsize_t odim[]={chunkstart};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};
  //  std::cout<<"Getting Data dimensions"<<std::endl;
  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim, odim,stridea,blocka);
  std::cout<<"Starting to write data"<<std::endl;
  dataset->write(&tdata[0],PredType::NATIVE_FLOAT,*mspace,*fdataspace);
  file->flush(H5F_SCOPE_GLOBAL);
  dataset->close();
  mspace->close();
  fdataspace->close();
  group->close();
  file->close();
  return(0);
}


size_t write_uint_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, arma::uvec &data, const int deflate_level){

  typedef std::vector<unsigned int> uistd;
  std::vector<unsigned int>tdata = arma::conv_to<uistd>::from(data);
  hsize_t chunkstart=0;
  size_t chunksize = tdata.size();
  std::vector<hsize_t> cumdim{0};
  std::vector<hsize_t> maxdim{Nsnps};
  std::vector<hsize_t> chunkdim{chunksize};
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,PredType::NATIVE_UINT,cumdim,maxdim,chunkdim,deflate_level);

  DataSpace* fdataspace = new DataSpace(dataset->getSpace());
  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);
  chunkstart = datadim[0];
  datadim[0]=datadim[0]+chunksize;
  dataset->extend(datadim);
  hsize_t memdim[]={chunksize};
  fdataspace->close();
  delete fdataspace;
  fdataspace = new DataSpace(dataset->getSpace());
  fdataspace->getSimpleExtentDims(datadim,NULL);

  DataSpace *mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  hsize_t odim[]={chunkstart};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};
  //  std::cout<<"Getting Data dimensions"<<std::endl;
  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim, odim,stridea,blocka);
  std::cout<<"Starting to write data"<<std::endl;
  dataset->write(&tdata[0],PredType::NATIVE_UINT,*mspace,*fdataspace);
  file->flush(H5F_SCOPE_GLOBAL);
  dataset->close();
  mspace->close();
  fdataspace->close();
  group->close();
  file->close();
  return(0);
}




size_t write_umat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, const hsize_t Nind, arma::umat &data,const  int deflate_level){

  typedef std::vector<unsigned int> uistd;
  std::vector<unsigned int>tdata = arma::conv_to<uistd>::from(vectorise(data));
  hsize_t chunkstart=0;
  size_t chunksize = data.size()/Nind;
  std::vector<hsize_t> cumdim{0};
  std::vector<hsize_t> maxdim{Nsnps};
  std::vector<hsize_t> chunkdim{chunksize};
  hsize_t adim[1];//dimension of each array element (ncols)
  adim[0]=Nind; //Assign array size dimension
  ArrayType ahtypew(PredType::NATIVE_UINT,1,adim); //Create array size type
  H5FilePtr file =create_or_open_file(h5file);
  H5GroupPtr group=create_or_open_group(file,groupname);
  H5DataSetPtr dataset = create_or_open_dataset(group,dataname,ahtypew,cumdim,maxdim,chunkdim,deflate_level);
  DataSpace* fdataspace = new DataSpace(dataset->getSpace());
  hsize_t datadim[1];
  fdataspace->getSimpleExtentDims(datadim,NULL);

  chunkstart = datadim[0];
  datadim[0]=datadim[0]+chunksize;
  dataset->extend(datadim);
  hsize_t memdim[]={chunksize};
  fdataspace->close();
  delete fdataspace;
  fdataspace = new DataSpace(dataset->getSpace());
  fdataspace->getSimpleExtentDims(datadim,NULL);

  DataSpace *mspace= new DataSpace(1,memdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)
  hsize_t odim[]={chunkstart};//dimension of each offset (current_chunk*chunksize)
  hsize_t stridea[]={1};
  hsize_t blocka[]={1};
  //  std::cout<<"Getting Data dimensions"<<std::endl;
  fdataspace->selectHyperslab( H5S_SELECT_SET, memdim, odim,stridea,blocka);
  std::cout<<"Starting to write data"<<std::endl;
  dataset->write(&tdata[0],ahtypew,*mspace,*fdataspace);
  file->flush(H5F_SCOPE_GLOBAL);

  dataset->close();
  fdataspace->close();
  mspace->close();
  group->close();
  file->close();
  delete fdataspace;
  delete mspace;

  return(0);
}



//[[Rcpp::export]]
int write_dmatrix_h5(Rcpp::String h5file,Rcpp::String groupname, Rcpp::String dataname, Rcpp::IntegerVector Nsnps, Rcpp::IntegerVector Nind, Rcpp::NumericMatrix data,const  int deflate_level){
  std::string th5 =h5file;
  std::string gn = groupname;
  std::string dn = dataname;
  hsize_t nsnp = Nsnps[0];
  hsize_t nid= Nind[0];
  arma::fmat tdat = Rcpp::as<arma::fmat>(data);
  Rcpp::Rcout<<" Writing "<<th5<<"Group name: "<<gn<<" data name:"<<dn<<" nsnp(cols): "<<nsnp<<" nid(rows)"<<nid<<" mat[0,0]:"<<tdat(0,0)<<std::endl;
  int ret = write_mat_h5(th5,gn,dn,nsnp,nid,tdat,deflate_level);
  return(ret);
}




//[[Rcpp::export]]
int write_Rint_h5(const std::string h5file, const std::string groupname, const std::string dataname, Rcpp::IntegerVector &data,const int deflate_level){
  arma::uvec tdat = Rcpp::as<arma::uvec>(data);
  Rcpp::Rcout<<" Writing "<<h5file<<"Group name: "<<groupname<<" data name:"<<dataname<<std::endl;
  int ret = write_int_h5(h5file,groupname,dataname,H5S_UNLIMITED,tdat,deflate_level);
  return(ret);
}

//[[Rcpp::export]]
int write_Rnumeric_h5(const std::string h5file, const std::string groupname, const std::string dataname, Rcpp::NumericVector &data,const int deflate_level){
  arma::fvec tdat = Rcpp::as<arma::fvec>(data);
  Rcpp::Rcout<<" Writing "<<h5file<<"Group name: "<<groupname<<" data name:"<<dataname<<std::endl;
  int ret = write_float_h5(h5file,groupname,dataname,H5S_UNLIMITED,tdat,deflate_level);
  return(ret);
}







