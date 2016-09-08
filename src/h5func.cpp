#include "RcppArmadillo.h"
#include <vector>
#include <memory>
#include <H5Cpp.h>
#include "h5func.hpp"


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
    file = new H5::H5File(fname.c_str(), H5F_ACC_TRUNC);
  }else{
    file = new H5::H5File(fname.c_str(), H5F_ACC_RDWR);
  }
  return H5FilePtr(file);
}


H5GroupPtr create_or_open_group(H5FilePtr &file,const std::string &groupname)
{
  Group* group;
  Group *rg = new Group(file->openGroup("/"));
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
    group = new Group(file->openGroup(groupname));
  }else{
    group = new Group(file->createGroup("/"+groupname));
  }
  return H5GroupPtr(group);
}

H5DataSetPtr create_or_open_dataset(H5GroupPtr &group,const std::string &dataname, const DataType &data_type,std::vector<hsize_t> &cdatadim,std::vector<hsize_t> &mdatadim,std::vector<hsize_t> &chunkdim,const unsigned int deflate_level)
{
  H5::Exception::dontPrint();
  DataSet* dataset;
  DataSpace* fdataspace;
  hsize_t objc= group->getNumObjs();
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
    dataset = new DataSet(group->openDataSet(dataname));
  }else{
    hsize_t *cumdima = new hsize_t[cdatadim.size()];
    hsize_t *chunkdima = new hsize_t[chunkdim.size()];
    hsize_t*mdima = new hsize_t[mdatadim.size()];

    for(int i=0; i<cdatadim.size();i++){
      cumdima[i]=cdatadim[i];
      mdima[i]=mdatadim[i];
      chunkdima[i]=chunkdim[i];
    }
    fdataspace= new DataSpace(cdatadim.size(),cumdima,mdima); //Create dataspace for dataset (on disk)
    DSetCreatPropList cparms; //Create chunksize file parameters
    cparms.setChunk(chunkdim.size(),chunkdima); //Set chunksize
    cparms.setDeflate(deflate_level);
    dataset = new DataSet(group->createDataSet(dataname,data_type,*fdataspace,cparms));

  }
  return H5DataSetPtr(dataset);
}


ArrayType read_arraytype(const DataSet *dataset, const PredType pt){
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
size_t get_rownum_h5(const std::string hap_h5file,const std::string groupname, const std::string dataname){

  H5File* file;
  DataSet* dataset;
  try{
    file= new H5File( hap_h5file.c_str(), H5F_ACC_RDONLY);
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
  }
  DataSpace fspace =dataset->getSpace();

  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  fspace.close();
  dataset->close();
  file->close();
  return(datadim[0]);
}



//[[Rcpp::export]]
arma::mat read_dmat_h5(const std::string hap_h5file,const std::string groupname, const std::string dataname, size_t offset, size_t chunksize){

  H5File* file;
  DataSet* dataset;
  try{
    file= new H5File( hap_h5file.c_str(), H5F_ACC_RDONLY);
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
  }

  //  std::cout<<"opened dataset"<<std::endl;
  ArrayType mem_arraytype=read_arraytype(dataset,PredType::NATIVE_DOUBLE);


  //  std::cout<<"Getting Array dimensions"<<std::endl;
  //  std::cout<<"Getting Dataspace"<<std::endl;
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


  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t offseta[1];
  offseta[0]=offset;
  fspace.selectHyperslab(H5S_SELECT_SET,datadim,offseta);
  std::cout<<"Allocating matrix of size:"<<matrix_dims[0]<<"x"<<matrix_dims[1]<<std::endl;
  arma::mat retmat(matrix_dims[0],matrix_dims[1]);
  DataSpace memspace(1,datadim);
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read(retmat.memptr(),mem_arraytype,memspace,fspace);
  mem_arraytype.close();
  memspace.close();
  fspace.close();
  haplogroup.close();
  delete dataset;
  delete file;
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




std::vector<int> read_int_h5(const std::string h5file, const std::string groupname, const std::string dataname){
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
  }
  DataType dt = dataset->getDataType();
  DataSpace fspace =dataset->getSpace();
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  hsize_t readrows=datadim[0];
  DataSpace memspace(1,datadim);
  std::vector<int> retdat(datadim[0]);

  std::cout<<"Reading in int dataset dataname"<<std::endl;

  dataset->read(&retdat[0],dt,memspace,fspace);
  fspace.close();
  memspace.close();
  dataset->close();
  file->close();

  std::cout<<"Last element of data is :"<<retdat[retdat.size()-1]<<std::endl;
  return(retdat);

}


std::vector<unsigned int> read_uint_h5(const std::string h5file, const std::string groupname, const std::string dataname){
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
  }
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




size_t write_mat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, const hsize_t Nind, arma::mat &data,const unsigned int deflate_level){

  Rcpp::Rcout<<"Initializing write for file:"<<h5file<<std::endl;
  Rcpp::Rcout<<"Initiating write_mat on:"<<dataname<<std::endl;
  hsize_t chunkstart=0;
  size_t chunksize = data.size()/Nind;
  if(chunksize!=data.n_cols){
    Rcpp::Rcerr<<"chunksize not equal to column number!"<<std::endl;
    Rcpp::stop("error int write_mat_h5");
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
  std::cout<<"Starting to write data"<<std::endl;
  dataset->write(data.memptr(),ahtypew,*mspace,*fdataspace);
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




size_t write_int_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, arma::uvec &data,const unsigned int deflate_level){


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
  dataset->close();
  mspace->close();
  fdataspace->close();
  group->close();
  file->close();
  return(0);
}



size_t write_float_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, arma::fvec &data,const unsigned int deflate_level){

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
  dataset->close();
  mspace->close();
  fdataspace->close();
  group->close();
  file->close();
  return(0);
}


size_t write_uint_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, arma::uvec &data, const unsigned int deflate_level){

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
  dataset->close();
  mspace->close();
  fdataspace->close();
  group->close();
  file->close();
  return(0);
}




size_t write_umat_h5(const std::string h5file, const std::string groupname, const std::string dataname,const hsize_t Nsnps, const hsize_t Nind, arma::umat &data,const unsigned int deflate_level){

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
int write_dmatrix_h5(Rcpp::String h5file,Rcpp::String groupname, Rcpp::String dataname, Rcpp::IntegerVector Nsnps, Rcpp::IntegerVector Nind, Rcpp::NumericMatrix data,const unsigned int deflate_level){
  std::string th5 =h5file;
  std::string gn = groupname;
  std::string dn = dataname;
  hsize_t nsnp = Nsnps[0];
  hsize_t nid= Nind[0];
  arma::mat tdat = Rcpp::as<arma::mat>(data);
  Rcpp::Rcout<<" Writing "<<th5<<"Group name: "<<gn<<" data name:"<<dn<<" nsnp(cols): "<<nsnp<<" nid(rows)"<<nid<<" mat[0,0]:"<<tdat(0,0)<<std::endl;
  int ret = write_mat_h5(th5,gn,dn,nsnp,nid,tdat,deflate_level);
  return(ret);
}


//[[Rcpp::export]]
int write_Rint_h5(const std::string h5file, const std::string groupname, const std::string dataname, Rcpp::IntegerVector &data,const unsigned int deflate_level){
  arma::uvec tdat = Rcpp::as<arma::uvec>(data);
  Rcpp::Rcout<<" Writing "<<h5file<<"Group name: "<<groupname<<" data name:"<<dataname<<std::endl;
  int ret = write_int_h5(h5file,groupname,dataname,H5S_UNLIMITED,tdat,deflate_level);
  return(ret);
}

//[[Rcpp::export]]
int write_Rnumeric_h5(const std::string h5file, const std::string groupname, const std::string dataname, Rcpp::NumericVector &data,const unsigned int deflate_level){
  arma::fvec tdat = Rcpp::as<arma::fvec>(data);
  Rcpp::Rcout<<" Writing "<<h5file<<"Group name: "<<groupname<<" data name:"<<dataname<<std::endl;
  int ret = write_float_h5(h5file,groupname,dataname,H5S_UNLIMITED,tdat,deflate_level);
  return(ret);
}







