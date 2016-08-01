#include "RcppArmadillo.h"
#include <H5Cpp.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <zlib.h>



// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//
using namespace H5;

// [[Rcpp::export]]
void write_haplotype_h5(const std::string hap_gzfile,const std::string hap_h5file,const size_t nrows,const size_t ncols,const size_t chunksize){

  gzFile gz_in=gzopen(hap_gzfile.c_str(),"rb");
  size_t rc=0;
  int* haplotypes = new int[chunksize*ncols];
  H5File* file;
  try{
    H5File* file = new H5File( hap_h5file.c_str(), H5F_ACC_TRUNC);
    hsize_t dim[1]; //dimensions of overall dataspace (nrows)
    dim[0]=nrows; //Assign total dataspace dimension
    DataSpace space(1,dim); //Create dataspace for dataset (on disk)

    hsize_t adim[1];//dimension of each array element (ncols)
    adim[0]=ncols; //Assign array size dimension
    std::cout<<"Size of haplotype element is "<<sizeof(int)<<std::endl;
    std::cout<<"Size of *haplotype element is "<<sizeof(int*)<<std::endl;
    ArrayType ahtypem(PredType::NATIVE_INT32,1,adim);
    ArrayType ahtypew(PredType::NATIVE_INT32,1,adim); //Create array size type

    hsize_t cdim[1];//dimension of each chunk (chunksize)
    cdim[0]=chunksize; //Assign chunksize dimension
    DSetCreatPropList cparms; //Create chunksize file parameters
    cparms.setChunk(1,cdim); //Set chunksize
    cparms.setDeflate(2);

    DataSpace mspace(1,cdim); //Size of first dataset (in memory, can be bigger or smaller than size on disk, depending on how much you're writing)

    hsize_t odim[1];//dimension of each offset (current_chunk*chunksize)
    odim[0]=0; //Assign first offset to zero
    hsize_t stridea[]={1};
    hsize_t blocka[]={1};
    Group haplogroup =file->createGroup("/Haplotype");
    DataSet dataset = file->createDataSet("/Haplotype/haplotype",ahtypem,space,cparms);

    DataSpace fspace = dataset.getSpace();
    size_t len= chunksize*ncols*2;
    char *buffer= new char[len];
    size_t bytes_read=0;

    for(;;){
      bytes_read=gzread(gz_in,buffer,len);
      if(bytes_read==0){
        break;
      }
      size_t chunkrs=bytes_read/(ncols*2);
//      std::cout<<"Read "<<chunkrs<<"rows"<<std::endl;
//      std::cout<<"Row "<<odim[0]<<" of "<<nrows<<std::endl;
      size_t bc=0;
      size_t i=0;
      for(size_t row=0; row<chunkrs; row++){
        if(row%10000==0){
          std::cout<<"row :"<<odim[0]+row<<" hap[row][0]:"<<haplotypes[i]<<std::endl;
        }
          //          std::cout<<"buffer[bc]:'"<<buffer[bc]<<"'="<<(int)(buffer[bc]-'0')<<" col="<<col<<" row="<<row<<std::endl;
          for(size_t col=0; col<ncols; col++){
            haplotypes[i]=(int)(buffer[bc]-'0');
            i++;
            bc+=2;
          }
      }
      //      std::cout<<"Finished processing haplotype array"<<std::endl;
      //      std::cout<<"Selected Hyperslab"<<std::endl;
      /*
      Args(
        H5S_SELECT_SET,
      number of "blocks"
      (chunksize),
      start of slab (offset),
      stride(1),
      size of block(1)
      Hyperslabs can be selected
      from memory as well as
      from file spaces, in this
      case we're doing it in the
      file because the data is
      larger than what we want
      in memory
      */
      std::cout<<"Selecting hyperslab of size: "<<cdim[0]<<" at offset of: "<<odim[0]<<" stride size: "<<stridea[0]<<"block size: "<<blocka[0]<<std::endl;
      fspace.selectHyperslab( H5S_SELECT_SET, cdim, odim,stridea,blocka);
      if(chunkrs!=chunksize){
        std::cout<<"Resizing "<<std::endl;
        hsize_t ncdim[1];
        ncdim[0]=chunkrs;
        DataSpace nmspace(1,ncdim);
        std::cout<<"Starting to write data"<<std::endl;
        fspace.selectHyperslab( H5S_SELECT_SET, ncdim, odim,stridea,blocka);
        dataset.write(haplotypes,ahtypew,nmspace,fspace);
        nmspace.close();
        std::cout<<"Data written"<<std::endl;
      }else{
        std::cout<<"Starting to write data"<<std::endl;
        dataset.write(haplotypes,ahtypew,mspace,fspace);
        std::cout<<"Data written"<<std::endl;
      }
      /*
       Data buffer(data),
       Memory datatype,
       memory dataspace(1xchunksize),
       file_space,Transfer property list
       */
      odim[0]+=chunkrs;
    }
    fspace.close();
    dataset.close();
    haplogroup.close();
    mspace.close();
    cparms.close();
    ahtypem.close();
    ahtypew.close();
    file->close();

  }
  catch( FileIException error )
  {
    error.printError();

  }
  // catch failure caused by the DataSet operations
  catch( DataSetIException error )
  {
    error.printError();

  }
  // catch failure caused by the DataSpace operations
  catch( DataSpaceIException error )
  {
    error.printError();

  }
  // catch failure caused by the DataSpace operations
  catch( DataTypeIException error )
  {
    error.printError();

  }
  std::cout<<"Data written (done)"<<std::endl;
}

//[[Rcpp::export]]
arma::Mat<int> read_haplotype_ind_h5(const std::string hap_h5file,Rcpp::IntegerVector indexes){

  H5File* file;
  DataSet* dataset;
  hid_t dtype;
  try{
    file= new H5File( hap_h5file.c_str(), H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
  }
  Group haplogroup = file->openGroup("Haplotype");
  try{
    dataset = new DataSet(haplogroup.openDataSet("haplotype"));
  }
  catch( DataSetIException error )
  {
    error.printError();
  }

  std::cout<<"opened dataset"<<std::endl;
  hsize_t adim[1];
  DataType dt = dataset->getDataType();
  DataType dtss = dt.getSuper();
  size_t dts=dt.getSize();
  size_t bcs =dtss.getSize();
  hsize_t dims= (dts/bcs);
  adim[0]=dims;
  ArrayType mem_arraytype(PredType::NATIVE_INT32,1,adim);


  std::cout<<"Getting Array dimensions"<<std::endl;
  size_t ncol=adim[0];
  std::cout<<"Getting Dataspace"<<std::endl;
  DataSpace fspace =dataset->getSpace();

  std::cout<<"Getting Data dimensions"<<std::endl;
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t readSNPs=indexes.length();
  if(readSNPs>datadim[0]){
    Rcpp::stop("indexes longer than extent of data");
  }
  hsize_t *hind = new hsize_t[readSNPs];
  for(size_t i=0; i<readSNPs;i++){
    hind[i]=indexes[i]-1;
    if(hind[i]>=datadim[0]){
      Rcpp::stop("Attempting to read outside extent of data");
    }
  }
  hsize_t rdim[1];
  rdim[0]=readSNPs;
  DataSpace memspace(1,rdim);
  fspace.selectElements(H5S_SELECT_SET,readSNPs,hind);
  std::cout<<"Allocating temp vector of size:"<<readSNPs<<"x"<<adim[0]<<std::endl;
  int *tdat= new int[readSNPs*adim[0]];
  std::cout<<"Reading data"<<std::endl;
  dataset->read(tdat,mem_arraytype,memspace,fspace);
  arma::Mat<int> retmat(tdat,adim[0],rdim[0]);
  delete [] tdat;
  mem_arraytype.close();
  memspace.close();
  fspace.close();
  haplogroup.close();
  delete dataset;
  delete file;
  return(retmat);
}


// [[Rcpp::export]]
arma::Mat<int> read_haplotype_h5(const std::string hap_h5file,const size_t readSNPs,const size_t skipSNPs=0){

  H5File* file;
  DataSet* dataset;
  hid_t dtype;
  try{
    file= new H5File( hap_h5file.c_str(), H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
  }
  Group haplogroup = file->openGroup("Haplotype");
  try{
    dataset = new DataSet(haplogroup.openDataSet("haplotype"));
  }
  catch( DataSetIException error )
  {
    error.printError();

  }

  std::cout<<"opened dataset"<<std::endl;
  hsize_t adim[1];
  DataType dt = dataset->getDataType();
  DataType dtss = dt.getSuper();
  size_t dts=dt.getSize();

  size_t bcs =dtss.getSize();
  hsize_t dims= (dts/bcs);
  adim[0]=dims;
  ArrayType mem_arraytype(PredType::NATIVE_INT32,1,adim);


  std::cout<<"Getting Array dimensions"<<std::endl;
  size_t ncol=adim[0];
  std::cout<<"Getting Dataspace"<<std::endl;
  DataSpace fspace =dataset->getSpace();

  std::cout<<"Getting Data dimensions"<<std::endl;
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  hsize_t rdim[1];
  std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  if(skipSNPs+readSNPs<=datadim[0]){
    rdim[0]=readSNPs;
  }else{
    if(skipSNPs<datadim[0]){
      //      Rcpp::warning("skiprows+readrows > total number or rows, using a smaller number for readrows");
      std::cerr<<"skiprows+readrows > total number or rows, using a smaller number for readrows"<<std::endl;
      rdim[0]=datadim[0]-skipSNPs;
    }
    else{
      std::cerr<<"skiprows greater than number of total rows"<<std::endl;
      //      Rcpp::stop("skiprows greater than number of total rows");
    }
  }
  hsize_t offset[1];
  offset[0]=skipSNPs;
  DataSpace memspace(1,rdim);
  fspace.selectHyperslab(H5S_SELECT_SET,rdim,offset);
  std::cout<<"Allocating temp vector of size:"<<rdim[0]<<"x"<<adim[0]<<std::endl;
  int *tdat= new int[rdim[0]*adim[0]];
  std::cout<<"Reading data"<<std::endl;
  dataset->read(tdat,mem_arraytype,memspace,fspace);
  arma::Mat<int> retmat(tdat,adim[0],rdim[0]);
  delete [] tdat;
  mem_arraytype.close();
  memspace.close();
  fspace.close();
  haplogroup.close();
  delete dataset;
  delete file;
  return(retmat);
}




