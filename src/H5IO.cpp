#include "RcppArmadillo.h"
#include <H5Cpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include <zlib.h>

#define ARMA_USE_CXX11
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

//[[Rcpp::export]]
arma::mat convertTSparse(Rcpp::S4 &mat){

  Rcpp::IntegerVector dims = mat.slot("Dim");
  std::cout<<"Dim of mat is "<<dims[0]<<"x"<<dims[1]<<std::endl;
  arma::irowvec i = Rcpp::as<arma::irowvec>(mat.slot("i"));
  arma::irowvec j = Rcpp::as<arma::irowvec>(mat.slot("j"));
  arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));
  arma::umat locations= arma::join_vert(arma::conv_to<arma::urowvec>::from(i),arma::conv_to<arma::urowvec>::from(j));


  size_t nrow = dims[0], ncol = dims[1];
  arma::sp_mat res(arma::conv_to<arma::umat>::from(locations),x,nrow,ncol);
  int bandwidth = arma::max(arma::abs(i-j));
  std::cout<<"bandwidth is "<<bandwidth<<std::endl;

  arma::mat bandmat(nrow, bandwidth+1,arma::fill::zeros);
  for(size_t ind=0; ind<=bandwidth; ind++){
    bandmat.col(ind).head(nrow-ind)= res.diag(ind);
  }
  return(bandmat);
}


//[[Rcpp::export]]
int findbandwidth(Rcpp::IntegerVector &i, Rcpp::IntegerVector &j, Rcpp::NumericVector &x,double cutoff){
  int bandwidth=0;
  for(size_t ind=0; ind<x.length(); ind++){
    int tdiff = j[ind]-i[ind];
    if(tdiff>bandwidth){
      if(fabs(x[ind])>cutoff){
        bandwidth=tdiff;
      }
    }
  }
  //std::cout<<i<<std::endl;
  return(bandwidth);
}

//[[Rcpp::export]]
int findcutoff(Rcpp::IntegerVector &i, Rcpp::IntegerVector &j, Rcpp::NumericVector &x,int bandwidth){
  double cutoff =0;
  for(size_t ind=0; ind<x.length(); ind++){
    int tdiff = j[ind]-i[ind];
    if(tdiff>=bandwidth){
      if(fabs(x[ind])>=cutoff){
        cutoff = fabs(x[ind]);
      }
    }
  }
  //std::cout<<i<<std::endl;
  return(bandwidth);
}





// [[Rcpp::export]]
size_t write_haplotype_h5(const std::string hap_gzfile,const std::string hap_h5file,const size_t nrows,const size_t ncols,size_t chunksize){

  gzFile gz_in=gzopen(hap_gzfile.c_str(),"rb");
  size_t badrows=0;
  if(chunksize>nrows){
    Rcpp::warning("chunksize is greater than number of rows!, resizing");
    chunksize=nrows;
  }
  int* haplotypes = new int[chunksize*ncols];
  //  size_t badrows;
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
      size_t i=0;
      size_t bc=0;
      for(size_t row=0; row<chunkrs; row++){
        size_t rowsum=0;
        if(row%10000==0){
          std::cout<<"row :"<<odim[0]+row<<" hap[row][0]:"<<haplotypes[i]<<std::endl;

        }

          for(size_t col=0; col<ncols; col++){
//            std::cout<<"buffer[bc]:'"<<buffer[bc]<<"'="<<(int)(buffer[bc]-'0')<<" col="<<col<<" row="<<row<<std::endl;
            haplotypes[i]=(int)(buffer[bc]-'0');
            rowsum+=haplotypes[i];
            i++;
            bc+=2;
          }
          if(rowsum==0||rowsum==ncols){
            badrows++;
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
  return(badrows);
}



//[[Rcpp::export]]
arma::uvec read_flip(const std::string hap_h5file,arma::uvec indexes){
  H5File* file;
  DataSet* dataset;
  try{
    file= new H5File( hap_h5file.c_str(), H5F_ACC_RDONLY);
  }
  catch( FileIException error )
  {
    error.printError();
  }
  Group haplogroup = file->openGroup("Legend");
  try{
    dataset = new DataSet(haplogroup.openDataSet("doFlip"));
  }
  catch( DataSetIException error )
  {
    error.printError();
  }

  DataType dt = dataset->getDataType();


  //  std::cout<<"Getting Array dimensions"<<std::endl;
  //  std::cout<<"Getting Dataspace"<<std::endl;
  DataSpace fspace =dataset->getSpace();

  //  std::cout<<"Getting Data dimensions"<<std::endl;
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t readSNPs=indexes.n_elem;
  if(readSNPs>datadim[0]){
    Rcpp::stop("indexes longer than extent of data");
  }
  hsize_t *hind = new hsize_t[readSNPs];
  for(size_t i=0; i<readSNPs;i++){
    hind[i]=indexes(i)-1;
    if(hind[i]>=datadim[0]){
      Rcpp::stop("Attempting to read outside extent of data");
    }
  }
  hsize_t rdim[1];
  rdim[0]=readSNPs;
  DataSpace memspace(1,rdim);
  fspace.selectElements(H5S_SELECT_SET,readSNPs,hind);
  //  std::cout<<"Allocating temp vector of size:"<<readSNPs<<"x"<<adim[0]<<std::endl;

  int *tdat= new int[readSNPs];
  //  std::cout<<"Reading data"<<std::endl;
  dataset->read(tdat,dt,memspace,fspace);
  arma::uvec retvec(readSNPs);
  for(size_t i=0; i<readSNPs; i++){
    retvec[i] =tdat[i];
  }
  delete [] tdat;
  memspace.close();
  fspace.close();
  haplogroup.close();
  delete dataset;
  delete file;
  return(retvec);
}

//[[Rcpp::export]]
arma::Mat<int> read_haplotype_ind_h5(const std::string hap_h5file,arma::uvec indexes){

  H5File* file;
  DataSet* dataset;
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

  //  std::cout<<"opened dataset"<<std::endl;
  hsize_t adim[1];
  DataType dt = dataset->getDataType();
  DataType dtss = dt.getSuper();
  size_t dts=dt.getSize();
  size_t bcs =dtss.getSize();
  hsize_t dims= (dts/bcs);
  adim[0]=dims;
  ArrayType mem_arraytype(PredType::NATIVE_INT32,1,adim);


  //  std::cout<<"Getting Array dimensions"<<std::endl;
  //  std::cout<<"Getting Dataspace"<<std::endl;
  DataSpace fspace =dataset->getSpace();

  //  std::cout<<"Getting Data dimensions"<<std::endl;
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  // std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
  hsize_t readSNPs=indexes.n_elem;
  if(readSNPs>datadim[0]){
    Rcpp::stop("indexes longer than extent of data");
  }
  hsize_t *hind = new hsize_t[readSNPs];
  for(size_t i=0; i<readSNPs;i++){
    hind[i]=indexes(i)-1;
    if(hind[i]>=datadim[0]){
      Rcpp::stop("Attempting to read outside extent of data");
    }
  }
  hsize_t rdim[1];
  rdim[0]=readSNPs;
  DataSpace memspace(1,rdim);
  fspace.selectElements(H5S_SELECT_SET,readSNPs,hind);
  //  std::cout<<"Allocating temp vector of size:"<<readSNPs<<"x"<<adim[0]<<std::endl;
  int *tdat= new int[readSNPs*adim[0]];
  //  std::cout<<"Reading data"<<std::endl;
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





void ip_dist(const arma::rowvec &cummapa,const arma::rowvec &cummapb,arma::mat &distmat,bool isDiag){
  std::cout<<"Doing p_dist"<<std::endl;
  distmat.set_size(cummapa.n_elem,cummapb.n_elem);
  if(isDiag){
    for(arma::uword ind=0; ind<distmat.n_rows; ind++){
      distmat.row(ind).tail(distmat.n_cols-ind-1)=cummapb.tail(cummapb.n_elem-ind-1)-cummapa(ind);
    }
  }
  else{
    for(arma::uword ind=0; ind<distmat.n_rows; ind++){
      distmat.row(ind)=cummapb-cummapa(ind);
    }
  }
}

//[[Rcpp::export]]
arma::mat ip_cov(const arma::mat &Hpanela, const arma::mat &Hpanelb, bool isDiag){
  std::cout<<"Doing p_cov"<<std::endl;
  if(isDiag){
    arma::mat covmat=trimatu(cov(Hpanela,Hpanelb));
    return(covmat);
  }
  else{
    arma::mat covmat=cov(Hpanela,Hpanelb);
    return(covmat);
  }
}


//[[Rcpp::export]]
arma::mat flip_hap(const std::string hap_h5file,arma::uvec index, const::arma::uword chunk, const arma::uword chunksize,const arma::uword nSNPs){
  arma::uword nchunks=ceil((double)nSNPs/(double)chunksize);

  if(chunk>nchunks){
    Rcpp::stop("chunk greater than nchunks! in flip_hap");
  }
  arma::uword istart=chunk*chunksize;
  arma::uword istop=std::min(((chunk+1)*chunksize)-1,nSNPs-1);
  std::cout<<"From :"<<istart<<" to: "<<istop<<std::endl;
  arma::uvec indexa= index(arma::span(istart,istop));
  arma::uvec doFlipa =read_flip(hap_h5file,indexa);
  arma::mat hmata =arma::conv_to<arma::mat>::from(read_haplotype_ind_h5(hap_h5file,indexa));
  if(hmata.n_cols!=indexa.n_elem){
    Rcpp::stop("Subsetting failed (indexa.length != mata.n_cols)");
  }
  for(size_t tcol=0; tcol<hmata.n_cols; tcol++){
    if(doFlipa(tcol)==1){
      hmata.col(tcol) =arma::abs(1-hmata.col(tcol));
    }
  }
  return(hmata);
}

//[[Rcpp::export]]
void cov_2_cor(arma::mat &covmat, arma::mat &rowvara, arma::mat &colvarb, const bool isDiag){

  rowvara = 1/sqrt(rowvara);
  colvarb=1/sqrt(colvarb);
  if(rowvara.n_rows!=covmat.n_rows){
    arma::inplace_trans(rowvara);
  }
  if(colvarb.n_cols!=covmat.n_cols){
    arma::inplace_trans(colvarb);
  }
  covmat.each_col() %=rowvara; //rowvara should have as many rows as covmat
  covmat.each_row() %=colvarb; //rowva
  if(isDiag){
    covmat.diag().ones();
  }
}



//[[Rcpp::export]]
void compute_shrinkage(arma::mat &distmat,arma::mat &S, const arma::mat &hmata ,const arma::mat &hmatb, const double theta, const double m, const double Ne,const double cutoff, const bool isDiag){
  distmat=4*Ne*distmat/100;
  distmat=exp(-distmat/(2*m));

  distmat.elem(find(distmat<cutoff)).zeros();
  distmat%=S;
  S.resize(0);
  if(isDiag){
    std::cout<<"Computing Diagonals"<<std::endl;
    distmat.diag() = arma::var(hmata);
    distmat=(1-theta)*(1-theta)*distmat+0.5*theta*(1-0.5*theta)*eye(size(distmat));
  }
  else{
    distmat*=(1-theta)*(1-theta);
  }
}
//[[Rcpp::export]]
void calcLD(arma::mat &hmata, arma::mat &hmatb, arma::rowvec &mapa, arma::rowvec &mapb,arma::mat &distmat, const double m, const double Ne,const double cutoff, const arma::uword aind, const arma::uword bind){

  std::cout<<"mata:"<<mapa.n_elem<<std::endl;
  std::cout<<"matb:"<<mapb.n_elem<<std::endl;
  double nmsum =arma::sum(1/arma::regspace<arma::vec>(1,(2*m-1)));
  double theta =(1/nmsum)/(2*m+1/nmsum);
  ip_dist(mapa,mapb,distmat,aind==bind);
  std::cout<<"Sum of distmat is "<<accu(distmat)<<std::endl;
  arma::mat S=ip_cov(hmata,hmatb,aind==bind);

  std::cout<<"Sum of covmat is "<<accu(S)<<std::endl;

  std::cout<<"Performing shrinkage"<<std::endl;
  compute_shrinkage(distmat,S, hmata , hmatb, theta, m, Ne,cutoff, aind==bind);
  std::cout<<"Sum of cormat is "<<accu(distmat)<<std::endl;
  arma::mat rowveca = arma::var(hmata)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  arma::mat colvecb= arma::var(hmatb)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  cov_2_cor(distmat,rowveca,colvecb,aind==bind);
}


//[[Rcpp::export]]
arma::sp_mat gen_sparsemat(arma::mat ldmat,const arma::uword istart,const arma::uword jstart,const arma::uword nSNPs){
  arma::umat indmat=arma::ind2sub(size(ldmat),arma::find(ldmat!=0));
  if(indmat.n_cols>0){
    //    std::cout<<"size of tmat:"<<size(tmat)<<std::endl;
    indmat.row(0)+=istart;
    indmat.row(1)+=jstart;
    std::cout<<"Allocating sparse matrix"<<std::endl;
    arma::sp_mat retS(indmat,ldmat.elem(arma::find(ldmat!=0)),nSNPs,nSNPs);
    return(retS);
  }
  else{
    arma::sp_mat retS(nSNPs,nSNPs);
    return(retS);
  }
}

//[[Rcpp::export]]
arma::sp_mat flip_hap_LD(const std::string hap_h5file, arma::uvec index,arma::rowvec map,const double m, const double Ne, const double cutoff,const arma::uword i, const arma::uword j, const arma::uword chunksize){
  std::cout<<"counting SNPs"<<std::endl;
  arma::uword nSNPs=map.n_elem;
  arma::uword nchunks=ceil((double)nSNPs/(double)chunksize);

  arma::uword istart=i*chunksize;
  arma::uword jstart=j*chunksize;

  arma::uword istop=std::min((i+1)*chunksize-1,map.n_elem-1);
  arma::uword jstop=std::min((j+1)*chunksize-1,map.n_elem-1);
  std::cout<<"subsetting map"<<std::endl;
  arma::rowvec mapa= map(arma::span(istart,istop));
  arma::rowvec mapb= map(arma::span(jstart,jstop));
  std::cout<<"subsetting haplotype data"<<std::endl;
  arma::mat hmata =flip_hap(hap_h5file,index,i,chunksize,nSNPs);
  arma::mat hmatb =flip_hap(hap_h5file,index,j,chunksize,nSNPs);
  if(hmata.n_cols!=mapa.n_elem){
    Rcpp::stop("Subsetting failed (mapa.length != mata.n_cols)");
  }
  if(hmatb.n_cols!=mapb.n_elem){
    std::cout<<"!!!!map is of length"<<mapb.n_elem<<" while hmatb is has col number of "<<hmatb.n_cols<<std::endl;
    Rcpp::stop("Subsetting failed (mapb.length != matb.n_cols)");
  }

  std::cout<<"Calculating LD"<<std::endl;
  arma::mat distmat(mapa.n_elem,mapb.n_elem,arma::fill::zeros);
  calcLD(hmata,hmatb,mapa,mapb,distmat,m,Ne,cutoff,i,j);
  std::cout<<"Checking distmat"<<std::endl;
  arma::uvec infvec =arma::find(abs(distmat)>1);
  if(infvec.n_elem>0){
    std::cout<<"incorrect correlation values found in chunk i:"<<i<<" j:"<<j<<std::endl;
    std::cout<<"Value is: "<<distmat(infvec)<<std::endl;
    arma::umat badel=arma::ind2sub(size(distmat),infvec);
    size_t overr=istart+badel(0,0);
    size_t overc=jstart+badel(1,0);
    std::cout<<"overall row is: istart+badel(0,0)="<<istart<<"+"<<badel(0,0)<<"="<<overr<<std::endl;
    std::cout<<"overall col is: jstart+badel(1,0)="<<jstart<<"+"<<badel(1,0)<<"="<<overc<<std::endl;
    Rcpp::stop("correlation values greater than 1 found in correlation matrix!");
  }
  distmat.elem(arma::find_nonfinite(distmat)).zeros();

  std::cout<<"Finding nonzero elements"<<std::endl;
  std::cout<<"Creating index matrix"<<std::endl;
  arma::sp_mat retS= gen_sparsemat( distmat, istart, jstart, nSNPs);
  return(retS);
}




// [[Rcpp::export]]
arma::Mat<int> read_haplotype_h5(const std::string hap_h5file,const size_t readSNPs,const size_t skipSNPs=0){

  H5File* file;
  DataSet* dataset;
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

//  std::cout<<"opened dataset"<<std::endl;
  hsize_t adim[1];
  DataType dt = dataset->getDataType();
  DataType dtss = dt.getSuper();
  size_t dts=dt.getSize();

  size_t bcs =dtss.getSize();
  hsize_t dims= (dts/bcs);
  adim[0]=dims;
  ArrayType mem_arraytype(PredType::NATIVE_INT32,1,adim);


//  std::cout<<"Getting Array dimensions"<<std::endl;
//  std::cout<<"Getting Dataspace"<<std::endl;
  DataSpace fspace =dataset->getSpace();

//  std::cout<<"Getting Data dimensions"<<std::endl;
  hsize_t datadim[1];
  fspace.getSimpleExtentDims(datadim,NULL);
  hsize_t rdim[1];
//  std::cout<<"Full data is of size "<<datadim[0]<<std::endl;
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
//  std::cout<<"Allocating temp vector of size:"<<rdim[0]<<"x"<<adim[0]<<std::endl;
  int *tdat= new int[rdim[0]*adim[0]];
//  std::cout<<"Reading data"<<std::endl;
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




