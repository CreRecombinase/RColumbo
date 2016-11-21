
#include "RcppArmadillo.h"
#include <H5Cpp.h>
#include "H5IO.hpp"
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <numeric>
#include "h5func.hpp"
#include "snp_exp.hpp"
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

// class Legend{
//arma::uvec chromposmap; //vector of length 22 giving cumulative position of last position of the chromomosme (from dbSNP)
// public:
//   arma::uvec chrom;
//   arma::uvec pos;
//   arma::uvec abspos;
//   Legend(const arma::uvec tchrom,const arma::uvec tpos, const arma::uvec chrompos);
//   Legend subset_index(const arma::uvec index);
//   Legend subset_chrom(const arma::uword chromi);
// };

//[[Rcpp::export]]
arma::uvec make_long(arma::uvec &vchrom, arma::uvec &vpos,arma::uvec &posmap){
  if(vchrom.n_elem!=vpos.n_elem){
    Rcpp::stop("length of vectors not equal in make_long");
  }
  arma::uvec vlong =posmap.elem(vchrom-1)+vpos;
  return(vlong);
}
// Legend::Legend(){
//   chrom.set_size(0);
//   pos.set_size(0);
//   abspos.set_size(0);
// }
//
// Legend::Legend(const arma::uvec tchrom,const arma::uvec tpos, const arma::uvec chrompos):chrom(tchrom),pos(tpos),chromposmap(chrompos){
//   abspos=chromposmap.elem(chrom-1)+pos;
//   n=pos.n_elem;
//   for(size_t i=0; i<n;i++){
//     absposmap[abspos[i]]=i;
//   }
// }
//
// Legend Legend::subset_index(const arma::uvec &index){
//   arma::uvec npos=pos.elem(index);
//   arma::uvec nchrom=chrom.elem(index);
//   return(Legend(nchrom,npos,chromposmap));
// }
//
// Legend Legend::subset_chrom(const arma::uword chromi){
//   arma::uvec indv=arma::find(chrom==chromi);
//   return(subset_index(indv));
// }
//
// Legend Legend::subset_abspos(const arma::uvec &absp){
//   std::vector<arma::uword> indexvec;
//   indexvec.reserve(absp.size());
//   arma::uword vpos=0;
//   std::unordered_map<arma::uword,arma::uword>::iterator fit;
//   for(arma::uvec::const_iterator rit=absp.begin(); rit!=absp.end(); rit++){
//     fit =absposmap.find(*rit);
//     if(fit!=absposmap.end()){
//       indexvec.push_back(fit->second);
//     }
//   }
//   arma::uvec idv(&indexvec[0],indexvec.size());
//   return(subset_index(idv));
// }
//
//
// dbsnp::dbsnp(std::string h5file){
//   arma::uvec posv= arma::conv_to<arma::uvec>::from(read_uint_h5(h5file,"dbSNP","pos"));
//   arma::uvec chromv= arma::conv_to<arma::uvec>::from(read_int_h5(h5file,"dbSNP","chrom"));
//   arma::uvec rsidv = arma::conv_to<arma::uvec>::from(read_uint_h5(h5file,"dbSNP","rsid"));
//   arma::uvec cs(23,arma::fill::zeros);
//   for(size_t i=0; i<chromv.size();i++){
//     rsid2posmap[rsidv[i]]=i;
//     if(cs[chromv[i]]<posv[i]){
//       cs[chromv[i]]= posv[i];
//     }
//   }
//   arma::uvec chromposmap=cumsum(cs);
//   legend = Legend(chromv,posv,chromposmap);
//   pos2rsidmap=make_map(chromv,posv,rsidv,chromposmap);
//
//
// }
//
// arma::umat dbsnp::find_rsid(const Legend &qlegend){
//   std::vector<arma::uword> frsvec;
//   std::vector<arma::uword> indexvec;
//   frsvec.reserve(qlegend.size());
//   indexvec.reserve(qlegend.size());
//   arma::uword vpos=0;
//   std::unordered_map<arma::uword,arma::uword>::iterator fit;
//   std::cout<<"mapping rsid"<<std::endl;
//   for(arma::uvec::const_iterator rit=qlegend.abspos.begin(); rit!=qlegend.abspos.end(); rit++){
//     fit =pos2rsidmap.find(*rit);
//     vpos = rit-qlegend.abspos.begin();
//     if(fit!=pos2rsidmap.end()){
//       indexvec.push_back(vpos);
//       frsvec.push_back(fit->second);
//     }
//   }
//   arma::uvec rsv(&frsvec[0],frsvec.size());
//   arma::uvec idv(&indexvec[0],indexvec.size());
//   return(arma::join_horiz(rsv,idv));
// }
//
// Legend dbsnp::find_legend(const arma::uvec &rsidind){
//
//   std::vector<arma::uword> indexvec;
//   indexvec.reserve(rsidind.size());
//   std::unordered_map<arma::uword,arma::uword>::iterator fit;
//   for(arma::uvec::const_iterator rit=rsidind.begin(); rit!=rsidind.end(); rit++){
//     fit =rsid2posmap.find(*rit);
//     if(fit!=rsid2posmap.end()){
//       indexvec.push_back(fit->second);
//     }
//   }
//   arma::uvec idv(&indexvec[0],indexvec.size());
//   return(legend.subset_index(idv));
// }
//


// class GWAS{
// private:
//   size_t n;
// public:
//   Legend legend;
//   arma::vec betahat;
//   arma::vec serr;
//   GWAS(const std::string gwash5);
//
// // };
// GWAS::GWAS(const std::string gwash5,const std::string dbsnph5){
//   betahat= arma::conv_to<arma::vec>::from(read_float_h5(gwash5,"GWAS","beta"));
//   serr= arma::conv_to<arma::vec>::from(read_float_h5(gwash5,"GWAS","serr"));
//   rsid= arma::conv_to<arma::uvec>::from(read_int_h5(gwash5,"GWAS","rsid"));
//   dbsnp dbm(dbsnph5);
//   legend =dbm.find_legend(rsid);
// }



//
//
// }
//





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
arma::fmat read_fmat_gz(const std::string gzfile,const size_t chunksize, const size_t nrows, const arma::uvec keeprow,const size_t ncols){
  size_t keeprown=sum(keeprow);
  gzFile gz_in=gzopen(gzfile.c_str(),"rb");
  arma::fmat datamat(keeprown,ncols);
  size_t bytes_read=0;
  size_t len=chunksize*ncols*2;
  size_t datarow=0;
  char *buffer= new char[len];
  size_t bc=0;
  size_t totrow=0;
  for(;;){
    bytes_read=gzread(gz_in,buffer,len);
    bc=0;
    if(bytes_read==0){
      break;
    }
    size_t chunkrs=bytes_read/(ncols*2);
    for(size_t row=0; row<chunkrs; row++){
      if(totrow%10000==0){
        std::cout<<"row :"<<totrow<<std::endl;
      }
      for(size_t col=0; col<ncols; col++){
        //            std::cout<<"buffer[bc]:'"<<buffer[bc]<<"'="<<(int)(buffer[bc]-'0')<<" col="<<col<<" row="<<row<<std::endl;
        if(keeprow[totrow]==1){
          datamat(datarow,col)=(float)(buffer[bc]-'0');
        }
        bc+=2;
      }
      if(keeprow[totrow]==1){
        datarow++;
      }
      totrow++;
    }
  }
  return(datamat);
}

// [[Rcpp::export]]
size_t write_haplotype_h5(const std::string hap_gzfile,const std::string hap_h5file,const size_t nrows,const size_t ncols,size_t chunksize,const unsigned int deflate_level){

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
    cparms.setDeflate(deflate_level);

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
    Rcpp::Rcerr<<"data is of length:"<<datadim[0]<<" and there are "<<readSNPs<<" in the indexes"<<std::endl;
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




//[[Rcpp::export]]
arma::mat read_dmat_ind_h5(const std::string hap_h5file,const std::string groupname, const std::string dataname,arma::uvec indexes){

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
  hsize_t adim[1];
  DataType dt = dataset->getDataType();
  DataType dtss = dt.getSuper();
  size_t dts=dt.getSize();
  size_t bcs =dtss.getSize();
  hsize_t dims= (dts/bcs);
  adim[0]=dims;
  ArrayType mem_arraytype(PredType::NATIVE_DOUBLE,1,adim);


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
      Rcpp::Rcerr<<"i is: "<<i<<", hind[i] is: "<<hind[i]<<"datadim[0] is:"<<datadim[0]<<std::endl;
      Rcpp::stop("Attempting to read outside extent of data");
    }
  }
  hsize_t rdim[1];
  rdim[0]=readSNPs;
  DataSpace memspace(1,rdim);
  fspace.selectElements(H5S_SELECT_SET,readSNPs,hind);
  //  std::cout<<"Allocating temp vector of size:"<<readSNPs<<"x"<<adim[0]<<std::endl;
  arma::mat retmat(adim[0],rdim[0]);
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

arma::uvec ind_lookup(const arma::uvec &queryvec, const arma::uvec &targetvec){
  arma::uvec result(queryvec.n_elem);
  for(size_t i = 0;i<queryvec.n_elem;i++){
    arma::uvec::const_iterator it=std::find(targetvec.begin(),targetvec.end(),queryvec[i]);
    if(it!=targetvec.end()){
    result[i]= it-targetvec.begin();
    }else{
      result[i]=-1;
    }
  }
  result=result.elem(arma::find(result>0));
  return(result);
}

//[[Rcpp::export]]
arma::mat read_dmat_rowname(const std::string h5file,const std::string annogroupname, const std::string annocolname,const std::string datagroupname, const std::string datacolname,arma::uvec queryvec){
  arma::uvec annocol = arma::conv_to<arma::uvec>::from(read_uint_h5(h5file,annogroupname,annocolname));
  arma::uvec index = ind_lookup(queryvec,annocol)+1;
  return(read_dmat_chunk_ind(h5file,datagroupname,datacolname,index));
}

//[[Rcpp::export]]
arma::fmat read_fmat_rowname(const std::string h5file,const std::string annogroupname, const std::string annocolname,const std::string datagroupname, const std::string datacolname,arma::uvec queryvec){
  arma::uvec annocol = arma::conv_to<arma::uvec>::from(read_uint_h5(h5file,annogroupname,annocolname));
  arma::uvec index = ind_lookup(queryvec,annocol)+1;
  return(read_fmat_chunk_ind(h5file,datagroupname,datacolname,index));
}


// void tbb_dist(const arma::frowvec &cummapa,const arma::frowvec &cummapb,arma::fmat &distmat,bool isDiag){
//   using namespace tbb;
// #pragma warning(disable: 588)
//   distmat.set_size(cummapa.n_elem,cummapb.n_elem);
//   size_t n=distmat.n_rows;
//   size_t itern=(n*n-n)/2;
//   if(isDiag){
//     parallel_for(size_t(0),itern,[&](size_t ind){
//       size_t i= n-2-floor(sqrt(-8*ind+4*n*(n-1)-7)/2.0-0.5);
//       size_t j= ind+i+1-n*(n-1)/2+(n-i)*((n-i)-1)/2;
//       distmat.at(i,j)=cummapb[j]-cummapa[i];
//       });
//   }
//   else{
//     parallel_for(size_t(0),n,[&](size_t ind){
//       distmat.row(ind)=cummapb-cummapa(ind);
//     });
//   }
// }


void ip_dist(const arma::frowvec &cummapa,const arma::frowvec &cummapb,arma::fmat &distmat,bool isDiag){
//  std::cout<<"Doing p_dist"<<std::endl;
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
  std::cout<<"Done with p_dist"<<std::endl;
}

//[[Rcpp::export]]
arma::fmat ip_cov(const arma::fmat &Hpanela, const arma::fmat &Hpanelb, bool isDiag){
//  std::cout<<"Doing p_cov"<<std::endl;
  if(isDiag){
    arma::fmat covmat=trimatu(cov(Hpanela,Hpanelb));
    return(covmat);
  }
  else{
    arma::fmat covmat=cov(Hpanela,Hpanelb);
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

  arma::uvec indexa= index(arma::span(istart,istop));
  std::cout<<"From:"<<istart<<" to: "<<istop<<"("<<indexa[0]<<":"<<indexa[indexa.size()-1]<<")"<<std::endl;
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





// [[Rcpp::export]]
arma::Mat<int> read_haplotype_h5(const std::string hap_h5file,const size_t readSNPs,const size_t skipSNPs){

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




