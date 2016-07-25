// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "H5Cpp.h"
#include<algorithm>
#include<vector>

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

using namespace H5;

//[[Rcpp::export]]
arma::mat pcov(const arma::mat &tH,const arma::uvec &c1, const arma::uvec &c2, const bool isDiag){
  arma::mat covmat = arma::cov(tH.cols(c1),tH.cols(c2));
  if(isDiag){
    return(arma::trimatu(covmat));
  }
  else{
    return(covmat);
  }
}
//[[Rcpp::export]]
arma::uvec gen_rows(const int i, const size_t nrows, const size_t chunksize){
  return(arma::regspace<arma::uvec>((i-1)*chunksize+1,std::min(nrows,(i)*(chunksize)))-1);
}

//[[Rcpp:export]]
arma::mat pdist(const arma::vec &distvec,const arma::uvec &c1, const arma::uvec &c2){
  arma::mat retmat(c1.n_elem,c2.n_elem,arma::fill::zeros);
  for(size_t j=0; j<retmat.n_cols;j++ ){
    for(size_t i=j+1;i<retmat.n_rows; i++){
      retmat(i,j)=distvec(c1(j))-distvec(c2(i));
    }
  }
  return(retmat);
}

//[[Rcpp::export]]
arma::Mat<short> read_hap_txt(const char* inhapfile){
  arma::Mat<short> hapdat;
  hapdat.load(inhapfile,arma::raw_ascii);
  arma::inplace_trans(hapdat);
  return(hapdat);
}


// In this example each column is a SNP, and each row is a sample. Chunking will be across the whole individual
//[[Rcpp::export]]
bool write_hap_h5(const arma::Mat<short> &hapdat, const char* outfile,const arma::uword chunksize=10000){
  H5File file(outfile,H5F_ACC_TRUNC);
  hsize_t dimsf[] = {hapdat.n_rows,hapdat.n_cols};
  DataSpace dataspace(2,dimsf);
  IntType datatype(PredType::NATIVE_INT);
  DSetCreatPropList cparms;
  hsize_t chunk_dims[]={hapdat.n_rows,1};
  cparms.setChunk(2,chunk_dims);
  DataSet dataset = file.createDataSet("Haplotype",datatype,dataspace,cparms);
  DataSpace fspace = dataset.getSpace();
  DataSpace mspace(2,chunk_dims);
  hsize_t offset[]={0,0};
  hsize_t block_count[]={hapdat.n_rows,1};
  for(arma::uword i=0; i<hapdat.n_cols; i++){
    fspace.selectHyperslab(H5S_SELECT_SET,block_count,offset);
    dataset.write(hapdat.colptr(i),PredType::NATIVE_INT,mspace,fspace);
    offset[1]++;
  }
  file.close();
  return(true);
}


//[[Rcpp::export]]
arma::mat read_hap_h5(const char* inhapfile){
  arma::mat hapdat;
  hapdat.load(inhapfile,arma::hdf5_binary);
  return(hapdat);
}




//[[Rcpp::export]]
arma::sp_mat sparse_LD(const arma::vec &cummap, const arma::mat &Hpanel, const double Ne, const int m,const double cutoff){
  double nmsum =arma::sum(1/arma::regspace<arma::vec>(1,(2*m-1)));
  double theta =(1/nmsum)/(2*m+1/nmsum);
  std::vector<arma::uword> row_ind;
  std::vector<arma::uword> col_ind;
  std::vector<double> covvec;
  arma::uword nSNPs=Hpanel.n_cols;
  arma::vec variances(nSNPs);
  row_ind.reserve(nSNPs*10);
  col_ind.reserve(nSNPs*10);
  covvec.reserve(nSNPs*10);

  double rho=0;
  double shrinkage=0;
  double tcov;
  for(arma::uword i=0; i < nSNPs; i++){
    for(arma::uword j = i; j<nSNPs; j++){
      if(i==j){
        tcov=arma::var(Hpanel.col(i));
        row_ind.push_back(i);
        col_ind.push_back(j);
        covvec.push_back(tcov*(1-0.5*theta));
        variances(i)=tcov*(1-0.5*theta);
      }
      else{
        rho=4*Ne*(cummap(j)-cummap(i))/100;
        shrinkage=exp(-rho/(2*m));
        if(shrinkage>=cutoff){
          tcov=as_scalar(arma::cov(Hpanel.col(j),Hpanel.col(i)));
          row_ind.push_back(i);
          col_ind.push_back(j);
          covvec.push_back((tcov*shrinkage)*((1-theta)*(1-theta)));
        }
      }
    }
  }
  for(arma::uword i=0; i<covvec.size();i++){
    covvec[i]=covvec[i]/(sqrt(variances[row_ind[i]])*sqrt(variances[col_ind[i]]));
  }
  arma::umat xymat=join_rows(arma::uvec(&row_ind[0],row_ind.size()),arma::uvec(&col_ind[0],col_ind.size()));
  arma::inplace_trans(xymat);
  arma::sp_mat sighat(xymat,arma::vec(&covvec[0],covvec.size(),false),nSNPs,nSNPs);
  return(sighat);
}


//
// struct bigCov: public Worker {
//   //input matrix
//   const arma::mat Hpanel;
//
//   //output matrix
//   arma::mat S;
//   size_t chunksize;
//   size_t nchunks;
//   arma::uvec rvec;
//   arma::uvec cvec;
//   size_t nrows;
//
//   bigCov(const arma::mat inpm,const arma::mat rmat,const int nchunk, const int chunks) :Hpanel(inpm),S(rmat){
//     nrows=inpm.n_rows;
//     chunksize=chunks;
//     nchunks=nchunk;
//     size_t diag_mats=ceil(((double)nrows)/((double)chunksize));
//     size_t total_mats = pow(diag_mats,2);
//     size_t matnum = (total_mats-diag_mats)/2+diag_mats;
//     rvec= arma::ones<arma::uvec>(matnum);
//     cvec= arma::ones<arma::uvec>(matnum);
//     int t=0;
//     for(int i=0;i<nchunks; i++){
//       for(int j=0; j<nchunks;j++){
//         if(j>=i){
//           rvec(t)=i;
//           cvec(t)=j;
//           t++;
//         }
//       }
//     }
//
//   }
//
//   void operator()(std::size_t begin, std::size_t end){
//
//   }
// };

