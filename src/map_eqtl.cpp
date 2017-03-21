#include "RcppArmadillo.h"
//[[Rcpp::depends(RcppArmadillo)]]
// #include "h5func.hpp"
//#include <cmath>
#include <vector>
#include <algorithm>

// [[Rcpp::interfaces(r,cpp)]]

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//[[Rcpp::export]]
arma::mat orthogonalize_covar(const arma::mat &Covariates){
  arma::mat ncovariates;
  arma::mat R;
  arma::qr(ncovariates,R,Covariates);
  return(ncovariates.head_cols(Covariates.n_cols));
}

//[[Rcpp::export]]
arma::mat orthogonalize_data(const arma::mat &Data, const arma::mat &covariates){
  return(Data-arma::trans((Data.t()*covariates)*covariates.t()));
}

//[[Rcpp::export]]
arma::vec calcAF(const arma::mat &Genotype){
  return(arma::conv_to<arma::vec>::from(arma::mean(Genotype)/2));
}

// [[Rcpp::export]]
arma::mat betaMatrix(const arma::mat &Genotype,const arma::mat &Expression) {
  using namespace Rcpp;
  if(Expression.n_rows!=Genotype.n_rows){
    Rcerr<<"sizes (rows) not all equal in betaMatrix()"<<std::endl;
    stop("error in betaMatrix: dimensions of matrices are incorrect!");
  }
  arma::vec genossq = arma::trans(arma::sum(arma::pow(Genotype,2)));
  Rcout<<"Calculating correlation"<<std::endl;
  arma::mat Beta = (Genotype.t()*Expression);
  Beta.each_col() /= genossq;
  return(Beta);
}


//[[Rcpp::export]]
arma::cube fastest_eQTL(const arma::mat &Genotype, const arma::mat &Expression){
  using namespace Rcpp;

  double n =Genotype.n_rows;
  Rcpp::Rcout<<"Computing correlation"<<std::endl;
  arma::mat rmat = arma::cor(Genotype,Expression);
  arma::mat Betas= betaMatrix(Genotype,Expression);
  arma::mat tstat = sqrt(n-2)*(rmat/arma::sqrt(1-arma::pow(rmat,2)));
  arma::mat semat = Betas/tstat;
  return(arma::join_slices(Betas,semat));
}

//[[Rcpp::export]]
arma::cube d_fastest_eQTL(const arma::mat &Genotype, const arma::mat &Expression){
  using namespace Rcpp;

  double n =Genotype.n_rows;
  Rcpp::Rcout<<"Computing correlation"<<std::endl;
  arma::mat rmat = arma::cor(Genotype,Expression);
  arma::mat Betas= betaMatrix(Genotype,Expression);
  arma::mat tstat = sqrt(n-2)*(rmat/arma::sqrt(1-arma::pow(rmat,2)));
  arma::mat semat = Betas/tstat;
  return(arma::join_slices(Betas,semat));
}



void fast_eqtl_lm (const arma::vec &Genotype, const arma::vec &Expression, const int n, arma::vec &retvec){

  double xtx= arma::dot(Genotype,Genotype);
  retvec[0] =dot(Expression,Genotype)/xtx;
  arma::vec resid =  Expression-Genotype*retvec[0];
  double s2 =  arma::dot(resid,resid)/(n-1);
  retvec[1]  = sqrt(s2/xtx);

}



//[[Rcpp::export]]
arma::mat eqtl_lm(const arma::mat &Genotype, const arma::vec &Expression){

  arma::mat retmat(2,Genotype.n_cols);
  int n=Expression.size();
  arma::vec resid(n);
  double xtx=0;
  double b=0;
  double s2=0;

  for(size_t i=0; i<retmat.n_cols;i++){
    xtx=arma::dot(Genotype.col(i),Genotype.col(i));
    b=arma::dot(Expression,Genotype.col(i))/xtx;
    resid=Expression-Genotype.col(i)*b;
    s2=arma::dot(resid,resid)/(n-2);
    retmat.at(0,i)=b;
    retmat.at(1,i)=sqrt(s2/xtx);
  }
  return(arma::trans(retmat));

}



// void orthogonalize_dataset(std::string h5filename,std::string newh5filename,std::string covar_h5file,std::string datagroup, std::string datasetname,std::string newdatasetname,size_t chunksize,const unsigned int deflate_level){
//   arma::mat ocovariates= arma::conv_to<arma::mat>::from(read_fmat_h5(covar_h5file,"Covardat","covariates",0,40));
//   ocovariates = join_horiz(arma::mat(ocovariates.n_rows,1,arma::fill::ones),ocovariates);
//   arma::mat Covariates = orthogonalize_covar(ocovariates);
//   size_t Nrows= get_rownum_h5(h5filename,datagroup,datasetname);
//   Rcpp::Rcout<<"total number of rows:"<<Nrows<<std::endl;
//   size_t nchunks=ceil(((double)Nrows)/(double)chunksize);
//
//   Rcpp::Rcout<<"total number of chunks:"<<nchunks<<std::endl;
//   for(size_t i=0;i<nchunks;i++){
//     Rcpp::Rcout<<"Starting chunk"<<i<<" of size:"<<chunksize<<std::endl;
//     size_t offset =i*chunksize;
//     Rcpp::Rcout<<"Reading data for chunk"<<i<<std::endl;
//     arma::mat Data = arma::conv_to<arma::mat>::from(read_fmat_h5(h5filename,datagroup,datasetname,offset,chunksize));
//     Rcpp::Rcout<<"Orthogonalizing data wrt (orthogonalized)covariates for chunk"<<i<<std::endl;
//     arma::mat oData =orthogonalize_data(Data,Covariates);
//     Rcpp::Rcout<<"oData is of dimensions:"<<oData.n_rows<<"x"<<oData.n_cols<<std::endl;
//     write_mat_h5(newh5filename,datagroup,newdatasetname,Nrows,oData.n_rows,oData,deflate_level);
//   }
//
// }
//
//



