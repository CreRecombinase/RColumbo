#include "RcppArmadillo.h"
//[[Rcpp::depends(RcppArmadillo)]]
#include "h5func.hpp"
//#include <cmath>
#include <vector>
#include <algorithm>


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
arma::fmat orthogonalize_covar(const arma::fmat &Covariates){
  arma::fmat ncovariates;
  arma::fmat R;
  arma::qr(ncovariates,R,Covariates);
  return(ncovariates.head_cols(Covariates.n_cols));
}

//[[Rcpp::export]]
arma::fmat orthogonalize_data(const arma::fmat &Data, const arma::fmat &covariates){
  return(Data-arma::trans((Data.t()*covariates)*covariates.t()));
}

// [[Rcpp::export]]
arma::fmat betaMatrix(const arma::fmat &Genotype,const arma::fmat &Expression) {
  using namespace Rcpp;
  if(Expression.n_rows!=Genotype.n_rows){
    Rcerr<<"sizes (rows) not all equal in betaMatrix()"<<std::endl;
    stop("error in betaMatrix: dimensions of matrices are incorrect!");
  }
  arma::fvec genossq = arma::trans(arma::sum(arma::pow(Genotype,2)));
  Rcout<<"Calculating correlation"<<std::endl;
  arma::fmat Beta = (Genotype.t()*Expression);
  Beta.each_col() /= genossq;
  return(Beta);
}


//[[Rcpp::export]]
arma::fcube fastest_eQTL(const std::string genotypef, const arma::uvec &snpinter, const std::string expressionf, const arma::uvec expinter){
  using namespace Rcpp;

  arma::fmat Genotype = read_fmat_chunk_ind(genotypef,"SNPdata","genotype",snpinter);
  arma::fmat Expression = read_fmat_chunk_ind(expressionf,"EXPdata","expression",expinter);

  double n =Genotype.n_rows;
  Rcpp::Rcout<<"Computing correlation"<<std::endl;
  arma::fmat rmat = arma::cor(Genotype,Expression);
  arma::fmat Betas= betaMatrix(Genotype,Expression);
  arma::fmat tstat = sqrt(n-2)*(rmat/arma::sqrt(1-arma::pow(rmat,2)));
  return(arma::join_slices(Betas,tstat));
}


//[[Rcpp::export]]
void orthogonalize_dataset(std::string h5filename,std::string newh5filename,std::string covar_h5file,std::string datagroup, std::string datasetname,std::string newdatasetname,size_t chunksize,const unsigned int deflate_level){
  arma::fmat ocovariates= read_fmat_h5(covar_h5file,"Covardat","covariates",0,40);
  ocovariates = join_horiz(arma::fmat(ocovariates.n_rows,1,arma::fill::ones),ocovariates);
  arma::fmat Covariates = orthogonalize_covar(ocovariates);
  size_t Nrows= get_rownum_h5(h5filename,datagroup,datasetname);
  Rcpp::Rcout<<"total number of rows:"<<Nrows<<std::endl;
  size_t nchunks=ceil(((double)Nrows)/(double)chunksize);

  Rcpp::Rcout<<"total number of chunks:"<<nchunks<<std::endl;
  for(size_t i=0;i<nchunks;i++){
    Rcpp::Rcout<<"Starting chunk"<<i<<" of size:"<<chunksize<<std::endl;
    size_t offset =i*chunksize;
    Rcpp::Rcout<<"Reading data for chunk"<<i<<std::endl;
    arma::fmat Data = read_fmat_h5(h5filename,datagroup,datasetname,offset,chunksize);
    Rcpp::Rcout<<"Orthogonalizing data wrt (orthogonalized)covariates for chunk"<<i<<std::endl;
    arma::fmat oData =orthogonalize_data(Data,Covariates);
    Rcpp::Rcout<<"oData is of dimensions:"<<oData.n_rows<<"x"<<oData.n_cols<<std::endl;
    write_mat_h5(newh5filename,datagroup,newdatasetname,Nrows,oData.n_rows,oData,deflate_level);
  }

}





