// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include <algorithm>
#include <iterator>
#include <vector>
#include <iostream>
#include <fstream>
#include <tuple>


#define ARMA_USE_CXX11
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


// [[Rcpp::interfaces(r,cpp)]]

//[[Rcpp::export]]
arma::mat wrap_ip_dist(const arma::rowvec &cummapa,const arma::rowvec &cummapb,bool isDiag){
  arma::mat distmat(cummapa.n_elem,cummapb.n_elem,arma::fill::zeros);

//  distmat.set_size(cummapa.n_elem,cummapb.n_elem);
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
  return(distmat);
}


  void ip_dist(const arma::rowvec &cummapa,const arma::rowvec &cummapb,arma::mat &distmat,bool isDiag){
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
void compute_shrinkage(arma::mat &distmat,arma::mat &S, const arma::mat &hmata , const double theta, const double m, const double Ne,const double cutoff, const bool isDiag){
//  std::cout<<"Transforming distmat"<<std::endl;
  distmat=4*Ne*distmat/100;
//  std::cout<<"Exponentiating"<<std::endl;
  distmat=exp(-distmat/(2*m));
//  std::cout<<"Zeroing out values below cutoff"<<std::endl;
  distmat.elem(find(distmat<cutoff)).zeros();
//  std::cout<<"Multiplying by Covariance matrix"<<std::endl;
  distmat%=S;
  S.resize(0);
  if(isDiag){
//    std::cout<<"Computing Diagonals"<<std::endl;
    distmat.diag() = arma::var(hmata);
    distmat=(1-theta)*(1-theta)*distmat+0.5*theta*(1-0.5*theta)*arma::eye<arma::fmat>(size(distmat));
  }
  else{
    distmat*=(1-theta)*(1-theta);
  }
}


//[[Rcpp::export]]
arma::sp_mat gen_sparsemat(arma::mat ldmat,const arma::uword istart,const arma::uword jstart,const arma::uword nSNPs,Rcpp::LogicalVector makeSymmetric){

  if((nSNPs==ldmat.n_cols)&&(istart==1)&&(jstart==1)){
    if(makeSymmetric[0]){
      arma::sp_mat retS(symmatu(ldmat));
      return(retS);
    }else{
      arma::sp_mat retS(ldmat);
      return(retS);
    }
  }


  arma::uvec nz=arma::find(ldmat!=0);
  arma::umat indmat=arma::ind2sub(size(ldmat),nz);
  if(indmat.n_cols>0){
    //    std::cout<<"size of tmat:"<<size(tmat)<<std::endl;
    indmat.row(0)+=istart-1;
    indmat.row(1)+=jstart-1;
    std::cout<<"Allocating sparse matrix"<<std::endl;
    arma::sp_mat retS(indmat,ldmat.elem(nz),nSNPs,nSNPs);
    std::cout<<"Returning sparse matrix"<<std::endl;
    return(retS);
  }
  else{
    arma::sp_mat retS(nSNPs,nSNPs);
    return(retS);
  }
}


//[[Rcpp::export]]
arma::mat calcLD(const arma::mat hmata, const arma::mat hmatb, const arma::rowvec mapa, const arma::rowvec mapb,
            const double m, const double Ne,const double cutoff, const Rcpp::LogicalVector isDiag){
  const bool isdiag=isDiag[0];
  std::cout<<"mata:"<<mapa.n_elem<<std::endl;
  std::cout<<"matb:"<<mapb.n_elem<<std::endl;
  arma::mat distmat(mapa.n_elem,mapb.n_elem);
  double nmsum =arma::sum(1/arma::regspace<arma::vec>(1,(2*m-1)));
  double theta =(1/nmsum)/(2*m+1/nmsum);
  ip_dist(mapa,mapb,distmat,isdiag);
//    std::cout<<"Sum of distmat is "<<accu(distmat)<<std::endl;
  arma::mat S=ip_cov(hmata,hmatb,isdiag);

//  std::cout<<"Sum of covmat is "<<accu(S)<<std::endl;

  std::cout<<"Performing shrinkage"<<std::endl;
  compute_shrinkage(distmat,S, hmata, theta, m, Ne,cutoff, isdiag);
//    std::cout<<"Sum of cormat is "<<accu(distmat)<<std::endl;
  arma::mat rowveca = arma::var(hmata)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  arma::mat colvecb= arma::var(hmatb)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  cov_2_cor(distmat,rowveca,colvecb,isdiag);
  return(distmat);
}









