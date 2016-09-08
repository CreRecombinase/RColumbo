// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "H5Cpp.h"
#include "H5IO.hpp"
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

using namespace H5;





//[[Rcpp::export]]
arma::uvec gen_rows(const int i, const size_t nrows, const size_t chunksize){
  return(arma::regspace<arma::uvec>((i-1)*chunksize+1,std::min(nrows,(i)*(chunksize)))-1);
}



//[[Rcpp::export]]
arma::Mat<short> read_hap_txt(const char* inhapfile){
  arma::Mat<short> hapdat;
  hapdat.load(inhapfile,arma::raw_ascii);
  arma::inplace_trans(hapdat);
  return(hapdat);
}




//[[Rcpp::export]]
arma::mat read_hap_h5(const char* inhapfile){
  arma::mat hapdat;
  hapdat.load(inhapfile,arma::hdf5_binary);
  return(hapdat);
}

void p_dist(const arma::rowvec &cummapa,const arma::rowvec &cummapb, arma::mat &distmat,bool isDiag){
  //  std::cout<<"Doing p_dist"<<std::endl;
  if((cummapa.n_elem!=distmat.n_rows)||(cummapb.n_elem!=distmat.n_cols)){
    distmat.set_size(cummapa.n_elem,cummapb.n_elem);
  }
  if(isDiag){
    for(arma::uword i=0; i<distmat.n_rows; i++){
      distmat.row(i).tail(distmat.n_cols-i-1)=cummapb.tail(cummapb.n_elem-i-1)-cummapa(i);
    }
  }
  else{
    for(arma::uword i=0; i<distmat.n_rows; i++){
      distmat.row(i)=cummapb-cummapa(i);
    }
  }

}


void p_cov(const arma::mat &Hpanela, const arma::mat &Hpanelb, arma::mat &covmat, bool isDiag){
  //  std::cout<<"Doing p_cov"<<std::endl;
  if((Hpanela.n_cols!=covmat.n_rows)||(Hpanelb.n_cols!=covmat.n_cols)){
    covmat.set_size(Hpanela.n_cols,Hpanelb.n_cols);
  }
  if(isDiag){
    covmat=trimatu(cov(Hpanela,Hpanelb));
  }
  else{
    covmat=cov(Hpanela,Hpanelb);
  }
}




//[[Rcpp::export]]
arma::sp_mat p_sparse_LD(const arma::rowvec cummap, const arma::mat Hpanel, const double Ne, const int m,const double cutoff,const arma::uword chunksize,arma::uword i,arma::uword j){
  double nmsum =arma::sum(1/arma::regspace<arma::vec>(1,(2*m-1)));
  double theta =(1/nmsum)/(2*m+1/nmsum);
  arma::uword nSNPs=Hpanel.n_cols;
  arma::uword nchunks=ceil((double)nSNPs/(double)chunksize);
  arma::mat distmat(chunksize,chunksize);
  arma::mat S(chunksize,chunksize);
  arma::umat indmat;
  arma::vec valvec;
  arma::uvec nonz;
  arma::vec variances(nSNPs);
  std::cout<<"Starting Computation("<<nchunks<<" chunks in total, and "<<nSNPs<<" SNPs in total)"<<std::endl;

  arma::uword istart=i*chunksize;
  arma::uword jstart=j*chunksize;
  arma::uword istop=std::min((i+1)*chunksize-1,cummap.n_elem-1);
  arma::uword jstop=std::min((j+1)*chunksize-1,cummap.n_elem-1);

  std::cout<<"i From:"<<istart<<" to "<<istop<<std::endl;
  std::cout<<"j From:"<<jstart<<" to "<<jstop<<std::endl;

  p_dist(cummap(arma::span(istart,istop)),cummap(arma::span(jstart,jstop)),distmat,i==j);
  p_cov(Hpanel.cols(istart,istop),Hpanel.cols(jstart,jstop),S,i==j);
  //      std::cout<<"Making shrinkage"<<std::endl;
  distmat=4*Ne*distmat/100;
  distmat=exp(-distmat/(2*m));
  distmat.elem(find(distmat<cutoff)).zeros();
  distmat=distmat%S;
  if(i==j){
    //        std::cout<<"Computing Diagonals"<<std::endl;
    //        arma::vec mvarvec=
    //        std::cout<<"size of mvarvec:"<<size(mvarvec)<<" size of distmat.diag():"<<size(distmat.diag())<<std::endl;
    distmat.diag()=arma::var(Hpanel.cols(istart,istop));
    distmat=(1-theta)*(1-theta)*distmat+0.5*theta*(1-0.5*theta)*eye(size(distmat));
  }
  else{
    distmat=(1-theta)*(1-theta)*distmat;
  }
  nonz=arma::find(distmat!=0);
  indmat=arma::ind2sub(size(distmat),nonz);
  valvec=distmat.elem(nonz);
  //        std::cout<<"size of nonz:"<<size(nonz)<<std::endl;
  if(nonz.n_elem>0){
    arma::umat tmat=arma::ind2sub(size(distmat),nonz);
    //          std::cout<<"size of tmat:"<<size(tmat)<<std::endl;
    indmat.row(0)=indmat.row(0)+istart;
    indmat.row(1)=indmat.row(1)+jstart;
    arma::sp_mat retS(indmat,valvec,nSNPs,nSNPs);
    return(retS);
  }
  arma::sp_mat retS(nSNPs,nSNPs);
  return(retS);
}








//[[Rcpp::export]]
arma::sp_mat sparse_LD(const arma::vec cummap, const arma::mat Hpanel, const double Ne, const int m,const double cutoff, const arma::uword report_every){
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
    if(i%report_every==0){
      std::cout<<"Row :"<<i<<" of "<<nSNPs<<std::endl;
      std::cout<<"Size is :"<<row_ind.size()<<std::endl;
      std::cout<<"total possible is :"<<i*nSNPs<<std::endl;
      std::cout<<"Sparsity is :"<<((double)row_ind.size())/((double)(i*nSNPs))<<std::endl;
    }
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




//[[Rcpp::export]]
arma::umat find_bwd(arma::mat &LDmat, const double LDcutoff){
  arma::umat bandmat(LDmat.n_rows,2);
  for(size_t snp=0;snp<LDmat.n_rows;snp++){
    bandmat(snp,0)=arma::as_scalar(snp-arma::find(LDmat.col(snp)>LDcutoff,1,"first"));
    bandmat(snp,1)=arma::as_scalar(arma::find(LDmat.col(snp)>LDcutoff,1,"last")-snp);
  }
  return(bandmat);
}





//[[Rcpp::export]]
void cov_2_cor(arma::fmat &covmat, arma::fmat &rowvara, arma::fmat &colvarb, const bool isDiag){

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
void compute_shrinkage(arma::fmat &distmat,arma::fmat &S, const arma::fmat &hmata ,const arma::fmat &hmatb, const double theta, const double m, const double Ne,const double cutoff, const bool isDiag){
  distmat=4*Ne*distmat/100;
  distmat=exp(-distmat/(2*m));

  distmat.elem(find(distmat<cutoff)).zeros();
  distmat%=S;
  S.resize(0);
  if(isDiag){
    std::cout<<"Computing Diagonals"<<std::endl;
    distmat.diag() = arma::var(hmata);
    distmat=(1-theta)*(1-theta)*distmat+0.5*theta*(1-0.5*theta)*arma::eye<arma::fmat>(size(distmat));
  }
  else{
    distmat*=(1-theta)*(1-theta);
  }
}


//[[Rcpp::export]]
arma::sp_fmat gen_sparsemat(arma::fmat ldmat,const arma::uword istart,const arma::uword jstart,const arma::uword nSNPs){
  arma::umat indmat=arma::ind2sub(size(ldmat),arma::find(ldmat!=0));
  if(indmat.n_cols>0){
    //    std::cout<<"size of tmat:"<<size(tmat)<<std::endl;
    indmat.row(0)+=istart;
    indmat.row(1)+=jstart;
    std::cout<<"Allocating sparse matrix"<<std::endl;
    arma::sp_fmat retS(indmat,ldmat.elem(arma::find(ldmat!=0)),nSNPs,nSNPs);
    return(retS);
  }
  else{
    arma::sp_fmat retS(nSNPs,nSNPs);
    return(retS);
  }
}

//[[Rcpp::export]]
void calcLD(arma::fmat &hmata, arma::fmat &hmatb, arma::frowvec &mapa, arma::frowvec &mapb,arma::fmat &distmat, const double m, const double Ne,const double cutoff, const arma::uword aind, const arma::uword bind){

  std::cout<<"mata:"<<mapa.n_elem<<std::endl;
  std::cout<<"matb:"<<mapb.n_elem<<std::endl;
  double nmsum =arma::sum(1/arma::regspace<arma::vec>(1,(2*m-1)));
  double theta =(1/nmsum)/(2*m+1/nmsum);
  ip_dist(mapa,mapb,distmat,aind==bind);
  std::cout<<"Sum of distmat is "<<accu(distmat)<<std::endl;
  arma::fmat S=ip_cov(hmata,hmatb,aind==bind);

  std::cout<<"Sum of covmat is "<<accu(S)<<std::endl;

  std::cout<<"Performing shrinkage"<<std::endl;
  compute_shrinkage(distmat,S, hmata , hmatb, theta, m, Ne,cutoff, aind==bind);
  std::cout<<"Sum of cormat is "<<accu(distmat)<<std::endl;
  arma::fmat rowveca = arma::var(hmata)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  arma::fmat colvecb= arma::var(hmatb)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  cov_2_cor(distmat,rowveca,colvecb,aind==bind);
}

//[[Rcpp::export]]
arma::fmat gen_dense_LD(const std::string hap_h5file, arma::uvec index,arma::frowvec map,const double m, const double Ne, const double cutoff,const arma::uword i, const arma::uword j, const arma::uword chunksize){
  std::cout<<"counting SNPs"<<std::endl;
  arma::uword nSNPs=map.n_elem;
  arma::uword nchunks=ceil((double)nSNPs/(double)chunksize);

  arma::uword istart=i*chunksize;
  arma::uword jstart=j*chunksize;

  arma::uword istop=std::min((i+1)*chunksize-1,map.n_elem-1);
  arma::uword jstop=std::min((j+1)*chunksize-1,map.n_elem-1);
  std::cout<<"subsetting map(nSNPs is: "<<nSNPs<<")"<<std::endl;
  arma::frowvec mapa= map(arma::span(istart,istop));
  arma::frowvec mapb= map(arma::span(jstart,jstop));
  std::cout<<"subsetting haplotype data"<<std::endl;
  arma::fmat hmata =arma::conv_to<arma::fmat>::from(flip_hap(hap_h5file,index,i,chunksize,nSNPs));
  arma::fmat hmatb =arma::conv_to<arma::fmat>::from(flip_hap(hap_h5file,index,j,chunksize,nSNPs));
  if(hmata.n_cols!=mapa.n_elem){
    Rcpp::stop("Subsetting failed (mapa.length != mata.n_cols)");
  }
  if(hmatb.n_cols!=mapb.n_elem){
    std::cout<<"!!!!map is of length"<<mapb.n_elem<<" while hmatb is has col number of "<<hmatb.n_cols<<std::endl;
    Rcpp::stop("Subsetting failed (mapb.length != matb.n_cols)");
  }

  std::cout<<"Calculating LD"<<std::endl;
  arma::fmat distmat(mapa.n_elem,mapb.n_elem,arma::fill::zeros);
  calcLD(hmata,hmatb,mapa,mapb,distmat,m,Ne,cutoff,i,j);
  return(distmat);
}

//[[Rcpp::export]]
arma::sp_fmat flip_hap_LD(const std::string hap_h5file, arma::uvec index,arma::frowvec map,const double m, const double Ne, const double cutoff,const arma::uword i, const arma::uword j, const arma::uword chunksize){
  arma::uword istart=i*chunksize;
  arma::uword jstart=j*chunksize;
  arma::uword nSNPs=map.n_elem;
  arma::fmat distmat=gen_dense_LD(hap_h5file,index,map,m,Ne,cutoff,i,j,chunksize);
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
  arma::sp_fmat retS= gen_sparsemat( distmat, istart, jstart, nSNPs);
  return(retS);
}





