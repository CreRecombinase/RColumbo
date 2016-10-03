// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "H5Cpp.h"
#include "H5IO.hpp"
#include "h5func.hpp"
#include "snp_exp.hpp"
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
arma::uvec greedy_ind_LD(const arma::fmat &LDmat, const float LDcutoff,const arma::uword offset){
  arma::uword nc=LDmat.n_cols;
  arma::uvec indvec(nc,arma::fill::zeros);
  arma::uword i=offset;
  indvec[0]=i;
  arma::uvec ti;
  arma::uvec nti;
  arma::frowvec cuLD;
  arma::fmat tmat=abs(arma::trimatl(LDmat));
  arma::fmat nfmat;
  arma::uword ielem=1;
  while(i<LDmat.n_cols){
    Rcpp::Rcout<<i<<std::endl;
    ti=arma::find(tmat.col(i).tail(nc-1-i)<LDcutoff);
    if(ti.n_elem>=1){
      cuLD = arma::min(tmat.submat(ti+i+1,indvec.head(ielem)),0);
      nti = arma::find(cuLD<LDcutoff,1);
      if(nti.n_elem==1){
        i=arma::as_scalar(nti+1+i);
        ielem++;
        indvec[ielem]=i;
      }
      else{
        break;
      }
    }else{
      break;
    }
  }
  return(arma::conv_to<arma::uvec>::from(indvec));
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
void compute_shrinkage(arma::fmat &distmat,arma::fmat &S, const arma::fmat &hmata , const double theta, const double m, const double Ne,const double cutoff, const bool isDiag){
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
arma::fmat slowLD(const arma::fmat &hmata, const arma::fmat &hmatb, const arma::frowvec &mapa, const arma::frowvec &mapb, const double m, const double Ne,const double cutoff, const arma::uword aind, const arma::uword bind){
  arma::fmat distmat(mapa.n_elem,mapb.n_elem);
  std::cout<<"mata:"<<mapa.n_elem<<std::endl;
  std::cout<<"matb:"<<mapb.n_elem<<std::endl;
  double nmsum =arma::sum(1/arma::regspace<arma::vec>(1,(2*m-1)));
  double theta =(1/nmsum)/(2*m+1/nmsum);
  ip_dist(mapa,mapb,distmat,aind==bind);
  //  std::cout<<"Sum of distmat is "<<accu(distmat)<<std::endl;
  arma::fmat S=ip_cov(hmata,hmatb,aind==bind);

  //  std::cout<<"Sum of covmat is "<<accu(S)<<std::endl;

  //  std::cout<<"Performing shrinkage"<<std::endl;
  compute_shrinkage(distmat,S, hmata, theta, m, Ne,cutoff, aind==bind);
  //  std::cout<<"Sum of cormat is "<<accu(distmat)<<std::endl;
  arma::fmat rowveca = arma::var(hmata)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  arma::fmat colvecb= arma::var(hmatb)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  cov_2_cor(distmat,rowveca,colvecb,aind==bind);
  return(distmat);
}


//[[Rcpp::export]]
void calcLD(const arma::fmat &hmata, const arma::fmat &hmatb, const arma::frowvec &mapa, const arma::frowvec &mapb,arma::fmat &distmat, const double m, const double Ne,const double cutoff, const arma::uword aind, const arma::uword bind){

  std::cout<<"mata:"<<mapa.n_elem<<std::endl;
  std::cout<<"matb:"<<mapb.n_elem<<std::endl;
  double nmsum =arma::sum(1/arma::regspace<arma::vec>(1,(2*m-1)));
  double theta =(1/nmsum)/(2*m+1/nmsum);
  ip_dist(mapa,mapb,distmat,aind==bind);
  //  std::cout<<"Sum of distmat is "<<accu(distmat)<<std::endl;
  arma::fmat S=ip_cov(hmata,hmatb,aind==bind);

  //  std::cout<<"Sum of covmat is "<<accu(S)<<std::endl;

  //  std::cout<<"Performing shrinkage"<<std::endl;
  compute_shrinkage(distmat,S, hmata, theta, m, Ne,cutoff, aind==bind);
  //  std::cout<<"Sum of cormat is "<<accu(distmat)<<std::endl;
  arma::fmat rowveca = arma::var(hmata)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  arma::fmat colvecb= arma::var(hmatb)*(1-theta)*(1-theta)+0.5*theta*(1-0.5*theta);
  cov_2_cor(distmat,rowveca,colvecb,aind==bind);
}



//[[Rcpp::export]]
arma::fmat gen_chunk_LD(const std::string hap_h5file,arma::frowvec map,const double m, const double Ne, const double cutoff,const arma::uword i, const arma::uword j, const arma::uword chunksize){
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
  arma::fmat hmata=read_fmat_h5(hap_h5file,"Haplotype","genotype",istart,mapa.n_elem);
  arma::fmat hmatb=read_fmat_h5(hap_h5file,"Haplotype","genotype",jstart,mapb.n_elem);
  // arma::fmat hmata =arma::conv_to<arma::fmat>::from(flip_hap(hap_h5file,index,i,chunksize,nSNPs));
  // arma::fmat hmatb =arma::conv_to<arma::fmat>::from(flip_hap(hap_h5file,index,j,chunksize,nSNPs));
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
arma::fmat chunk_LD(const std::string h5file,const arma::frowvec &tmap,const size_t offset,const size_t chunksize,const double m, const double Ne, const double cutoff){
  std::cout<<"counting SNPs"<<std::endl;
  // arma::frowvec map=arma::conv_to<arma::frowvec>::from(read_float_h5(h5file,"Legend","cummap",offset,chunksize));
  arma::uword istart=offset;
  arma::uword istop=std::min(offset+(arma::uword)chunksize-1,tmap.n_elem-1);
  arma::frowvec map=tmap(arma::span(istart,istop));
  // arma::uword nSNPs=map.n_elem;
  // arma::uword nchunks=ceil((double)nSNPs/(double)chunksize);

//  std::cout<<"subsetting haplotype data"<<std::endl;
arma::fmat hmat=read_fmat_h5(h5file,"Haplotype","haplotype",offset,map.n_elem);
  if(hmat.n_cols!=map.n_elem){
    Rcpp::stop("Subsetting failed (mapa.length != mata.n_cols)");
  }
//  std::cout<<"Calculating LD"<<std::endl;
  arma::fmat distmat(map.n_elem,map.n_elem,arma::fill::zeros);
  calcLD(hmat,hmat,map,map,distmat,m,Ne,cutoff,0,0);
  return(distmat);
}




//[[Rcpp::export]]
arma::fmat fslide_LD(const std::string h5file, const size_t chunksize,const size_t offset,const float cutoff){
    
    float m=85;
    float Ne=11490.672741;
    size_t dchunksize=chunksize*2;
    LD_dataset lddat(h5file);
    lddat.set_offset(offset);
    arma::frowvec mapa=arma::trans(lddat.get_map_chunk(chunksize));
    arma::fmat mata=lddat.read_chunk(chunksize);
    if(mapa.n_elem!=mata.n_cols){
        Rcpp::stop("mapa n_elem!=mata.n_cols");
    }
    Rcpp::Rcout<<"distmatrix being generated"<<std::endl;
    arma::fmat distmata(mapa.n_elem,mapa.n_elem);
    Rcpp::Rcout<<"Calculating LD..."<<std::endl;
    calcLD(mata,mata,mapa,mapa,distmata,m,Ne,cutoff,0,0);
    size_t noffset=lddat.increment_offset(chunksize);
    if(noffset==lddat.P){
        return(distmata);
    }
    arma::frowvec mapb=arma::trans(lddat.get_map_chunk(chunksize));
    arma::fmat matb=lddat.read_chunk(chunksize);
    if(mapb.n_elem!=matb.n_cols){
        Rcpp::stop("mapb n_elem!=matb.n_cols");
    }
    arma::fmat distmatb(mapa.n_elem,mapb.n_elem);
    calcLD(mata,matb,mapa,mapb,distmatb,m,Ne,cutoff,0,1);
    return(arma::join_horiz(distmata,distmatb));
}

size_t write_cov_LD(const std::string oh5file,arma::fmat &rect_covmat,const size_t rowoffset,const size_t coloffset,const size_t dimension,const size_t rowchunksize,const size_t colchunksize){
    write_covmat_h5(oh5file,"LD_mat","LD",dimension,rect_covmat,rowoffset,coloffset,rowchunksize,colchunksize);
  return(rect_covmat.n_cols);
}

//[[Rcpp::export]]
int Rwrite_cov_LD(Rcpp::String th5file, Rcpp::IntegerVector tdimension,Rcpp::IntegerVector Offset, Rcpp::NumericMatrix data,Rcpp::IntegerVector chunkvec){
  std::string h5file=th5file;
  size_t dimension=tdimension[0];
  size_t rowoffset=Offset[0];
  size_t coloffset=Offset[0];
  std::vector<size_t> chunkdims(2);
  chunkdims[0]=chunkvec[0];
  chunkdims[1]=chunkvec[1];
  arma::fmat tdat = Rcpp::as<arma::fmat>(data);
  size_t ret=write_cov_LD(h5file,tdat,rowoffset,coloffset,dimension,chunkdims[0],chunkdims[1]);
  return(ret);
}


//[[Rcpp::export]]
int Rwrite_blosc_cov_LD(Rcpp::String th5file, Rcpp::IntegerVector tdimension,Rcpp::IntegerVector Offset, Rcpp::NumericMatrix data){
  std::string h5file=th5file;
  size_t dimension=tdimension[0];
  size_t rowoffset=Offset[0];
  size_t coloffset=Offset[0];
  arma::fmat tdat = Rcpp::as<arma::fmat>(data);
  size_t ret=write_covmat_h5(h5file,"LD_mat","LD",dimension,tdat,rowoffset,coloffset,1000,1000);
  return(ret);
}
//[[Rcpp::export]]
arma::fmat slide_LD(const std::string h5file,const arma::frowvec &tmap,const size_t offset,const size_t chunksize,const double m, const double Ne, const double cutoff){


  arma::uword astart=offset;
  arma::uword astop=std::min(astart+(arma::uword)chunksize-1,tmap.n_elem-1);
  if(astop==tmap.n_elem-1){
    return(chunk_LD(h5file,tmap,offset,chunksize,m,Ne,cutoff));
  }
  arma::uword bstart=offset+chunksize;
  arma::uword bstop=std::min(bstart+chunksize-1,tmap.n_elem-1);

  arma::frowvec mapa=tmap(arma::span(astart,astop));
  arma::frowvec mapb=tmap(arma::span(bstart,bstop));
  if(mapa.n_elem!=mapb.n_elem){
    Rcpp::Rcout<<"mapa has: "<<mapa.n_elem<<std::endl;
    Rcpp::Rcout<<"mapb has: "<<mapb.n_elem<<std::endl;
    Rcpp::warning("Mapa has a different number of elements than mapb");
  }

  //  std::cout<<"subsetting haplotype data"<<std::endl;
  arma::fmat hmat=read_fmat_h5(h5file,"Haplotype","haplotype",offset,mapa.n_elem+mapb.n_elem);
  //  std::cout<<"Calculating LD"<<std::endl;
  arma::fmat distmata(mapa.n_elem,mapa.n_elem,arma::fill::zeros);
  calcLD(hmat.head_cols(mapa.n_elem),hmat.head_cols(mapa.n_elem),mapa,mapa,distmata,m,Ne,cutoff,0,0);
  //  arma::umat keepinda= arma::ind2sub(arma::size(distmat),arma::find(abs(distmat)>LDcutoff));
  //  keepinda+=offset;
  arma::fmat distmatb(mapa.n_elem,mapb.n_elem);
  calcLD(hmat.head_cols(mapa.n_elem),hmat.tail_cols(mapb.n_elem),mapa,mapb,distmatb,m,Ne,cutoff,0,1);
  // arma::umat keepindb= arma::ind2sub(arma::size(distmat),arma::find(abs(distmat)>LDcutoff));
  // keepindb+=offset;
  // keepindb.row(1)+=chunksize;
  return(arma::join_horiz(distmata,distmatb));
}


//[[Rcpp::export]]
arma::fmat slide_LD_ind(const std::string h5file,const arma::frowvec &tmap,const arma::uvec index,const size_t offset,const size_t chunksize,const double m, const double Ne, const double cutoff){

  if(index.n_elem!=tmap.n_elem){
    Rcpp::stop("index contains a different number of elements than tmap");
  }
  arma::uword astart=offset;
  arma::uword astop=std::min(astart+(arma::uword)chunksize-1,tmap.n_elem-1);
  if(astop==tmap.n_elem-1){
    arma::frowvec mapa=tmap(arma::span(astart,astop));
    arma::uvec inda=index(arma::span(astart,astop));
    arma::fmat hmat=read_fmat_chunk_ind(h5file,"Haplotype","haplotype",inda);
    arma::fmat distmata(mapa.n_elem,mapa.n_elem,arma::fill::zeros);
    calcLD(hmat,hmat,mapa,mapa,distmata,m,Ne,cutoff,0,0);
    return(distmata);
  }
  arma::uword bstart=offset+chunksize;
  arma::uword bstop=std::min(bstart+chunksize-1,tmap.n_elem-1);

  arma::frowvec mapa=tmap(arma::span(astart,astop));
  arma::frowvec mapb=tmap(arma::span(bstart,bstop));
  arma::uvec inda=index(arma::span(astart,astop));
  arma::uvec indb=index(arma::span(bstart,bstop));
  arma::uvec full_ind =join_vert(inda,indb);
  if(mapa.n_elem!=mapb.n_elem){
    Rcpp::Rcout<<"mapa has: "<<mapa.n_elem<<std::endl;
    Rcpp::Rcout<<"mapb has: "<<mapb.n_elem<<std::endl;
    Rcpp::warning("Mapa has a different number of elements than mapb");
  }

  //  std::cout<<"subsetting haplotype data"<<std::endl;
  arma::fmat hmat=read_fmat_chunk_ind(h5file,"Haplotype","haplotype",full_ind);
  //  std::cout<<"Calculating LD"<<std::endl;
  arma::fmat distmata(mapa.n_elem,mapa.n_elem,arma::fill::zeros);
  calcLD(hmat.head_cols(mapa.n_elem),hmat.head_cols(mapa.n_elem),mapa,mapa,distmata,m,Ne,cutoff,0,0);
  //  arma::umat keepinda= arma::ind2sub(arma::size(distmat),arma::find(abs(distmat)>LDcutoff));
  //  keepinda+=offset;
  arma::fmat distmatb(mapa.n_elem,mapb.n_elem);
  calcLD(hmat.head_cols(mapa.n_elem),hmat.tail_cols(mapb.n_elem),mapa,mapb,distmatb,m,Ne,cutoff,0,1);
  // arma::umat keepindb= arma::ind2sub(arma::size(distmat),arma::find(abs(distmat)>LDcutoff));
  // keepindb+=offset;
  // keepindb.row(1)+=chunksize;
  return(arma::join_horiz(distmata,distmatb));
}


arma::umat cutoff_LD(const arma::fmat &LD_mat, const size_t offset,const float LD_cutoff){
    return(arma::trans(arma::ind2sub(arma::size(LD_mat),arma::find(abs(LD_mat)>LD_cutoff))+offset));
}

//[[Rcpp::export]]
arma::umat slide_cutoff_LD(const std::string h5file,const arma::frowvec &tmap,const size_t offset,const size_t chunksize,const double m, const double Ne, const double cutoff,const float LD_cutoff){
  return(cutoff_LD(slide_LD(h5file,tmap,offset,chunksize,m,Ne,cutoff),offset,LD_cutoff));
}

// arma::umat network_chrom_LDmat(const std::string hap_h5file,arma::frowvec &tmap, const size_t chunksize, const double m, const double Ne,const double cutoff){
//
// }

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








