#include "RcppArmadillo.h"
#include "h5func.hpp"
#include <cmath>
#include <vector>
#include <algorithm>
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

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
arma::fmat rMatrix(const arma::fmat &Genotype,const arma::fmat &Expression){
  return(arma::cor(Genotype,Expression));
}

//[[Rcpp::export]]
arma::uvec isCis(const arma::ivec snp_chrom,const arma::ivec snp_pos, const arma::sword exp_chrom,const arma::sword exp_start, const arma::sword exp_stop,const arma::uword cisdist_cutoff){
  arma::uvec chromvec(snp_chrom.n_elem,arma::fill::zeros);
  chromvec.elem(arma::find(snp_chrom==exp_chrom)).ones();
  if(sum(chromvec)==0){
    return(chromvec);
  }
  arma::ivec distvec=arma::min(arma::join_horiz(arma::abs(snp_pos-exp_start),arma::abs(snp_pos-exp_stop)), 1);
  distvec.elem(arma::find((exp_start-snp_pos)%(snp_pos-exp_stop)>0)).zeros();
  arma::uvec retvec(snp_chrom.n_elem,arma::fill::zeros);
  retvec.elem(arma::find(((distvec<cisdist_cutoff)+(chromvec==1))==2)).ones();
  return(retvec);
}


//[[Rcpp::export]]
arma::uvec addLD(const arma::uvec &snpind, const arma::fmat &LDmat, const float LDcutoff){
  arma::umat posadd=arma::trans(arma::ind2sub(arma::size(LDmat.n_rows,snpind.n_elem),arma::find(abs(LDmat.cols(snpind)>LDcutoff))));
  //      Rcout<<"found "<<posadd.n_cols<<"additional  eQTL from LD gene"<<std::endl;
return(arma::unique(arma::join_vert(snpind,posadd.col(0))));
}


//[[Rcpp::export]]
Rcpp::DataFrame extract_stats(const arma::fmat &Genotype,const Rcpp::DataFrame snpanno,const arma::fmat &Expression, const Rcpp::DataFrame expanno,const arma::fmat &LDmat,const arma::fmat &rmat,const double tcutoff,const double LDcutoff,const arma::uword cisdist,const bool display_progress,bool doCis){
  using namespace Rcpp;
  double n =Genotype.n_rows;
  double rcutoff =tcutoff/sqrt(n+tcutoff*tcutoff-2);
  std::vector<float> betavec;
  std::vector<float> serrvec;
  std::vector<int> snpvec;
  std::vector<int> genevec;
  std::vector<int> cistransvec;
  Rcpp::Rcout<<"Reading SNP annotations"<<std::endl;
  arma::uvec rsids = arma::conv_to<arma::uvec>::from(Rcpp::as<std::vector<int>>(snpanno["rsid"]));
  arma::ivec snp_chrom = Rcpp::as<arma::ivec>(snpanno["chrom"]);
  arma::ivec snp_pos = Rcpp::as<arma::ivec>(snpanno["pos"]);

  Rcpp::Rcout<<"Reading Exp annotations"<<std::endl;
  arma::uvec fgeneids = arma::conv_to<arma::uvec>::from(Rcpp::as<std::vector<int>>(expanno["fgeneid"]));
  arma::ivec exp_chrom = Rcpp::as<arma::ivec>(expanno["chrom"]);
  arma::ivec exp_start = Rcpp::as<arma::ivec>(expanno["TSStart"]);
  arma::ivec exp_stop = Rcpp::as<arma::ivec>(expanno["TSStart"]);

  Rcpp::Rcout<<"Reserving Memory"<<std::endl;
  betavec.reserve(rmat.n_elem);
  serrvec.reserve(rmat.n_elem);
  snpvec.reserve(rmat.n_elem);
  genevec.reserve(rmat.n_elem);
  cistransvec.reserve(rmat.n_elem);
  arma::fvec snpsd=arma::trans(arma::stddev(Genotype));
  arma::fvec expsd=arma::trans(arma::stddev(Expression));
  arma::fvec tr(rmat.n_rows);
  Progress p(rmat.n_cols, display_progress);
  bool skipSearch=false;
  if(doCis){
    std::vector<arma::uword> aschroms=arma::conv_to<std::vector<arma::uword>>::from(arma::unique(snp_chrom));
    std::vector<arma::uword> agchroms=arma::conv_to<std::vector<arma::uword>>::from(arma::unique(exp_chrom));
    std::vector<arma::uword> bchroms;
    std::set_intersection(aschroms.begin(),aschroms.end(),agchroms.begin(),agchroms.end(),std::back_inserter(bchroms));
    if(bchroms.size()==0){
      skipSearch=true;
    }
  }
  if(!skipSearch){
    for(arma::uword i=0; i<rmat.n_cols;i++){
      p.increment();
      tr=rmat.col(i);
      arma::uvec snpind =arma::find(abs(tr)>rcutoff);
      arma::uvec cid=isCis(snp_chrom.elem(snpind),snp_pos.elem(snpind),exp_chrom[i],exp_start[i],exp_stop[i],cisdist);
      if(cid.n_elem!=snpind.n_elem){
        Rcpp::Rcerr<<"differing number of elements in cid and snpind"<<std::endl;
        Rcpp::stop("error in finding cis SNPs");
      }
      if(doCis){
        snpind = snpind.elem(arma::find(cid==1));
      }else{
        snpind = snpind.elem(arma::find(cid==0));
      }

      if(snpind.n_elem>0){
        arma::uvec allind=addLD(snpind,LDmat,LDcutoff);
        //      Rcout<<"found "<<snpind.n_elem<<" sig eQTL for gene"<<i<<std::endl;
        //      Rcout<<"Total of "<<snpind.n_elem<<" eQTL will be recorded for gene "<<i<<", the largest of which is at position: "<<arma::max(snpind)<<std::endl;
        arma::fvec Betas=tr.elem(allind)%(expsd(i)/snpsd.elem(allind));
        arma::fvec tv= sqrt(n-2)*(tr.elem(allind)/arma::sqrt(1-arma::pow(tr.elem(allind),2)));
        arma::fvec serrv=Betas/tv;
        snpind=rsids.elem(allind);
        betavec.insert(betavec.end(),Betas.begin(),Betas.end());
        serrvec.insert(serrvec.end(),serrv.begin(),serrv.end());
        snpvec.insert(snpvec.end(),snpind.begin(),snpind.end());
        std::fill_n(std::back_inserter(genevec),snpind.n_elem,fgeneids(i));
        if(doCis){
          std::fill_n(std::back_inserter(cistransvec),snpind.n_elem,1);
        }else{
          std::fill_n(std::back_inserter(cistransvec),snpind.n_elem,0);
        }
      }
    }
    arma::uvec sizes={betavec.size(),serrvec.size(),snpvec.size(),genevec.size(),cistransvec.size()};
    if(arma::any(sizes!=betavec.size())){
      Rcerr<<"sizes not equal"<<std::endl;
      sizes.print();
      stop("Error: cannot create dataframe from vectors of unequal size");
    }
  }
  NumericVector theta(betavec.begin(),betavec.end());
  NumericVector serr(serrvec.begin(),serrvec.end());
  IntegerVector snpid(snpvec.begin(),snpvec.end());
  IntegerVector geneid(genevec.begin(),genevec.end());
  IntegerVector cistrans(cistransvec.begin(),cistransvec.end());

  return(DataFrame::create(_["theta"]= theta,
                           _["serr"]= serr,
                           _["rsid"]=snpid,
                           _["fgeneid"]=geneid,
                           _["cistrans"]=cistrans));
}









//[[Rcpp::export]]
arma::mat serrMatrix(const arma::mat &Genotype, const arma::mat &Expression, const arma::mat &Betas){
  arma::vec ssqm=arma::trans(arma::sum(arma::pow(Genotype,2))*(double)Genotype.n_rows);
  arma::mat serr(Betas.n_rows,Betas.n_cols);
  arma::mat Ehat(Genotype.n_rows,Genotype.n_cols);
  for(size_t i=0;i<Betas.n_cols;i++){
    Ehat=Genotype*arma::diagmat(Betas.col(i));
    serr.col(i)=arma::trans(arma::sum(arma::pow(Ehat.each_col()-Expression.col(i),2)));
  }
  serr.each_row() /=ssqm;
  return(serr);
}


//[[Rcpp::export]]
void orthogonalize_dataset(std::string h5filename,std::string datagroup, std::string datasetname,std::string newdatasetname,size_t chunksize,const unsigned int deflate_level){
  arma::mat ocovariates= read_dmat_h5(h5filename,"Covardat","covariates",0,40);
  ocovariates = join_horiz(arma::mat(ocovariates.n_rows,1,arma::fill::ones),ocovariates);
  arma::mat Covariates = orthogonalize_covar(ocovariates);
  size_t Nrows= get_rownum_h5(h5filename,datagroup,datasetname);
  Rcpp::Rcout<<"total number of rows:"<<Nrows<<std::endl;
  size_t nchunks=ceil(((double)Nrows)/(double)chunksize);

  Rcpp::Rcout<<"total number of chunks:"<<nchunks<<std::endl;
  for(size_t i=0;i<nchunks;i++){
    Rcpp::Rcout<<"Starting chunk"<<i<<" of size:"<<chunksize<<std::endl;
    size_t offset =i*chunksize;
    Rcpp::Rcout<<"Reading data for chunk"<<i<<std::endl;
    arma::mat Data = read_dmat_h5(h5filename,datagroup,datasetname,offset,chunksize);
    Rcpp::Rcout<<"Orthogonalizing data wrt (orthogonalized)covariates for chunk"<<i<<std::endl;
    arma::mat oData =orthogonalize_data(Data,Covariates);
    Rcpp::Rcout<<"oData is of dimensions:"<<oData.n_rows<<"x"<<oData.n_cols<<std::endl;
    write_mat_h5(h5filename,datagroup,newdatasetname,Nrows,oData.n_rows,oData,deflate_level);
  }

}





