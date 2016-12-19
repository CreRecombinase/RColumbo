#include "RcppArmadillo.h"
#include "h5func.hpp"
#include <cmath>
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
Rcpp::DataFrame cor_h5(const std::string h5file, const std::string groupname, const std::string dataname, const arma::uvec indvec,const float LDcutoff,const bool cutBelow){
  using namespace Rcpp;
  arma::fmat tdata= read_fmat_chunk_ind(h5file,groupname,dataname,indvec);
  arma::fmat rmat=arma::abs(cor(tdata));
  arma::umat sigind;
  if(cutBelow){
    sigind=arma::ind2sub(arma::size(rmat),arma::find(rmat<LDcutoff));
  }else{
    sigind=arma::ind2sub(arma::size(rmat),arma::find(rmat>LDcutoff));
  }
  arma::uvec rowind=indvec.elem(arma::trans(sigind.row(0)));
  arma::uvec colind=indvec.elem(arma::trans(sigind.row(1)));
  Rcpp::IntegerVector Row(rowind.begin(),rowind.end());
  Rcpp::IntegerVector Col(colind.begin(),colind.end());
  return(Rcpp::DataFrame::create(_["rowind"]= Row,
                                 _["colind"]= Col));
}


//[[Rcpp::export]]
arma::fmat rMatrix(const arma::fmat &Genotype,const arma::fmat &Expression){
  return(arma::cor(Genotype,Expression));
}

//[[Rcpp::export]]
arma::uvec isCis_mat(const arma::ivec snp_chrom,const arma::ivec snp_pos, const arma::ivec exp_chrom,const arma::ivec exp_start, const arma::ivec exp_stop,const arma::uword cisdist_cutoff){

  arma::uvec retvec(snp_chrom.n_elem,arma::fill::zeros);
  arma::uvec cisc=arma::find(snp_chrom==exp_chrom);
  if(cisc.n_elem==0){
    return(retvec);
  }
  arma::ivec distvec=arma::min(arma::join_horiz(arma::abs(snp_pos.elem(cisc)-exp_start.elem(cisc)),arma::abs(snp_pos.elem(cisc)-exp_stop.elem(cisc))),1);
  cisc=cisc.elem(arma::find(distvec<cisdist_cutoff));
  if(cisc.n_elem==0){
    return(retvec);
  }
  retvec.elem(cisc).ones();
  return(retvec);
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
  //Add SNPs that are adjacent to significant SNPs due to being in LD with them
  arma::umat posadd=arma::trans(arma::ind2sub(arma::size(LDmat.n_rows,snpind.n_elem),arma::find(abs(LDmat.cols(snpind)>LDcutoff))));
  //      Rcout<<"found "<<posadd.n_cols<<"additional  eQTL from LD gene"<<std::endl;
return(arma::unique(arma::join_vert(snpind,posadd.col(0))));
}


arma::uvec chunk_index(const arma::uvec indexvec, const size_t chunksize,size_t i){
  size_t nchunks =ceil(((double)indexvec.n_elem)/(double)chunksize);
  if(i>nchunks){
    Rcpp::stop("i greater than nchunks!");
  }
    size_t istart=i*chunksize;
    size_t istop=std::min((i+1)*chunksize-1,(size_t)indexvec.n_elem-1);
   return(indexvec(arma::span(istart,istop)));
}

//[[Rcpp::export]]
Rcpp::DataFrame fast_eQTL(const arma::fmat &Genotype, const Rcpp::DataFrame snpanno, const arma::fmat &Expression, const Rcpp::DataFrame expanno, const double cis_tcutoff, const double trans_tcutoff,const arma::uword cisdist, const bool doTrans,const bool doCis){
  using namespace Rcpp;
  double n =Genotype.n_rows;
  double cis_rcutoff =cis_tcutoff/sqrt(n+cis_tcutoff*cis_tcutoff-2);
  double trans_rcutoff =trans_tcutoff/sqrt(n+trans_tcutoff*trans_tcutoff-2);

  if(!doCis&&!doTrans){
    Rcpp::stop("at least one of doCis or doTrans must be true");
  }
  std::vector<int> cistransvec;
  Rcpp::Rcout<<"Reading SNP annotations"<<std::endl;
  // arma::uvec rsids = arma::conv_to<arma::uvec>::from(Rcpp::as<std::vector<int>>(snpanno["rsid"]));
  arma::ivec snp_chrom = Rcpp::as<arma::ivec>(snpanno["chrom"]);
  arma::ivec snp_pos = Rcpp::as<arma::ivec>(snpanno["pos"]);

  Rcpp::Rcout<<"Reading Exp annotations"<<std::endl;
  arma::uvec fgeneids = arma::conv_to<arma::uvec>::from(Rcpp::as<std::vector<int>>(expanno["fgeneid"]));
  arma::ivec exp_chrom = Rcpp::as<arma::ivec>(expanno["chrom"]);
  arma::ivec exp_start = Rcpp::as<arma::ivec>(expanno["start"]);
  arma::ivec exp_stop = Rcpp::as<arma::ivec>(expanno["end"]);
  Rcpp::Rcout<<"Computing correlation"<<std::endl;
  arma::fmat rmat =rMatrix(Genotype,Expression);

  Rcpp::Rcout<<"Reserving Memory"<<std::endl;

  // betavec.reserve(rmat.n_elem);
  // serrvec.reserve(rmat.n_elem);
  // snppvec.reserve(rmat.n_elem);
  // snpcvec.reserve(rmat.n_elem);
  // genevec.reserve(rmat.n_elem);
  // cistransvec.reserve(rmat.n_elem);
  arma::fvec Betas;
  arma::fvec serrv;
  arma::uvec cid;
  Rcpp::Rcout<<"Computing std error"<<std::endl;
  arma::fvec snpsd=arma::trans(arma::stddev(Genotype));
  arma::fvec expsd=arma::trans(arma::stddev(Expression));
  double rcutoff;
  arma::fvec rcutoff_vec(2);
  rcutoff_vec[0]=trans_rcutoff;
  rcutoff_vec[1]=cis_rcutoff;

  bool skipSearch=false;

  if(!doTrans){
    std::vector<arma::uword> aschroms=arma::conv_to<std::vector<arma::uword>>::from(arma::unique(snp_chrom));
    std::vector<arma::uword> agchroms=arma::conv_to<std::vector<arma::uword>>::from(arma::unique(exp_chrom));
    std::vector<arma::uword> bchroms;
    std::set_intersection(aschroms.begin(),aschroms.end(),agchroms.begin(),agchroms.end(),std::back_inserter(bchroms));
    if(bchroms.size()==0){
      skipSearch=true;
    }
  }
  if(!skipSearch){
    Rcpp::Rcout<<"Mapping eQTL"<<std::endl;
    double rcutoff=std::min(cis_rcutoff,trans_rcutoff);
    arma::uvec sigr=arma::find(abs(rmat)>rcutoff);
    arma::fvec tr = rmat.elem(sigr);
    arma::umat sigmat= arma::ind2sub(arma::size(rmat),sigr);
    rmat.clear();
    arma::uvec snpind = arma::trans(sigmat.row(0));
    arma::uvec geneind = arma::trans(sigmat.row(1));
    Rcpp::Rcout<<"ID-ing cis-trans relationships"<<std::endl;
    cid=isCis_mat(snp_chrom.elem(snpind),snp_pos.elem(snpind),exp_chrom.elem(geneind),exp_start.elem(geneind),exp_stop.elem(geneind),cisdist);
    Rcpp::Rcout<<"Subsetting Annotations"<<std::endl;
    if(!doCis){
      arma::uvec tid=arma::find(cid==0);
      snpind = snpind.elem(tid);
      geneind = geneind.elem(tid);
      cid = cid.elem(tid);
      tr=tr.elem(tid);
    }
    if(!doTrans){
      arma::uvec tid=arma::find(cid==1);
      snpind = snpind.elem(tid);
      geneind = geneind.elem(tid);
      cid = cid.elem(tid);
      tr=tr.elem(tid);
    }
    arma::uvec sigid=arma::find(abs(tr)>rcutoff_vec.elem(cid));
    snpind=snpind.elem(sigid);
    geneind=geneind.elem(sigid);
    cid = cid.elem(sigid);
    tr=tr.elem(sigid);
    if(snpind.n_elem>0){
      //        arma::uvec allind=addLD(snpind,LDmat,LDcutoff);
      //      Rcout<<"found "<<snpind.n_elem<<" sig eQTL for gene"<<i<<std::endl;
      //      Rcout<<"Total of "<<snpind.n_elem<<" eQTL will be recorded for gene "<<i<<", the largest of which is at position: "<<arma::max(snpind)<<std::endl;
      Betas=tr%(expsd(geneind)/snpsd.elem(snpind));
      arma::fvec tv= sqrt(n-2)*(tr/arma::sqrt(1-arma::pow(tr,2)));
      serrv=Betas/tv;
      snp_pos=snp_pos.elem(snpind);
      snp_chrom=snp_chrom.elem(snpind);
      fgeneids=fgeneids.elem(geneind);
    }
    arma::uvec sizes={Betas.n_elem,serrv.n_elem,snp_pos.n_elem,fgeneids.n_elem,cid.n_elem};
    if(arma::any(sizes!=Betas.n_elem)){
      Rcerr<<"sizes not equal"<<std::endl;
      sizes.print();
      stop("Error: cannot create dataframe from vectors of unequal size");
    }
  }
  NumericVector theta(Betas.begin(),Betas.end());
  NumericVector serr(serrv.begin(),serrv.end());
  IntegerVector snpchrom(snp_chrom.begin(),snp_chrom.end());
  IntegerVector snppos(snp_pos.begin(),snp_pos.end());
  IntegerVector geneid(fgeneids.begin(),fgeneids.end());
  IntegerVector cistrans(cid.begin(),cid.end());

  return(DataFrame::create(_["theta"]= theta,
                           _["serr"]= serr,
                           _["chrom"]=snpchrom,
                           _["pos"]=snppos,
                           _["fgeneid"]=geneid,
                           _["cistrans"]=cistrans));
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
  //  Progress p(rmat.n_cols, display_progress);
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





