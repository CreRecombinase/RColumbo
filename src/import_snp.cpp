#include "RcppArmadillo.h"
//[[Rcpp::depends(BH,RcppArmadillo)]]
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/qi_match.hpp>
#include <boost/spirit/include/qi_uint.hpp>
#include <iterator>
#include <unordered_map>
#include <numeric>
#include <fstream>
#include <string>
#include <H5Cpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include "h5func.hpp"

using namespace H5;


size_t makeFlip(arma::uvec &isFlipv, arma::fmat &genomat){
  size_t ct;
  if(isFlipv.n_elem!=genomat.n_cols){
    Rcpp::Rcerr<<"isFlipv has:"<<isFlipv.n_elem<<" and genomat has colnum:"<<genomat.n_cols<<std::endl;
    Rcpp::stop("error in makeFlip!");
  }
  for(size_t i=0; i<isFlipv.n_elem;i++){
    if(isFlipv[i]==1){
      genomat.col(i)=arma::abs(2-genomat.col(i));
      ct++;
    }
  }
  return(ct);
}



arma::uvec isFlip(std::vector<std::string> &ref,std::vector<std::string> &alt){
  if(ref.size()!=alt.size()){
    Rcpp::Rcout<<"ref is of size:"<<ref.size()<<" alt is of size:"<<alt.size()<<std::endl;
    Rcpp::stop("size of ref not equal to size of alt!");
  }
  arma::uvec ret(ref.size());
  for(size_t i=0;i<ret.size();i++){
    if(ref[i]<alt[i]){
      ret[i]=1;
    }else{
      ret[i]=0;
    }
  }
  return(ret);
}





size_t read_genotype_gz(boost::iostreams::filtering_istream &fs, const size_t Nsnps, const size_t Nind,const size_t chunksize,arma::uvec &mchroms,arma::uvec &mposs, std::vector<std::string> &refs, std::vector<std::string> &alts,arma::fmat &mgenotypes){
  using namespace Rcpp;
  using namespace boost::spirit::qi;
  //  std::vector<double> genotypes;

  std::vector<int> chroms;
  std::vector<unsigned int> poss;
  std::vector<float> genotypes;

  chroms.reserve(chunksize);
  poss.reserve(chunksize);
  genotypes.reserve(chunksize*Nind);

  refs.clear();
  refs.reserve(chunksize);
  alts.clear();
  alts.reserve(chunksize);
  size_t ct=0;
  std::string line;
  while(getline(fs,line)){

    std::string::const_iterator sbeg = line.begin();
    std::string::const_iterator send = line.end();
    phrase_parse(sbeg,send,int_>>'_'>>uint_>>'_'>>as_string[+char_("ACTGN")]>>"_">>as_string[+char_("ACTGN")]>>"_">>"b37">>*float_,space,chroms,poss,refs,alts,genotypes);

    if(refs.size()!=alts.size()){
      Rcpp::Rcerr<<"refs and alts different sizes! at line "<<ct<<" ("<<refs.size()<<" "<<alts.size()<<")"<<std::endl;
      Rcpp::stop("error in read_genotype_gz");
    }
    ct++;
    if(ct==chunksize){
      break;
    }
  }

  mchroms=arma::conv_to<arma::uvec>::from(chroms);
  mposs = arma::conv_to<arma::uvec>::from(poss);
  mgenotypes = arma::fmat(&genotypes[0],Nind,ct);
  arma::uvec sizes={mchroms.n_elem,mposs.n_elem,mgenotypes.n_cols,(arma::uword)refs.size(),(arma::uword)alts.size()};
  if(any(sizes!=mchroms.n_elem)){
    sizes.print();
    Rcpp::Rcerr<<"not all sizes are equal!:"<<std::endl;
    Rcpp::stop("error in read_genotype_gz");
  }

  return(ct);
}

// [[Rcpp::export]]
size_t lineno_gz(const char* filename){
  using namespace Rcpp;
  Rcout<<"Counting line number"<<std::endl;
  std::ifstream file(filename,std::ios_base::in | std::ios_base::binary);
  size_t ln=0;
  try{
    boost::iostreams::filtering_istream fs;
    fs.push(boost::iostreams::gzip_decompressor());
    fs.push(file);
    std::string tline;
    while(getline(fs,tline)){
      ln++;
    }
  }
  catch(const boost::iostreams::gzip_error& e){
    Rcerr<<e.what()<<std::endl;
  }
  return(ln-1);
}


// [[Rcpp::export]]
int write_genotype_h5(const char* snpdatmat,size_t Nind,size_t chunksize, const std::string h5file, bool doFlip,const unsigned int deflate_level){
  using namespace Rcpp;
  size_t Nsnps = lineno_gz(snpdatmat);
  std::cout<<"Starting to map file"<<std::endl;
  boost::iostreams::mapped_file_source mapfile(snpdatmat);
  boost::iostreams::stream<boost::iostreams::mapped_file_source> textstream(mapfile);
  boost::iostreams::filtering_istream fs;
  fs.push(boost::iostreams::gzip_decompressor{});
  fs.push(textstream);
  std::string title;
  getline(fs,title);
  size_t sr=0;
  size_t scum=0;
  std::cout<<"Starting to read genotype data"<<std::endl;

  arma::uvec chroms;
  arma::uvec poss;
  arma::uvec retdoFlip;
  arma::fmat genodat;
  std::vector<std::string> refs;
  std::vector<std::string> alts;
  size_t count=0;

  while(scum<Nsnps){
    sr = read_genotype_gz(fs, Nsnps, Nind,chunksize,chroms,poss,refs,alts,genodat);

    scum=sr+scum;
    Rcout<<"Line "<<scum<<" of "<<Nsnps<<std::endl;


    size_t retn =genodat.n_cols;
    retdoFlip=isFlip(refs,alts);
    if(doFlip){
      std::cout<<"Flipping alleles"<<std::endl;
      makeFlip(retdoFlip,genodat);
    }

  //  std::cout<<"Writing genotype matrix"<<std::endl;
    write_mat_h5(h5file, "SNPdata", "genotype", Nsnps,Nind,genodat,deflate_level);
//    std::cout<<"Writing retchrom"<<std::endl;
    write_int_h5(h5file,"SNPinfo","chrom",Nsnps,chroms,deflate_level);
//    std::cout<<"Writing retpos "<<std::endl;
    write_uint_h5(h5file,"SNPinfo","pos",Nsnps,poss,deflate_level);
//    std::cout<<"Writing doFlip"<<std::endl;
    write_int_h5(h5file,"SNPinfo","doFlip",Nsnps,retdoFlip,deflate_level);

    if(doFlip){
      std::cout<<"Writing doFlip"<<std::endl;
    }
    count++;
  }
  return(count);
}


