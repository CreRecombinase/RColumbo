#ifndef H5IO_HPP
#define H5IO_HPP
#include "RcppArmadillo.h"
#include <vector>
#include <memory>
#include <H5Cpp.h>

arma::mat convertTSparse(Rcpp::S4 &mat);

int findbandwidth(Rcpp::IntegerVector &i, Rcpp::IntegerVector &j, Rcpp::NumericVector &x,double cutoff);

int findcutoff(Rcpp::IntegerVector &i, Rcpp::IntegerVector &j, Rcpp::NumericVector &x,int bandwidth);

size_t write_haplotype_h5(const std::string hap_gzfile,const std::string hap_h5file,const size_t nrows,const size_t ncols,size_t chunksize,const unsigned int deflate_level);

arma::uvec read_flip(const std::string hap_h5file,arma::uvec indexes);

arma::Mat<int> read_haplotype_ind_h5(const std::string hap_h5file,arma::uvec indexes);

arma::mat read_dmat_ind_h5(const std::string hap_h5file,const std::string groupname, const std::string dataname,arma::uvec indexes);

arma::uvec ind_lookup(const arma::uvec &queryvec, const arma::uvec &targetvec);

arma::mat read_dmat_rowname(const std::string h5file,const std::string annogroupname, const std::string annocolname,const std::string datagroupname, const std::string datacolname,arma::uvec queryvec);

void ip_dist(const arma::frowvec &cummapa,const arma::frowvec &cummapb,arma::fmat &distmat,bool isDiag);

arma::fmat ip_cov(const arma::fmat &Hpanela, const arma::fmat &Hpanelb, bool isDiag);

arma::mat flip_hap(const std::string hap_h5file,arma::uvec index, const::arma::uword chunk, const arma::uword chunksize,const arma::uword nSNPs);

arma::Mat<int> read_haplotype_h5(const std::string hap_h5file,const size_t readSNPs,const size_t skipSNPs=0);

#endif
