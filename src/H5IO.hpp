#ifndef H5IO_HPP
#define H5IO_HPP
#include "RcppArmadillo.h"
#include "tbb/tbb.h"
#include "tbb/concurrent_hash_map.h"
#include "tbb/blocked_range.h"
#include "tbb/parallel_for.h"
#include <unordered_map>
#include <vector>
#include <tuple>
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
arma::fmat read_fmat_rowname(const std::string h5file,const std::string annogroupname, const std::string annocolname,const std::string datagroupname, const std::string datacolname,arma::uvec queryvec);


void ip_dist(const arma::frowvec &cummapa,const arma::frowvec &cummapb,arma::fmat &distmat,bool isDiag);
void tbb_dist(const arma::frowvec &cummapa,const arma::frowvec &cummapb,arma::fmat &distmat,bool isDiag);

arma::fmat ip_cov(const arma::fmat &Hpanela, const arma::fmat &Hpanelb, bool isDiag);

arma::mat flip_hap(const std::string hap_h5file,arma::uvec index, const::arma::uword chunk, const arma::uword chunksize,const arma::uword nSNPs);

arma::Mat<int> read_haplotype_h5(const std::string hap_h5file,const size_t readSNPs,const size_t skipSNPs=0);
arma::uvec make_long(arma::uvec &vchrom, arma::uvec &vpos,arma::uvec &posmap);


//
// class Legend{
// private:
//   arma::uvec chromposmap; //vector of length 22 giving cumulative position of last position of the chromomosme (from dbSNP)
//   size_t n;
//   std::unordered_map<arma::uword,arma::uword> absposmap;
// public:
//   arma::uvec chrom;
//   arma::uvec pos;
//   arma::uvec abspos;
//   Legend();
//   Legend(const arma::uvec tchrom,const arma::uvec tpos, const arma::uvec chrompos);
//   Legend subset_index(const arma::uvec &index);
//   Legend subset_chrom(const arma::uword chrom);
//   Legend subset_abspos(const arma::uvec &absp);
//   size_t size()const{return(n);}
// };
//
// struct HashCompare{
//   static size_t hash(const arma::uword& x){
//     return(x);
//   }
//   static bool equal(const arma::uword& x, const arma::uword& y){
//     return(x==y);
//   }
// };
//
//
// class GWAS{
// private:
//   size_t n;
// public:
//   Legend legend;
//   arma::vec betahat;
//   arma::vec serr;
//   arma::uvec rsid;
//   GWAS(const std::string gwash5,const std::string dbsnph5);
// };
//
//
// class GTEx{
// public:
//
// };

//
// typedef tbb::concurrent_hash_map<arma::uword,arma::uword,HashCompare> lposmap;
// using namespace tbb;
//
// class indmap{
//   lposmap lpm;
//
//   rshash(const arma::uvec trsidvec,const arma::uvec tlposvec):rsidvec(trsidvec),lposvec(tlposvec){}
//   void operator()(const blocked_range<arma::uword> &r) const{
//     for(arma::uword i=r.begin(); i!=r.end();++i){
//       rsidmap::accessor a;
//       arma::uword tind=lposvec[i];
//       rstlpos.insert(a,tind);
//       a->second =rsidvec[i];
//     }
//   }
// };

// class dbsnp{
// private:
//   size_t n;
// public:
//   dbsnp(std::string h5file);
//   std::unordered_map<arma::uword,arma::uword> pos2rsidmap;
//   std::unordered_map<arma::uword,arma::uword> rsid2posmap;
//   arma::uvec chromposmap;
//   Legend legend;
//   arma::uvec rsid;
//   Legend find_legend(const arma::uvec &rsidind);
//   arma::umat find_rsid(const Legend &queryind);
//
//   }
// };




#endif
