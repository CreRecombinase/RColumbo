// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// write_haplotype_h5
void write_haplotype_h5(const std::string hap_gzfile, const std::string hap_h5file, const size_t nrows, const size_t ncols, size_t chunksize);
RcppExport SEXP RColumbo_write_haplotype_h5(SEXP hap_gzfileSEXP, SEXP hap_h5fileSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP, SEXP chunksizeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const std::string >::type hap_gzfile(hap_gzfileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type hap_h5file(hap_h5fileSEXP);
    Rcpp::traits::input_parameter< const size_t >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< size_t >::type chunksize(chunksizeSEXP);
    write_haplotype_h5(hap_gzfile, hap_h5file, nrows, ncols, chunksize);
    return R_NilValue;
END_RCPP
}
// read_haplotype_ind_h5
arma::Mat<int> read_haplotype_ind_h5(const std::string hap_h5file, arma::uvec indexes);
RcppExport SEXP RColumbo_read_haplotype_ind_h5(SEXP hap_h5fileSEXP, SEXP indexesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const std::string >::type hap_h5file(hap_h5fileSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type indexes(indexesSEXP);
    __result = Rcpp::wrap(read_haplotype_ind_h5(hap_h5file, indexes));
    return __result;
END_RCPP
}
// flip_hap
arma::mat flip_hap(const std::string hap_h5file, arma::uvec index, arma::uvec doFlip, const ::arma::uword chunk, const arma::uword chunksize, const arma::uword nSNPs);
RcppExport SEXP RColumbo_flip_hap(SEXP hap_h5fileSEXP, SEXP indexSEXP, SEXP doFlipSEXP, SEXP chunkSEXP, SEXP chunksizeSEXP, SEXP nSNPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const std::string >::type hap_h5file(hap_h5fileSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type index(indexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type doFlip(doFlipSEXP);
    Rcpp::traits::input_parameter< const ::arma::uword >::type chunk(chunkSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type chunksize(chunksizeSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type nSNPs(nSNPsSEXP);
    __result = Rcpp::wrap(flip_hap(hap_h5file, index, doFlip, chunk, chunksize, nSNPs));
    return __result;
END_RCPP
}
// flip_hap_LD
arma::sp_mat flip_hap_LD(const std::string hap_h5file, arma::uvec index, arma::uvec doFlip, arma::rowvec map, const int m, const double Ne, const double cutoff, const arma::uword i, const arma::uword j, const arma::uword chunksize);
RcppExport SEXP RColumbo_flip_hap_LD(SEXP hap_h5fileSEXP, SEXP indexSEXP, SEXP doFlipSEXP, SEXP mapSEXP, SEXP mSEXP, SEXP NeSEXP, SEXP cutoffSEXP, SEXP iSEXP, SEXP jSEXP, SEXP chunksizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const std::string >::type hap_h5file(hap_h5fileSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type index(indexSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type doFlip(doFlipSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type map(mapSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type i(iSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type j(jSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type chunksize(chunksizeSEXP);
    __result = Rcpp::wrap(flip_hap_LD(hap_h5file, index, doFlip, map, m, Ne, cutoff, i, j, chunksize));
    return __result;
END_RCPP
}
// read_haplotype_h5
arma::Mat<int> read_haplotype_h5(const std::string hap_h5file, const size_t readSNPs, const size_t skipSNPs);
RcppExport SEXP RColumbo_read_haplotype_h5(SEXP hap_h5fileSEXP, SEXP readSNPsSEXP, SEXP skipSNPsSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const std::string >::type hap_h5file(hap_h5fileSEXP);
    Rcpp::traits::input_parameter< const size_t >::type readSNPs(readSNPsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type skipSNPs(skipSNPsSEXP);
    __result = Rcpp::wrap(read_haplotype_h5(hap_h5file, readSNPs, skipSNPs));
    return __result;
END_RCPP
}
// gen_rows
arma::uvec gen_rows(const int i, const size_t nrows, const size_t chunksize);
RcppExport SEXP RColumbo_gen_rows(SEXP iSEXP, SEXP nrowsSEXP, SEXP chunksizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const int >::type i(iSEXP);
    Rcpp::traits::input_parameter< const size_t >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type chunksize(chunksizeSEXP);
    __result = Rcpp::wrap(gen_rows(i, nrows, chunksize));
    return __result;
END_RCPP
}
// read_hap_txt
arma::Mat<short> read_hap_txt(const char* inhapfile);
RcppExport SEXP RColumbo_read_hap_txt(SEXP inhapfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const char* >::type inhapfile(inhapfileSEXP);
    __result = Rcpp::wrap(read_hap_txt(inhapfile));
    return __result;
END_RCPP
}
// read_hap_h5
arma::mat read_hap_h5(const char* inhapfile);
RcppExport SEXP RColumbo_read_hap_h5(SEXP inhapfileSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const char* >::type inhapfile(inhapfileSEXP);
    __result = Rcpp::wrap(read_hap_h5(inhapfile));
    return __result;
END_RCPP
}
// p_sparse_LD
arma::sp_mat p_sparse_LD(const arma::rowvec cummap, const arma::mat Hpanel, const double Ne, const int m, const double cutoff, const arma::uword chunksize, arma::uword i, arma::uword j);
RcppExport SEXP RColumbo_p_sparse_LD(SEXP cummapSEXP, SEXP HpanelSEXP, SEXP NeSEXP, SEXP mSEXP, SEXP cutoffSEXP, SEXP chunksizeSEXP, SEXP iSEXP, SEXP jSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::rowvec >::type cummap(cummapSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Hpanel(HpanelSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type chunksize(chunksizeSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type i(iSEXP);
    Rcpp::traits::input_parameter< arma::uword >::type j(jSEXP);
    __result = Rcpp::wrap(p_sparse_LD(cummap, Hpanel, Ne, m, cutoff, chunksize, i, j));
    return __result;
END_RCPP
}
// sparse_LD
arma::sp_mat sparse_LD(const arma::vec cummap, const arma::mat Hpanel, const double Ne, const int m, const double cutoff, const arma::uword report_every);
RcppExport SEXP RColumbo_sparse_LD(SEXP cummapSEXP, SEXP HpanelSEXP, SEXP NeSEXP, SEXP mSEXP, SEXP cutoffSEXP, SEXP report_everySEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< const arma::vec >::type cummap(cummapSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type Hpanel(HpanelSEXP);
    Rcpp::traits::input_parameter< const double >::type Ne(NeSEXP);
    Rcpp::traits::input_parameter< const int >::type m(mSEXP);
    Rcpp::traits::input_parameter< const double >::type cutoff(cutoffSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type report_every(report_everySEXP);
    __result = Rcpp::wrap(sparse_LD(cummap, Hpanel, Ne, m, cutoff, report_every));
    return __result;
END_RCPP
}
