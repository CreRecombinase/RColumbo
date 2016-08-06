#ifndef LD_HPP
#define LD_HPP
#include <Rcpp.h>
#include "RcppArmadillo.h"

void p_dist(const arma::rowvec &cummapa,const arma::rowvec &cummapb, arma::mat &distmat,bool isDiag);
void p_cov(const arma::mat &Hpanela, const arma::mat &Hpanelb, arma::mat &covmat, bool isDiag);


#endif
