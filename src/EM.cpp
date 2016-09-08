#include "RcppArmadillo.h"
#include "tbb/tbb.h"
#include <vector>
using namespace tbb;
#pragma warning(disable: 588)

//[[Rcpp::export]]
arma::vec pMu(const arma::vec &betahat, const arma::vec &serr,const double pi,const double tau){
  arma::vec output(betahat.size());
  size_t n=betahat.size();
  parallel_for(size_t(0),n,[&](size_t i){
    output[i]=1/(1+((1-pi)*R::dnorm4(betahat[i],0,serr[i],0))/(pi*R::dnorm4(betahat[i],0,sqrt(tau*tau+serr[i]*serr[i]),0)));
                 }
  );
  return(output);
}




//[[Rcpp::export]]
Rcpp::NumericVector sslab_em(const arma::vec &p, const arma::vec bh,const arma::vec &si){
  using namespace Rcpp;
  double pi = p[0];
  double tau = p[1];
  arma::vec mu= pMu(bh,si,pi,tau);
  Rcpp::NumericVector newp(2);
  double smu=arma::sum(mu);
  newp[0]=smu/mu.n_elem;
  newp[1]=sqrt(arma::sum(mu%(bh%bh))/smu-arma::sum((si%si)%mu)/smu);
  return(newp);
}

// sslab.lik <- function(p,y){
//   stopifnot(ceiling(length(y)/2)==floor(length(y)/2))
//   bh <- y[1:(length(y)/2)]
//   si <- y[(length(y)/2+1):length(y)]
//   pi <- p[1]
//   tau <- p[2]
//   uzin <- pi*dnorm(bh,mean=0,sd=sqrt(tau^2+si^2))
//   uizd <- uzin+(1-pi)*dnorm(bh,mean=0,sd=si,log=F)
//   uiz <- uzin/uizd
//   return(-sum(uiz*dnorm(bh,mean=0,sd=sqrt(tau^2+si^2))+(1-uiz)*dnorm(bh,mean=0,sd=si)))
// }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

