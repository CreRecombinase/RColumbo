#include "RcppArmadillo.h"
#include "tbb/tbb.h"
#include <vector>
using namespace tbb;
#pragma warning(disable: 588)

//[[Rcpp::export]]
arma::vec pMu(const arma::vec betahat, const arma::vec serr,const double pi,const double tau){
  arma::vec output(betahat.size());
  size_t n=betahat.size();
  parallel_for(size_t(0),n,[&](size_t i){
    output[i]=1/(1+((1-pi)*R::dnorm4(betahat[i],0,serr[i],0))/(pi*R::dnorm4(betahat[i],0,sqrt(tau*tau+serr[i]*serr[i]),0)));
  }
  );
  return(output);
}

//[[Rcpp::export]]
Rcpp::NumericVector sslab_em(const arma::vec p, const arma::vec bh,const arma::vec si){
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

//[[Rcpp::export]]
double sslab_lik(const arma::vec p,const arma::vec bh, const arma::vec si){
  using namespace Rcpp;
  double pi=p[0];
  double tau=p[1];
  arma::vec mu = pMu(bh,si,pi,tau);
  arma::vec iprobs(mu.n_elem);
  for(size_t i=0;i<iprobs.n_elem; i++){
    iprobs[i] =mu[i]*log(R::dnorm4(bh[i],0,sqrt(tau*tau+si[i]*si[i]),0)*pi)+(1-mu[i])*log(R::dnorm4(bh[i],0,si[i],0)*(1-pi));
  }
  return(-arma::sum(iprobs));
}


  // sslab.lik <- function(p,bh,si){
  //   temp_pi <- p[1]
  //   temp_tau <- p[2]
  //   uzin <- temp_pi*dnorm(bh,mean=0,sd=sqrt(temp_tau^2+si^2))
  //   uizd <- uzin+(1-temp_pi)*dnorm(bh,mean=0,sd=si,log=F)
  //   uiz <- uzin/uizd
  //   fiprobs <- uiz*dnorm(bh,mean=0,sd=sqrt(temp_tau^2+si^2))
  //   siprobs <-(1-uiz)*dnorm(bh,mean=0,sd=si)
  //   iprobs <- uiz*log(dnorm(bh,mean=0,sd=sqrt(temp_tau^2+si^2))*temp_pi)+(1-uiz)*log(dnorm(bh,mean=0,sd=si)*(1-temp_pi))
  //   return(-sum(iprobs))
  // }



