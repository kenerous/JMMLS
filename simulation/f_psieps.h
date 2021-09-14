#ifndef _F_PSIEPS_H
#define _F_PSIEPS_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector f_psieps(NumericVector y, NumericVector xi, NumericVector eta, NumericVector b, NumericVector b_label, NumericVector Queue, NumericVector latent_label, NumericVector psi_eps, int N, int SNT, int NY, double alpha_eps0, double beta_eps0) {
  NumericVector psieps_out = psi_eps;
  int cal = -2; // location of latent variables
  
  for(int i=0; i<NY; ++i){
    double sum = 0;
    if(b_label[i]==0){
      cal +=1;
    }
    if(cal==(-1)){// psi_eps corresponding to eta
      for(int j=0; j<SNT; ++j){
        sum += pow(y[i*SNT+j]-b[i]*eta[j],2);
      }
      psieps_out[i] = 1/(rgamma(1,alpha_eps0+0.5*SNT,1/(beta_eps0+0.5*sum))[0]);
    }else{// psi_eps corresponding to xi
      for(int j=0; j<N; ++j){
        sum += pow(y[i*SNT+Queue[j]]-b[i]*xi[cal*N+j],2);
      }
      psieps_out[i] = 1/(rgamma(1,alpha_eps0+0.5*N,1/(beta_eps0+0.5*sum))[0]);        
    }

  }
  //return List::create(Rcpp::Named("psi_eps") = psieps_out);  
  return psieps_out;
}

#endif
