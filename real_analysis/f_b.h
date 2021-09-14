#ifndef _F_B_H
#define _F_B_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector f_b(NumericVector u, NumericVector alpha, NumericVector B, NumericVector y, NumericVector b_label, NumericVector xi, NumericVector eta, NumericVector b, NumericVector Queue, NumericVector latent_label, NumericVector psi_eps, int N, int SNT, int NY, double b0, double h0) {
  NumericVector b_out = b;
  int cal = -2; // location of latent variables
  
  for(int i=0; i<NY; ++i){
    double mean_b = b0/(psi_eps[i]*h0);
    double sigma_b = 1/(psi_eps[i]*h0);
    if(b_label[i]==0){
      cal +=1;
    }else{
      if(cal==(-1)){// b corresponding to eta
        for(int j=0; j<SNT; ++j){
          sigma_b += pow(eta[j],2)/psi_eps[i];
          mean_b += y[i*SNT+j]*eta[j]/psi_eps[i];
        }
        sigma_b = 1/sigma_b;
        mean_b = mean_b*sigma_b;
        b_out[i] = rnorm(1,mean_b, sqrt(sigma_b))[0];
      }else{
        // b corresponding to xi
        for(int j=0; j<N; ++j){
          sigma_b += pow(xi[cal*N+j],2)/psi_eps[i];
          mean_b += y[i*SNT+Queue[j]]*xi[cal*N+j]/psi_eps[i];
        }
        sigma_b = 1/sigma_b;
        mean_b = mean_b*sigma_b;
        b_out[i] = rnorm(1,mean_b, sqrt(sigma_b))[0];        
      }
    }
  }
  return b_out;  
}

#endif
