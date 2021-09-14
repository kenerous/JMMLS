#ifndef _F_PSIDEL_H
#define _F_PSIDEL_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double f_psidel(int NX, NumericVector x, NumericVector alpha, NumericVector B, NumericVector eta, NumericVector u, NumericVector Time, NumericVector Queue, NumericVector NT, int N, int SNT, int NU, double alpha_del0, double beta_del0, double psi_del) {
  double psidel_out = psi_del;
  NumericVector ft(SNT);
  NumericVector sum1(SNT);
  for(int i=0; i<N; ++i){
    int t_temp = NT[i];
    NumericVector second_level(NU);
    for(int j=0; j<NU; ++j){
    	for(int l=0 ;l<(NX+1); ++l){
      		second_level[j] += alpha[j*(NX+1)+l]*x[l*N+i];
		}
	}
    for(int j=0; j<t_temp; ++j){    	
      for(int k=0; k<NU; ++k){
        ft[Queue[i]+j] += B[k*SNT+Queue[i]+j]*(u[k*N+i]+second_level[k]);
      }
      
      sum1[Queue[i]+j] = pow(eta[Queue[i]+j]-ft[Queue[i]+j],2);
    }
  }
  psidel_out = 1/(rgamma(1,alpha_del0+0.5*SNT,1/(beta_del0+0.5*sum(sum1)))[0]);

  return psidel_out;  
}
#endif
