#ifndef _F_ETA_H
#define _F_ETA_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector f_eta(NumericVector x, NumericVector alpha, int NX, NumericVector B, NumericVector u, NumericVector NT, NumericVector Queue, int N, int NU, int SNT) {
  NumericVector eta_out(SNT);
  for(int i=0; i<N; ++i){
  	int t_temp = NT[i];
  	NumericVector second_level(NU);
  	for(int j=0; j<NU; ++j){
  		for(int k=0; k<(NX+1); ++k){
  			second_level[j] += alpha[j*(NX+1)+k] * x[k*N+i];
		  }
	  }
  	for(int j=0; j<t_temp; ++j){  		
  		for(int k=0; k<NU; ++k){
  			eta_out[Queue[i]+j] += (second_level[k] + u[k*N+i])*B[k*SNT+Queue[i]+j]; 
	  	}
	  }

  }
  return eta_out;  
}

#endif
