#ifndef _F_C_H
#define _F_C_H

#include "f_mmult.h"
#include "f_integral.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List f_c(int IF_QUADRATIC, int NX, NumericVector x, NumericVector alpha, NumericVector Queue, int SNT, NumericVector OT_label, NumericVector B, NumericVector kns, double c0, double sigma_c0, NumericVector xi, NumericVector u, NumericVector DELTA, NumericVector z, NumericVector OT, NumericVector Fail_label, NumericVector J, int at_c, NumericVector gamma, NumericVector beta, NumericVector lambda, int N, int NU, int NXI, int NZ, double c, int K, double c_c){
  double c_out = c;
  NumericVector u1(NU);
  NumericVector x1(NX+1);
  NumericMatrix beta1(1,NZ);
  NumericMatrix z1(NZ,1);
  NumericMatrix gamma1(1,NXI);
  NumericMatrix xi1(NXI,1);
  
  for(int i=0; i<NZ; ++i){
    beta1(0,i) = beta[i];
  }
  for(int i=0; i<NXI; ++i){
    gamma1(0,i) = gamma[i];
  }
  
  double c_star = rnorm(1,c,c_c)[0];
  
  double ratio_c = pow(c_star-c0,2)/(-2) - pow(c-c0,2)/(-2);
  
  for(int i=0; i<N; ++i){
    double ftu = 0;
    NumericVector second_level(NU);
    for(int j=0; j<(NX+1); ++j){
    	x1[j] = x[j*N+i];
	}
    for(int j=0; j<NU; ++j){
    	for(int k=0; k<(NX+1); ++k){
    		second_level[j] += alpha[j*(NX+1)+k]*x1[k];
		}
      u1[j] = u[j*N+i];
      //ftu += pow(OT[i],j)*u1[j];
      ftu += B[j*SNT+Queue[i]+OT_label[i]]*(u1[j]+second_level[j]);
    }
    for(int j=0; j<NZ; ++j){
      z1(j,0) = z[j*N+i];
    }
    for(int j=0; j<NXI; ++j){
      xi1(j,0) = xi[j*N+i];
    }
    
    
    if(DELTA[i]==1){
      ratio_c += (c_star-c)*ftu;
    }
    
    ratio_c += (f_integral(IF_QUADRATIC, NX, x1, alpha, kns, u1, J, lambda, Fail_label[i], K, NU, OT[i], c, 0, Fail_label[i])-f_integral(IF_QUADRATIC, NX, x1, alpha, kns, u1, J, lambda, Fail_label[i], K, NU, OT[i], c_star, 0, Fail_label[i]))*exp(f_mmult(beta1,z1)(0,0)+f_mmult(gamma1,xi1)(0,0));
    
  }
  
  ratio_c = exp(ratio_c);
  
  double randnum = runif(1)[0];
  if(randnum<ratio_c){
    c_out = c_star;
    at_c += 1;
  }
  
  
  
  return List::create(Rcpp::Named("c") = c_out,
                      Rcpp::Named("ratio_c") = ratio_c,
                      Rcpp::Named("at_c") = at_c);
}
#endif
