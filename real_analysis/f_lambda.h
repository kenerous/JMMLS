#ifndef _F_LAMBDA_H
#define _F_LAMBDA_H

#include "f_mmult.h"
#include "f_integral.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List f_lambda(int IF_QUADRATIC, int NX, NumericVector x, NumericVector alpha, NumericVector kns, NumericVector J, NumericVector lambda, NumericVector xi, NumericVector z, NumericVector OT, NumericVector DELTA, NumericVector u, NumericVector Fail_label, NumericVector beta, NumericVector gamma, int N, int NU, int K, int NZ, int NXI, double alpha1, double alpha2, double c) {
  NumericVector lambda_out = lambda;
  int L = J.size()-1; // number of intervals for estimating baseline hazard rate
  NumericVector lambda_null(L);
  NumericMatrix beta1(1,NZ);
  NumericMatrix gamma1(1,NXI);
  NumericVector sum1(L);
  NumericVector sum2(L);
  for(int i=0; i<L; ++i){
    lambda_null[i] = 1;
  }
  for(int i=0; i<NZ; ++i){
    beta1(0,i) = beta[i];
  }
  for(int i=0; i<NXI; ++i){
    gamma1(0,i) = gamma[i];
  }
  
  for(int i=0; i<L; ++i){

    for(int j=0; j<N; ++j){
      NumericVector u1(NU);
      NumericMatrix z1(NZ,1);
      NumericMatrix xi1(NXI,1);
      NumericVector x1(NX+1);
      
      for(int k=0; k<NU; ++k){
        u1[k] = u[k*N+j];
      }
      for(int k=0; k<NZ; ++k){
        z1(k,0) = z[k*N+j];
      }
      for(int k=0; k<NXI; ++k){
        xi1(k,0) = xi[k*N+j];
      }
      for(int k=0; k<(NX+1); ++k){
      	x1[k] = x[k*N+j];
	  }
      
      if((Fail_label[j]==(i+1))&(DELTA[j]==1)){
        sum1[i] += 1;
      }
      if(Fail_label[j]>=(i+1)){
        sum2[i] += f_integral(IF_QUADRATIC, NX, x1, alpha, kns, u1, J, lambda_null, Fail_label[j], K, NU, OT[j], c, i, i+1)*exp(f_mmult(beta1,z1)(0,0)+f_mmult(gamma1,xi1)(0,0));
      }
            
    }
    lambda_out[i] = rgamma(1,alpha1+sum1[i], 1/(alpha2+sum2[i]))[0];
  }
  return List::create(Rcpp::Named("lambda") = lambda_out,
                      Rcpp::Named("sum1") = sum1,
                      Rcpp::Named("sum2") = sum2);  
}
#endif
