#ifndef _F_GAMMA_H
#define _F_GAMMA_H


#include "f_inv.h"
#include "f_mmult.h" 
#include "f_madd.h"
#include "f_mvnorm.h"
#include "f_integral.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List f_gamma(NumericVector kns, NumericVector gamma0, NumericMatrix Sigma_gamma0, NumericVector xi, NumericVector u, NumericVector DELTA, NumericVector z, NumericVector OT, NumericVector Fail_label, NumericVector J, int at_gamma, NumericVector gamma, NumericVector beta, NumericVector Time, NumericVector lambda, int N, int NU, int NXI, int NZ, double c, int K, double c_gamma){
  NumericVector gamma_out = gamma;
  NumericMatrix sigma_gamma = f_inv(Sigma_gamma0);
  NumericMatrix xi1(NXI,1);
  NumericMatrix xi2(1,NXI);
  NumericVector u1(NU);
  NumericMatrix beta1(1,NZ);
  NumericMatrix z1(NZ,1);
  NumericMatrix gamma_star1(1,NXI);
  NumericMatrix gamma_star2(NXI,1);
  NumericVector mean_gamma = gamma;
  NumericMatrix gamma01(1,NXI);
  NumericMatrix gamma02(NXI,1);
  NumericMatrix gamma_old1(1,NXI);
  NumericMatrix gamma_old2(NXI,1);
  NumericVector integral(N);
  
  for(int i=0; i<NZ; ++i){
    beta1(0,i) = beta[i];
  }
  
  
  
  for(int i=0; i<N; ++i){
    for(int j=0; j<NXI; ++j){
      xi1(j,0) = xi[j*N+i];
      xi2(0,j) = xi[j*N+i];
    }
    for(int j=0; j<NU; ++j){
      u1[j] = u[j*N+i];
    }
    for(int j=0; j<NZ; ++j){
      z1(j,0) = z[j*N+i];
    }
    integral[i] = exp(f_mmult(beta1,z1)(0,0))*f_integral(kns, u1, J, lambda, Fail_label[i], K, NU, OT[i], c, 0, Fail_label[i]);
    sigma_gamma = f_madd(sigma_gamma,f_mmult(xi1,xi2)*integral[i]);
  }
  
  sigma_gamma = f_inv(sigma_gamma);
  
  gamma_star2 = f_mvnorm(1,mean_gamma,sigma_gamma*c_gamma);
  
  for(int i=0; i<NXI; ++i){
    gamma_old1(0,i) = gamma[i];
    gamma_old2(i,0) = gamma[i];
    gamma01(0,i) = gamma0[i];
    gamma02(i,0) = gamma0[i];
    gamma_star1(0,i) = gamma_star2(i,0);
  }
  
  double ratio_gamma = -0.5*f_mmult(f_mmult(f_madd(gamma_star1,gamma01*(-1)),f_inv(Sigma_gamma0)),f_madd(gamma_star2,gamma02*(-1)))(0,0)+0.5*f_mmult(f_mmult(f_madd(gamma_old1,gamma01*(-1)),f_inv(Sigma_gamma0)),f_madd(gamma_old2,gamma02*(-1)))(0,0);
  
  for(int i=0; i<N; ++i){
    for(int j=0; j<NXI; ++j){
      xi1(j,0) = xi[j*N+i];
    }
    ratio_gamma += integral[i]*(exp(f_mmult(gamma_old1,xi1)(0,0))-exp(f_mmult(gamma_star1,xi1)(0,0)))+f_mmult(f_madd(gamma_star1,gamma_old1*(-1)),xi1)(0,0)*DELTA[i];
  }
  
  ratio_gamma = exp(ratio_gamma);
  
  double randnum = runif(1)[0];
  if(randnum<ratio_gamma){
    for(int i=0; i<NXI; ++i){
      gamma_out[i] = gamma_star2(i,0);
    }
    at_gamma += 1;
  }
  
  
  
  return List::create(Rcpp::Named("gamma") = gamma_out,
                      Rcpp::Named("at_gamma") = at_gamma,
                      Rcpp::Named("gamma_star1") = gamma_star1);
}
#endif
