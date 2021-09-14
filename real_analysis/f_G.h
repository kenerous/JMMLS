#ifndef _F_G_H
#define _F_G_H


#include "f_inv.h"
#include "f_mmult.h"
#include "f_madd.h"
#include "f_wishart.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix f_G(NumericVector u_num, NumericVector u, int N, int NU, int rho0, NumericMatrix G0) {
  NumericMatrix G_out(NU,NU);
  NumericMatrix temp_u1(NU,N);
  NumericMatrix temp_u2(N,NU);
  
  
  for(int i=0; i<N; ++i){
    for(int j=0; j<NU; ++j){
      temp_u1(j,i)=u[j*N+i];
      temp_u2(i,j)=u[j*N+i];
    }
  }
  NumericMatrix sigma_u = f_mmult(temp_u1,temp_u2);
  if(NU>2){
    NumericMatrix sigma_u1(2,2);
    NumericMatrix G01(2,2);
    for(int i=0; i<2; ++i){
      for(int j=0; j<2; ++j){
        sigma_u1(i,j) = sigma_u(i,j);
        G01(i,j) = G0(i,j);
      }
    }
    
    
    NumericMatrix G_out1 = f_inv(f_wishart(N+rho0,f_inv(f_madd(G01,sigma_u1))));
    for(int i=0; i<2; ++i){
      for(int j=0; j<2; ++j){
        G_out(i,j) = G_out1(i,j);
      }
    }
    for(int i=2; i<NU; ++i){
      G_out(i,i) = 1/(rgamma(1,9+0.5*u_num[i],1/(4+0.5*sigma_u(i,i)))[0]);
    }
  }else{
    // NumericMatrix sigma_u1(NU,NU);
    // for(int i=0; i<NU; ++i){
    //   for(int j=0; j<NU; ++j){
    //     sigma_u1(i,j) = sigma_u(i,j);
    //   }
    // }
    G_out = f_inv(f_wishart(N+rho0,f_inv(f_madd(G0,sigma_u))));
  }
  
  
  
  return G_out;  
}

#endif
