#ifndef _F_PHI_H
#define _F_PHI_H


#include "f_inv.h"
#include "f_mmult.h"
#include "f_madd.h"
#include "f_wishart.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix f_phi(NumericMatrix phi0, NumericMatrix phi, NumericVector xi, int N, int NXI, double rho_phi0) {
  NumericMatrix phi_out = phi;
  NumericMatrix temp_xi1(NXI,N);
  NumericMatrix temp_xi2(N,NXI);
  
  for(int i=0; i<N; ++i){
    for(int j=0; j<NXI; ++j){
      temp_xi1(j,i)=xi[j*N+i];
      temp_xi2(i,j)=xi[j*N+i];
    }
  }

  phi_out = f_inv(f_wishart(N+rho_phi0,f_inv(f_madd(phi0,f_mmult(temp_xi1,temp_xi2)))));

  
  return phi_out;  
}
#endif

