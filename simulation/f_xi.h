#ifndef _F_XI_H
#define _F_XI_H

 
#include "f_inv.h"
#include "f_mmult.h" 
#include "f_madd.h"
#include "f_mvnorm.h"
#include "f_integral.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector f_xi(NumericVector kns, NumericVector y, NumericVector xi, NumericVector u, NumericVector b, NumericVector Queue, NumericVector latent_label, NumericVector psi_eps, NumericVector DELTA, NumericVector phi, NumericVector z, NumericVector OT, NumericVector lambda, NumericVector Fail_label, NumericVector J, NumericVector at_xi, NumericVector gamma, NumericVector beta, int N, int NU, int NY, int NXI, int NZ, double c, int K, int SNT, double c_xi){
  NumericVector xi_out = xi;
  NumericVector ratio_xi(N);
  NumericMatrix sigma_xi(NXI,NXI);
  NumericMatrix sigma_phi(NXI,NXI);
  NumericMatrix sigma_eps(NY-latent_label[0],NY-latent_label[0]);
  NumericMatrix B1(NXI,NY-latent_label[0]);
  NumericMatrix B2(NY-latent_label[0],NXI);
  NumericMatrix gamma1(NXI,1);
  NumericMatrix gamma2(1,NXI);
  NumericMatrix sigma_gamma(NXI,NXI);
  NumericMatrix sigma_temp(NXI,NXI);
  NumericMatrix beta1(1,NZ);
  NumericMatrix z1(NZ,1);
  NumericVector mean_xi(NXI);
  NumericMatrix xi_star1(NXI,1);
  NumericMatrix xi_star2(1,NXI);
  NumericMatrix xi_old1(NXI,1);
  NumericMatrix xi_old2(1,NXI);
  NumericMatrix y1(NY-latent_label[0],1);
  NumericMatrix y2(1,NY-latent_label[0]);
  
  int temp = NY-latent_label[0];
  for(int i=0; i<NXI; ++i){
    sigma_phi(i,i) = phi[i];
  }
  for(int i=0; i<temp; ++i){
    sigma_eps(i,i) =psi_eps[i+latent_label[0]];
  }
  int temp1 = latent_label[0];
  for(int i=0; i<NXI; ++i){
    int temp2 = temp1+latent_label[i+1];
    for(int j=temp1; j<temp2; ++j){
      B1(i,j-latent_label[0]) = b[j];
      B2(j-latent_label[0],i) = b[j];
    }
    temp1 = temp2;
  }

  sigma_temp = f_madd(f_mmult(f_mmult(B1,sigma_eps),B2),f_inv(sigma_phi));

  for(int i=0; i<NXI; ++i){
    gamma1(i,0) = gamma[i];
    gamma2(0,i) = gamma[i];
  }
  sigma_gamma = f_mmult(gamma1,gamma2);

  for(int i=0; i<N; ++i){
    NumericVector u_temp(NU);
    for(int j=0; j<NU; ++j){
      u_temp[j] = u[j*N+i];
    }
    for(int j=0; j<NZ; ++j){
      beta1(0,j) = beta[j];
      z1(j,0) = z[j*N+i];
    }
    double integral = f_integral(kns, u_temp, J, lambda, Fail_label[i], K, NU, OT[i], c, 0, Fail_label[i]);
    sigma_xi = f_inv(f_madd(sigma_temp,sigma_gamma*(integral*exp(f_mmult(beta1,z1)(0,0)))));
    for(int j=0; j<NXI; ++j){
      mean_xi[j] = xi[j*N+i];
    }
    xi_star1 = f_mvnorm(1,mean_xi,sigma_xi*c_xi);
    
    for(int j=0; j<NXI; ++j){
      xi_star2(0,j) = xi_star1(j,0);
      xi_old1(j,0) = mean_xi[j];
      xi_old2(0,j) = mean_xi[j];
    }
    for(int j=0; j<temp; ++j){
      y1(j,0) = y[(latent_label[0]+j)*SNT+Queue[i]];
      y2(0,j) = y[(latent_label[0]+j)*SNT+Queue[i]];
    }
    
    double ratio_new = f_mmult(f_mmult(xi_star2,f_inv(sigma_phi)),xi_star1)(0,0)*(-0.5)+f_mmult(f_mmult(f_madd(y2,f_mmult(xi_star2,B1)*(-1)),f_inv(sigma_eps)),f_madd(y1,f_mmult(B2,xi_star1)*(-1)))(0,0)*(-0.5)+f_mmult(gamma2,xi_star1)(0,0)*DELTA[i]-integral*exp(f_mmult(beta1,z1)(0,0)+f_mmult(gamma2,xi_star1)(0,0));
    double ratio_old = f_mmult(f_mmult(xi_old2,f_inv(sigma_phi)),xi_old1)(0,0)*(-0.5)+f_mmult(f_mmult(f_madd(y2,f_mmult(xi_old2,B1)*(-1)),f_inv(sigma_eps)),f_madd(y1,f_mmult(B2,xi_old1)*(-1)))(0,0)*(-0.5)+f_mmult(gamma2,xi_old1)(0,0)*DELTA[i]-integral*exp(f_mmult(beta1,z1)(0,0)+f_mmult(gamma2,xi_old1)(0,0));
     
    ratio_xi[i] = exp(ratio_new-ratio_old);
    
  
    double randnum = runif(1)[0];
    if(randnum<ratio_xi[i]){
      for(int j=0; j<NXI; ++j){
        xi_out[j*N+i] = xi_star1(j,0);
      }
      at_xi[i] += 1;
    }

  }
  
  
  return xi_out;
}


#endif
