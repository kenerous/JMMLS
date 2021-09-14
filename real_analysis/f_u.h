#ifndef _F_U_H
#define _F_U_H

#include "f_inv.h"
#include "f_mmult.h" 
#include "f_utemp.h"
#include "f_madd.h"
#include "f_mvnorm.h"
#include "f_integral.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List f_u(int IF_QUADRATIC, NumericVector u_label, NumericVector psi_eps, NumericVector y, NumericVector b, NumericVector latent_label, int NX, NumericVector x, NumericVector alpha, NumericVector OT_label, NumericVector B, NumericVector kns, NumericMatrix G, NumericVector xi, NumericVector u, NumericVector NT, NumericVector Queue, NumericVector DELTA, NumericVector z, NumericVector OT, NumericVector Fail_label, NumericVector J, NumericVector at_u, NumericVector gamma, NumericVector beta, NumericVector lambda, int N, int NU, int NXI, int NZ, double c, int K, double c_u, int SNT){
  NumericVector u_out = u;
  NumericVector ratio_u(N);
  NumericVector mean_u(NU);
  NumericMatrix u_star(NU,1);
  NumericMatrix u_star2(1,NU);
  NumericMatrix u_old(NU,1);
  NumericMatrix u_old2(1,NU);
  NumericMatrix beta1(1,NZ);
  NumericMatrix z1(NZ,1);
  NumericVector u_star1(NU);
  NumericMatrix gamma1(1,NXI);
  NumericMatrix xi1(NXI,1);
  NumericMatrix multi(NU,NU);
  NumericVector x1(NX+1);
  int dim_eta = latent_label[0];
  
  double eta_temp = 0;
  for(int k=0; k<dim_eta; ++k){
	eta_temp += pow(b[k],2)/psi_eps[k];
  }
  
  for(int i=0; i<N; ++i){
  	int u_scale = u_label[i];
    NumericMatrix sigma_u = f_inv(G);
    int t_temp = NT[i];
    NumericMatrix ft1(NU,1);
    NumericMatrix ft2(1,NU);
    for(int j=0; j<NZ; ++j){
      beta1(0,j) = beta[j];
      z1(j,0) = z[j*N+i];
    }
    for(int j=0; j<NXI; ++j){
      gamma1(0,j) = gamma[j];
      xi1(j,0) = xi[j*N+i];
    }
    for(int j=0; j<(NX+1); ++j){
    	x1[j] = x[j*N+i];
	}
    
    for(int j=0; j<t_temp; ++j){
    	
    	for(int k=0; k<NU; ++k){
        	ft1(k,0) = B[k*SNT+Queue[i]+j];
        	ft2(0,k) = B[k*SNT+Queue[i]+j];
      	}

      sigma_u = f_madd(sigma_u,f_mmult(ft1,ft2)*eta_temp);
    }
    
    //multi = f_utemp(x1, alpha, NX, kns, J, lambda, Fail_label[i], K, NU, OT[i], c, 0, Fail_label[i])*exp(f_mmult(gamma1,xi1)(0,0)+f_mmult(beta1,z1)(0,0));
       
    //sigma_u = f_inv(f_madd(sigma_u,multi));
    sigma_u = f_inv(sigma_u);
    
    for(int j=0; j<NU; ++j){
      mean_u[j] = u[j*N+i];
      u_old(j,0) = u[j*N+i];
      u_old2(0,j) = u[j*N+i];
    }
    
//    NumericMatrix one(NU,NU);
//    for(int j=0; j<NU; ++j){
//    	one(j,j) = 1;
//	}
//    
//    u_star = f_mvnorm(1,mean_u,one*c_u);
    
    u_star = f_mvnorm(1,mean_u,sigma_u*c_u);
    if(u_scale!=NU){
    	for(int j=u_scale; j<NU; ++j){
    		u_star(j,0) = 0;
		}
	}
    
    for(int j=0; j<NU; ++j){
      u_star1[j] = u_star(j,0);
      u_star2(0,j) = u_star(j,0);
    }
    
    NumericMatrix second_level(NU,1);
    for(int j=0; j<NU; ++j){
    	for(int k=0; k<(NX+1); ++k){
    		second_level(j,0) += alpha[j*(NX+1)+k]*x1[k];	
		}       	
  	}
    
    double ratio_new = f_mmult(f_mmult(u_star2,f_inv(G)),u_star)(0,0)/(-2);
    double ratio_old = f_mmult(f_mmult(u_old2,f_inv(G)),u_old)(0,0)/(-2);
    for(int j=0; j<t_temp; ++j){
      for(int k=0; k<NU; ++k){
        ft2(0,k) = B[k*SNT+Queue[i]+j];
      }
      for(int k=0; k<dim_eta; ++k){
      	ratio_new += pow(y[k*SNT+Queue[i]+j]-f_mmult(ft2,f_madd(u_star,second_level))(0,0)*b[k],2)/(-2*psi_eps[k]);
        ratio_old += pow(y[k*SNT+Queue[i]+j]-f_mmult(ft2,f_madd(u_old,second_level))(0,0)*b[k],2)/(-2*psi_eps[k]);
	  }
    }
    
    
    for(int j=0; j<NU; ++j){
      ft2(0,j) = B[j*SNT+Queue[i]+OT_label[i]];
    }// failure time if DELTA=1
    
    
    ratio_new += f_mmult(ft2,u_star)(0,0)*c*DELTA[i]-exp(f_mmult(gamma1,xi1)(0,0)+f_mmult(beta1,z1)(0,0))*f_integral(IF_QUADRATIC, NX, x1, alpha, kns, u_star1, J, lambda, Fail_label[i], K, NU, OT[i], c, 0, Fail_label[i]);
    ratio_old += f_mmult(ft2,u_old)(0,0)*c*DELTA[i]-exp(f_mmult(gamma1,xi1)(0,0)+f_mmult(beta1,z1)(0,0))*f_integral(IF_QUADRATIC, NX, x1, alpha, kns, mean_u, J, lambda, Fail_label[i], K, NU, OT[i], c, 0, Fail_label[i]);
    
    ratio_u[i] = exp(ratio_new-ratio_old);
    
    double randnum = runif(1)[0];
    if(randnum<ratio_u[i]){
      for(int j=0; j<NU; ++j){
        u_out[j*N+i] = u_star(j,0);
      }
      at_u[i] += 1;
    }
      
  }

  return List::create(Rcpp::Named("u") = u_out,
                      Rcpp::Named("at_u") = at_u);
}
#endif
