#ifndef _F_ALPHA_H
#define _F_ALPHA_H

#include "f_inv.h"
#include "f_mmult.h" 
#include "f_utemp.h"
#include "f_madd.h"
#include "f_mvnorm.h"
#include "f_integral.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List f_alpha(int IF_QUADRATIC, NumericVector psi_eps, NumericVector y, NumericVector latent_label, NumericVector b, int at_alpha, NumericVector alpha0, NumericMatrix Sigma_alpha0, int NX, NumericVector x, NumericVector alpha, NumericVector OT_label, NumericVector B, NumericVector kns, NumericVector xi, NumericVector u, NumericVector NT, NumericVector Queue, NumericVector DELTA, NumericVector z, NumericVector OT, NumericVector Fail_label, NumericVector J, NumericVector gamma, NumericVector beta, NumericVector lambda, int N, int NU, int NXI, int NZ, double c, int K, double c_alpha, int SNT){
  NumericVector alpha_out = alpha;
  double ratio_alpha;
  NumericVector mean_alpha(NU*(NX+1));
  NumericMatrix alpha_star(NU*(NX+1),1);
  NumericMatrix alpha_star2(1,NU*(NX+1));
  NumericVector alpha_star1(NU*(NX+1));
  NumericMatrix alpha_old(NU*(NX+1),1);
  NumericMatrix alpha_old2(1,NU*(NX+1));
  NumericMatrix alpha01(1,NU*(NX+1));
  NumericMatrix alpha02(NU*(NX+1),1);
  
  NumericMatrix beta1(1,NZ);
  NumericMatrix z1(NZ,1);
  NumericMatrix u1(NU,1);
  NumericVector u2(NU);
  NumericMatrix gamma1(1,NXI);
  NumericMatrix xi1(NXI,1);
  NumericVector x1(NX+1);
  int dim_eta = latent_label[0];
  
  double eta_temp = 0;
  for(int k=0; k<dim_eta; ++k){
	eta_temp += pow(b[k],2)/psi_eps[k];
  }
  for(int k=0; k<(NU*(NX+1)); ++k){
  	alpha01(0,k) = alpha0[k];
  	alpha02(k,0) = alpha0[k];
  }
  
  NumericMatrix sigma_alpha = f_inv(Sigma_alpha0);
  
  for(int i=0; i<N; ++i){
  	int t_temp = NT[i];
  	NumericMatrix ft1(NU*(NX+1),1);
    NumericMatrix ft2(1,NU*(NX+1));
  	for(int j=0; j<t_temp; ++j){
  		
  		for(int k=0; k<NU; ++k){	
  			for(int l=0; l<(NX+1); ++l){
				ft1(k*(NX+1)+l,0) = x[l*N+i]*B[k*SNT+Queue[i]+j];
				ft2(0,k*(NX+1)+l) = x[l*N+i]*B[k*SNT+Queue[i]+j];
		  	}
	  	}

	  	sigma_alpha = f_madd(sigma_alpha,f_mmult(ft1,ft2)*eta_temp);
  	}
  }
  
  sigma_alpha = f_inv(sigma_alpha);
  
  for(int j=0; j<(NU*(NX+1)); ++j){
    mean_alpha[j] = alpha[j];
    alpha_old(j,0) = alpha[j];
    alpha_old2(0,j) = alpha[j];
  }
  
  alpha_star = f_mvnorm(1,mean_alpha,sigma_alpha*c_alpha);
  for(int j=0; j<(NU*(NX+1)); ++j){
    alpha_star1[j] = alpha_star(j,0);
    alpha_star2(0,j) = alpha_star(j,0);
  }
  
  double ratio_new = f_mmult(f_mmult(f_madd(alpha_star2,alpha01*(-1)),f_inv(Sigma_alpha0)),f_madd(alpha_star,alpha02*(-1)))(0,0)/(-2);
  double ratio_old = f_mmult(f_mmult(f_madd(alpha_old2,alpha01*(-1)),f_inv(Sigma_alpha0)),f_madd(alpha_old,alpha02*(-1)))(0,0)/(-2);
  for(int i = 0; i<N ; ++i){
  	int t_temp = NT[i];
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
	NumericMatrix second_level_star(NU,1);
	NumericMatrix second_level_old(NU,1);
    for(int j=0; j<NU; ++j){
    	u1(j,0) = u[j*N+i];
    	u2[j] = u[j*N+i];
    	for(int k=0; k<(NX+1); ++k){
    		second_level_star(j,0) += alpha_star1[j*(NX+1)+k]*x1[k];
			second_level_old(j,0) += mean_alpha[j*(NX+1)+k]*x1[k];	
		}       	
  	}
  	
  	for(int j=0; j<t_temp; ++j){
      for(int k=0; k<NU; ++k){
        ft2(0,k) = B[k*SNT+Queue[i]+j];
      }
      for(int k=0; k<dim_eta; ++k){
      	ratio_new += pow(y[k*SNT+Queue[i]+j]-f_mmult(ft2,f_madd(u1,second_level_star))(0,0)*b[k],2)/(-2*psi_eps[k]);
        ratio_old += pow(y[k*SNT+Queue[i]+j]-f_mmult(ft2,f_madd(u1,second_level_old))(0,0)*b[k],2)/(-2*psi_eps[k]);
	  }
    }
    
    for(int j=0; j<NU; ++j){
      ft2(0,j) = B[j*SNT+Queue[i]+OT_label[i]];
    }// failure time if DELTA=1
    
    
    ratio_new += f_mmult(ft2,second_level_star)(0,0)*c*DELTA[i]-exp(f_mmult(gamma1,xi1)(0,0)+f_mmult(beta1,z1)(0,0))*f_integral(IF_QUADRATIC, NX, x1, alpha_star1, kns, u2, J, lambda, Fail_label[i], K, NU, OT[i], c, 0, Fail_label[i]);
    ratio_old += f_mmult(ft2,second_level_old)(0,0)*c*DELTA[i]-exp(f_mmult(gamma1,xi1)(0,0)+f_mmult(beta1,z1)(0,0))*f_integral(IF_QUADRATIC, NX, x1, mean_alpha, kns, u2, J, lambda, Fail_label[i], K, NU, OT[i], c, 0, Fail_label[i]);
    
    
  }
  ratio_alpha = exp(ratio_new-ratio_old);
  double randnum = runif(1)[0];
  if(randnum<ratio_alpha){
    for(int j=0; j<(NU*(NX+1)); ++j){
    	alpha_out[j] = alpha_star(j,0);
      }
    at_alpha += 1;
  }
 
  return List::create(Rcpp::Named("alpha") = alpha_out,
                      Rcpp::Named("at_alpha") = at_alpha);
}
#endif
