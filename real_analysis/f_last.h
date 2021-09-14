#ifndef _F_LAST_H
#define _F_LAST_H

#include "f_inv.h"
#include "f_mmult.h" 
#include "f_madd.h"
#include "f_mvnorm.h"
#include "f_integral.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List f_last(int indi1_start, int indi1_end, int indi2_start, int indi2_end, NumericVector Fail_label_test, NumericVector OT_label_test, int IF_QUADRATIC, NumericMatrix Y_last, NumericMatrix eta_last, NumericMatrix Fail_label_last, double c, NumericVector kns, NumericMatrix u_label_last, NumericMatrix xi_test, NumericVector alpha, NumericMatrix x_test, NumericMatrix y_test, NumericVector b, NumericVector latent_label, NumericVector psi_eps, NumericVector time_test, NumericMatrix B_test, NumericMatrix u_last, NumericVector NT_test, NumericVector Queue_test, NumericVector DELTA_test, NumericMatrix z_test, NumericVector OT_test, NumericVector J, NumericMatrix at_u_last, NumericVector gamma, NumericVector beta, NumericVector lambda, double c_u_last, NumericMatrix G, int K, int S){
  NumericMatrix u_last_out = u_last;
  int indi1_length = indi1_end - indi1_start;
  int indi2_length = indi2_end - indi2_start;
      
  const int N_test = DELTA_test.size();
  const int NU = G.nrow();
  const int NX = x_test.nrow();
  const int NZ = z_test.nrow();
  const int NXI = xi_test.nrow();
  
  
  
  NumericMatrix ratio_u(indi1_end-indi1_start,N_test);
  
  NumericMatrix gamma1(1,NXI);
  NumericMatrix beta1(1,NZ);
  gamma1(0,_) = gamma;
  beta1(0,_) = beta;
  
  NumericMatrix sigma_u0 = f_inv(G);
  

  
  int num_eta = latent_label[0];
  
  double eta_temp = 0;
  for(int k=0; k<num_eta; ++k){
    eta_temp += pow(b[k],2)/psi_eps[k];
  }
  NumericMatrix alpha1(NU,NX);
  for(int k=0; k<NU; ++k){
    alpha1(k,_) = alpha[Range(k*NX,(k+1)*NX-1)];
  }
  NumericMatrix second_level = f_mmult(alpha1,x_test);
  
  NumericVector gammaxi = f_mmult(gamma1,xi_test)(0,_);
  NumericVector betaz = f_mmult(beta1,z_test)(0,_);


        NumericMatrix ft1(NU,1);
        NumericMatrix ft2(1,NU);
        NumericMatrix second_level1(NU,1);
        NumericMatrix xi1(NXI,1);
        NumericMatrix z1(NZ,1);
        
        for(int indi1=indi1_start; indi1<indi1_end; ++indi1){
          for(int i=0; i<N_test; ++i){
            
            if(NT_test[i]>(indi1+1)){ // at least indi1+2 observations (So we can compare the prediction with the true value of the indi1+2-th observation)
              
              double t_temp1 = time_test[Queue_test[i]+indi1];
              
              //if(!((OT_test[i]<=t_temp1)&(DELTA_test[i]==1))){
              
              NumericMatrix sigma_u = sigma_u0;
              NumericVector u_mean_temp = u_last(_,i);
              NumericVector u_mean = u_mean_temp[Range((indi1-indi1_start)*NU,(indi1-indi1_start+1)*NU-1)];
              // NumericVector u_mean = u_last(_,i);
              NumericMatrix u_old1(NU,1);
              NumericMatrix u_old2(1,NU);
              u_old1(_,0) = u_mean;
              u_old2(0,_) = u_mean;
              
              int label_temp = 0;
              double t_temp = time_test[Queue_test[i]];
              
              
              while(t_temp<=t_temp1){
                ft1(_,0) = B_test(_,Queue_test[i]+label_temp);
                ft2(0,_) = B_test(_,Queue_test[i]+label_temp);
                
                
                sigma_u = f_madd(sigma_u,f_mmult(ft1,ft2)*eta_temp);
                
                label_temp += 1;
                t_temp = time_test[Queue_test[i]+label_temp];
                
              }
              sigma_u = f_inv(sigma_u);
              
              NumericMatrix u_star1 = f_mvnorm(1,u_mean,sigma_u*c_u_last);
              
              int u_scale = u_label_last((indi1-indi1_start),i);
              if(u_scale!=NU){
                for(int j=u_scale; j<NU; ++j){
                  u_star1(j,0) = 0;
                }
              }
              NumericMatrix u_star2(1,NU);
              u_star2(0,_) = u_star1(_,0);
              
              // calculate the acceptance
              double ratio_new = f_mmult(u_star2,f_mmult(sigma_u0,u_star1))(0,0)*(-0.5);
              double ratio_old = f_mmult(u_old2,f_mmult(sigma_u0,u_old1))(0,0)*(-0.5);
              
              label_temp = 0;
              t_temp = time_test[Queue_test[i]];
              
              second_level1(_,0) = second_level(_,i);
              while(t_temp<=t_temp1){
                
                ft2(0,_) = B_test(_,Queue_test[i]+label_temp);
                for(int k=0; k<num_eta; ++k){
                  ratio_new += pow(y_test(k,Queue_test[i]+label_temp)-f_mmult(ft2,f_madd(u_star1,second_level1))(0,0)*b[k],2)/(-2*psi_eps[k]);
                  ratio_old += pow(y_test(k,Queue_test[i]+label_temp)-f_mmult(ft2,f_madd(u_old1,second_level1))(0,0)*b[k],2)/(-2*psi_eps[k]);
                }
                
                label_temp += 1;
                t_temp = time_test[Queue_test[i]+label_temp];
                
              }
              
              if((DELTA_test[i]==1)&(OT_label_test[i]<=indi1)){
                
                ft2(0,_) = B_test(_,Queue_test[i]+OT_label_test[i]);
                ratio_new += f_mmult(ft2,u_star1)(0,0)*c - exp(gammaxi[i]+betaz[i])*f_integral(IF_QUADRATIC, NX-1, x_test(_,i), alpha, kns, u_star1(_,0), J, lambda, Fail_label_test[i], K, NU, OT_test[i], c, 0, Fail_label_test[i]);
                ratio_old += f_mmult(ft2,u_old1)(0,0)*c - exp(gammaxi[i]+betaz[i])*f_integral(IF_QUADRATIC, NX-1, x_test(_,i), alpha, kns, u_mean, J, lambda, Fail_label_test[i], K, NU, OT_test[i], c, 0, Fail_label_test[i]);
                
              }else{
                ratio_new += (-1)*exp(gammaxi[i]+betaz[i])*f_integral(IF_QUADRATIC, NX-1, x_test(_,i), alpha, kns, u_star1(_,0), J, lambda, Fail_label_last(indi1-indi1_start,i), K, NU, t_temp1, c, 0, Fail_label_last(indi1-indi1_start,i));
                ratio_old += (-1)*exp(gammaxi[i]+betaz[i])*f_integral(IF_QUADRATIC, NX-1, x_test(_,i), alpha, kns, u_mean, J, lambda, Fail_label_last(indi1-indi1_start,i), K, NU, t_temp1, c, 0, Fail_label_last(indi1-indi1_start,i));
                
              }
              
              
              ratio_u(indi1-indi1_start,i) = exp(ratio_new-ratio_old);
              
              double randnum = runif(1)[0];
              if(randnum<ratio_u(indi1-indi1_start,i)){
                NumericVector u_last_temp(u_last_out(_,i));
                u_last_temp[Range((indi1-indi1_start)*NU,(indi1-indi1_start+1)*NU-1)] = u_star1(_,0);
                u_last_out(_,i) = u_last_temp;
                // u_last_out(_,i) = u_star1(_,0);
                at_u_last(indi1-indi1_start,i) += 1;
              }
              
              NumericVector u_temp_temp(u_last_out(_,i));
              NumericMatrix u_temp(NU,1);
              u_temp(_,0) = u_temp_temp[Range((indi1-indi1_start)*NU,(indi1-indi1_start+1)*NU-1)];
              
              for(int indi2=indi2_start; indi2<indi2_end; ++indi2){
                if(NT_test[i]>(indi1+indi2)){
                    
                    ft2(0,_) = B_test(_,Queue_test[i]+indi1+indi2);
                  double latent = f_mmult(ft2,f_madd(u_temp,second_level1))(0,0)/S;
                  eta_last((indi1-indi1_start)*indi2_length+indi2-indi2_start,i) += latent;
                  
                  for(int j=0; j<num_eta; ++j){
                    Y_last((indi1-indi1_start)*indi2_length*num_eta+(indi2-indi2_start)*num_eta+j,i) += b[j]*latent;
                  } 
                  
                  
                }
              }
              
              
              
              
              
              //}
              
            }
            
            
          }
        }
        

  
  NumericVector u_last_out_vec(NU*indi1_length*N_test);
  NumericVector Y_last_vec(indi1_length*indi2_length*num_eta*N_test);
  NumericVector at_u_last_vec(indi1_length*N_test);
  NumericVector eta_last_vec(indi1_length*indi2_length*N_test);
  for(int k=0; k<(indi1_length*NU); ++k){
    u_last_out_vec[Range(k*N_test,(k+1)*N_test-1)] = u_last_out(k,_);
  }
  for(int k=0; k<(indi1_length*indi2_length*num_eta); ++k){
    Y_last_vec[Range(k*N_test,(k+1)*N_test-1)] = Y_last(k,_);
  }
  for(int k=0; k<indi1_length; ++k){
    at_u_last_vec[Range(k*N_test,(k+1)*N_test-1)] = at_u_last(k,_);
  }
  for(int k=0; k<(indi1_length*indi2_length); ++k){
    eta_last_vec[Range(k*N_test,(k+1)*N_test-1)] = eta_last(k,_);
  }
  return List::create(Rcpp::Named("u_last_vec") = u_last_out_vec,
                      Rcpp::Named("at_u_last_vec") = at_u_last_vec,
                      Rcpp::Named("eta_last_vec") = eta_last_vec,
                      Rcpp::Named("Y_last_vec") = Y_last_vec
                        
  );
}
#endif
