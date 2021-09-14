#include <RcppArmadillo.h>

#include "f_det.h"
#include "f_mmult.h"
#include "f_utemp.h"
#include "f_eta.h"
#include "f_pre.h"
#include "f_last.h"
#include "f_u.h"
#include "f_G.h"
#include "f_b.h"
#include "f_psieps.h"
#include "f_lambda.h"
#include "f_beta.h"
#include "f_c.h"



using namespace Rcpp;

// [[Rcpp::export]]
List mcmc(List observed_data, List test_data, List para, List mh, List prior, int K, int iter, int rep, int T, int S, int IF_PRE, int IF_LAST, int IF_QUADRATIC, int IF_SEM){
  int rho0 = prior[0];
  NumericMatrix G0 = prior[1];
  double b0 = prior[2];
  double h0 = prior[3];
  double alpha_eps0 = prior[4];
  double beta_eps0 = prior[5];
  double alpha1 = prior[6];
  double alpha2 = prior[7];
  NumericVector beta0 = prior[8];
  NumericMatrix Sigma_beta0 = prior[9];
  double c0 = prior[10];
  double sigma_c0 = prior[11];
  NumericVector alpha0 = prior[12];
  NumericMatrix Sigma_alpha0 = prior[13];

  double c_u = mh[0];
  double c_beta = mh[1];
  double c_c = mh[2];
  double c_u_test = mh[4];
  double c_u_last = mh[5];

  
  NumericVector u = para[0];
  NumericMatrix G = para[1];
  NumericVector b = para[2];
  NumericVector psi_eps = para[3];
  NumericVector lambda = para[4];
  NumericVector beta = para[5];
  NumericVector alpha = para[6];
  double c = para[7];
  NumericVector xi = para[8];
  NumericVector gamma = para[9];

  
  NumericVector x = observed_data[0];
  NumericVector OT_label = observed_data[1];
  NumericVector B = observed_data[2];
  NumericVector kns = observed_data[3];
  NumericVector y = observed_data[4];
  NumericVector b_label = observed_data[5];
  NumericVector J =observed_data[6];
  NumericVector Fail_label = observed_data[7];
  NumericVector z = observed_data[8];
  NumericVector OT = observed_data[9];
  NumericVector DELTA = observed_data[10];
  NumericVector NT = observed_data[11];
  NumericVector Queue = observed_data[12];
  NumericVector latent_label = observed_data[13];
  NumericVector u_label = observed_data[14];
  NumericVector u_num = observed_data[15];
  
  NumericMatrix x_test = test_data[0];
  NumericVector OT_label_test = test_data[1];
  NumericMatrix B_test = test_data[2];
  NumericMatrix Fail_label_tdelta = test_data[3];
  NumericMatrix y_test = test_data[4];
  NumericVector t_test = test_data[5];
  NumericVector delta_test = test_data[6];
  NumericVector Fail_label_test = test_data[7];
  NumericMatrix z_test = test_data[8];
  NumericVector OT_test = test_data[9];
  NumericVector DELTA_test = test_data[10];
  NumericVector NT_test = test_data[11];
  NumericVector Queue_test = test_data[12];
  NumericVector time_test = test_data[13];
  NumericVector u_label_test = test_data[14];
  NumericMatrix xi_test = test_data[15];
  NumericMatrix u_test = test_data[16];
  NumericMatrix u_label_last = test_data[17];
  NumericMatrix u_last = test_data[18];
  NumericMatrix Fail_label_last = test_data[19];
  int indi1_start = test_data[20];
  int indi1_end = test_data[21];
  int indi2_start = test_data[22];
  int indi2_end = test_data[23];
  
  
  const int num_eta = latent_label[0];
  int t_length = t_test.size();
  int delta_length = delta_test.size();
  
  int NY = b.size();
  int NU = G.ncol();
  int SNT = B.size()/NU;
  int N = u.size()/NU;
  int NX = alpha.size()/NU - 1;
  int NZ = beta.size();
  int NXI = gamma.size();
  int J_len = J.size();
  
  NumericVector eta_str(S*SNT);
  NumericVector u_str(S*N*NU);
  NumericVector G_str(S*NU*NU);
  NumericVector b_str(S*NY);
  NumericVector psieps_str(S*NY);
  NumericVector lambda_str(S*(J_len-1));
  NumericVector beta_str(S*NZ);
  NumericVector gamma_str(S*NXI);
  NumericVector c_str(S);
  NumericVector phi_str(S*NXI*NXI);
  NumericVector xi_str(S*N*NXI);
  NumericVector alpha_str(S*NU*(NX+1));
 
  
  NumericVector eta(SNT);
  
  NumericVector at_u(N);
  NumericVector at_xi(N);
  int at_beta = 0;
  // int at_gamma = 0;
  int at_c = 0;
  int at_alpha = 0;
  
  
  NumericMatrix numerator_pair(t_length*delta_length,4);
  NumericMatrix denominator_pair(t_length*delta_length,4);
  NumericVector numerator_pair_vec(t_length*delta_length*4);
  NumericVector denominator_pair_vec(t_length*delta_length*4);
  NumericMatrix AUC(S,t_length*delta_length);
  const int N_test = DELTA_test.size();
  NumericMatrix pi_test(t_length*delta_length,N_test);
  NumericMatrix at_u_test(t_length*delta_length,N_test);
  NumericMatrix all_u_test(t_length*delta_length,N_test);
  NumericVector pi_test_vec(t_length*delta_length*N_test);
  NumericVector at_u_test_vec(t_length*delta_length*N_test);
  NumericVector all_u_test_vec(t_length*delta_length*N_test);
  NumericVector u_test_vec(t_length*delta_length*NU*N_test);
  NumericVector AUC_vec(t_length*delta_length);
  
  
  int indi1_length = indi1_end - indi1_start;
  int indi2_length = indi2_end - indi2_start;
  NumericMatrix Y_last(indi1_length*indi2_length*num_eta,N_test);
  NumericMatrix eta_last(indi1_length*indi2_length,N_test);
  NumericMatrix at_u_last(indi1_length,N_test);
  NumericVector Y_last_vec(indi1_length*indi2_length*num_eta*N_test);
  NumericVector u_last_vec(indi1_length*NU*N_test);
  NumericVector at_u_last_vec(indi1_length*N_test);
  NumericVector eta_last_vec(indi1_length*indi2_length*N_test);
  
  
  for(int t=0; t<T; ++t){
    
    
    List u_all = f_u(IF_QUADRATIC, u_label, psi_eps, y, b, latent_label, NX, x, alpha, OT_label, B, kns, G, xi, u, NT, Queue, DELTA, z, OT, Fail_label, J, at_u, gamma, beta, lambda, N, NU, NXI, NZ, c, K, c_u, SNT);
    u = u_all[0];
    at_u = u_all[1];

    G = f_G(u_num, u, N, NU, rho0, G0);

    eta = f_eta(x, alpha, NX, B, u, NT, Queue, N, NU, SNT);
    
    if(IF_SEM==1){
      b = f_b(u, alpha, B, y, b_label, xi, eta, b, Queue, latent_label, psi_eps, N, SNT, NY, b0, h0);
    }
    
    psi_eps = f_psieps(y, xi, eta, b, b_label, Queue, latent_label, psi_eps, N, SNT, NY, alpha_eps0, beta_eps0);

    List lambda_all = f_lambda(IF_QUADRATIC, NX, x, alpha, kns, J, lambda, xi, z, OT, DELTA, u, Fail_label, beta, gamma, N, NU, K, NZ, NXI, alpha1, alpha2, c);
    lambda = lambda_all[0];

    List beta_all = f_beta(IF_QUADRATIC, NX, x, alpha, kns, beta0, Sigma_beta0, xi, u, DELTA, z, OT, Fail_label, J, at_beta, gamma, beta, lambda, N, NU, NXI, NZ, c, K, c_beta);
    beta = beta_all[0];
    at_beta = beta_all[1];


    List c_all = f_c(IF_QUADRATIC, NX ,x, alpha, Queue, SNT, OT_label, B, kns, c0, sigma_c0, xi, u, DELTA, z, OT, Fail_label, J, at_c, gamma, beta, lambda, N, NU, NXI, NZ, c, K, c_c);
    c = c_all[0];
    at_c = c_all[2];

    
    if(t>(T-S-1)){
      int ind = t-(T-S);
      for(int i=0; i<NU; ++i){

        for(int j=0; j<NU; ++j){
          G_str[ind*NU*NU+i*NU+j] = G(i,j);
        }
        for(int j=0; j<N; ++j){
          u_str[ind*N*NU+i*N+j] = u[i*N+j];
        }
      }
      
      for(int i=0; i<N; ++i){
        for(int j=0; j<NXI; ++j){
          xi_str[ind*N*NXI+j*N+i] = xi[j*N+i]; 
        }
      }
      for(int i=0; i<SNT; ++i){
        eta_str[ind*SNT+i] = eta[i];
      }
      for(int i=0; i<NY; ++i){
        b_str[ind*NY+i] = b[i];
        psieps_str[ind*NY+i] = psi_eps[i];
      }
      c_str[ind] = c;
      for(int i=0; i<(J_len-1); ++i){
        lambda_str[ind*(J_len-1)+i] = lambda[i];
      }
      for(int i=0; i<NZ; ++i){
        beta_str[ind*NZ+i] = beta[i];
      }
      for(int i=0; i<NXI; ++i){
        gamma_str[ind*NXI+i] = gamma[i];
      }
      for(int i=0; i<(NU*(NX+1)); ++i){
        alpha_str[ind*NU*(NX+1)+i] = alpha[i];
      }
      
      // dynamic prediction
      if(IF_PRE==1){

        List pre_all = f_pre(IF_QUADRATIC, all_u_test, c, kns, u_label_test, xi_test, alpha, x_test, y_test, b, latent_label, psi_eps, Fail_label_tdelta, t_test, delta_test, time_test, B_test, u_test, NT_test, Queue_test, DELTA_test, z_test, OT_test, J, at_u_test, gamma, beta, lambda, c_u_test, G, K, S);
        u_test_vec = pre_all[0];
        at_u_test_vec = pre_all[1];
        all_u_test_vec = pre_all[2];
        numerator_pair_vec = pre_all[3];
        denominator_pair_vec = pre_all[4];
        AUC_vec = pre_all[5];
        pi_test_vec = pre_all[6];

        for(int i=0; i<t_length; ++i){
          for(int j=0; j<delta_length; ++j){
            for(int k=0; k<NU; ++k){
              u_test((i*delta_length+j)*NU+k,_) = u_test_vec[Range((i*delta_length+j)*NU*N_test+k*N_test,(i*delta_length+j)*NU*N_test+(k+1)*N_test-1)];

            }
            at_u_test(i*delta_length+j,_) = at_u_test_vec[Range((i*delta_length+j)*N_test,(i*delta_length+j+1)*N_test-1)];
            all_u_test(i*delta_length+j,_) = all_u_test_vec[Range((i*delta_length+j)*N_test,(i*delta_length+j+1)*N_test-1)];

            pi_test(i*delta_length+j,_) = pi_test_vec[Range((i*delta_length+j)*N_test,(i*delta_length+j+1)*N_test-1)];
            numerator_pair(i*delta_length+j,_) = numerator_pair_vec[Range((i*delta_length+j)*4,(i*delta_length+j+1)*4-1)];
            denominator_pair(i*delta_length+j,_) = denominator_pair_vec[Range((i*delta_length+j)*4,(i*delta_length+j+1)*4-1)];

          }
        }
        AUC(ind,_) = AUC_vec;



      }
      // last-visit prediction
      if(IF_LAST==1){
        // given the first indi1_start+1,...,indi1_end longitudinal observations, predict the next indi2_start-1,...,indi2_end-1 observations  
        
        List last_all = f_last(indi1_start, indi1_end, indi2_start, indi2_end, Fail_label_test, OT_label_test, IF_QUADRATIC, Y_last, eta_last, Fail_label_last, c, kns, u_label_last, xi_test, alpha, x_test, y_test, b, latent_label, psi_eps, time_test, B_test, u_last, NT_test, Queue_test, DELTA_test, z_test, OT_test, J, at_u_last, gamma, beta, lambda, c_u_last, G, K, S);
        u_last_vec = last_all[0];
        at_u_last_vec = last_all[1];
        eta_last_vec = last_all[2];
        Y_last_vec = last_all[3];
        
        for(int k=0; k<(indi1_length*NU); ++k){
          u_last(k,_) = u_last_vec[Range(k*N_test,(k+1)*N_test-1)];
        }
        for(int k=0; k<(indi1_length*indi2_length*num_eta); ++k){
          Y_last(k,_) = Y_last_vec[Range(k*N_test,(k+1)*N_test-1)];
        }
        for(int k=0; k<indi1_length; ++k){
          at_u_last(k,_) = at_u_last_vec[Range(k*N_test,(k+1)*N_test-1)];
        }
        for(int k=0; k<(indi1_length*indi2_length); ++k){
          eta_last(k,_) = eta_last_vec[Range(k*N_test,(k+1)*N_test-1)];
        }
            
      }
      

    }
    
    
    
    if((t+1)%100==0){
      Rprintf("\rRunning MCMC, %d/%d replication, %f%% completed...",iter,rep,(t+1)*100.0/T);
    }
    
    
  }
  
  return List::create(Rcpp::Named("u_str") = u_str,
                      Rcpp::Named("G_str") = G_str,
                      Rcpp::Named("b_str") = b_str,
                      Rcpp::Named("AUC") = AUC,
                      Rcpp::Named("pi_test") = pi_test,
                      Rcpp::Named("u_test") = u_test,
                      Rcpp::Named("at_u_test") = at_u_test,
                      Rcpp::Named("Y_last") = Y_last,
                      Rcpp::Named("eta_last") = eta_last,
                      Rcpp::Named("u_last") = u_last,
                      Rcpp::Named("psieps_str") = psieps_str,
                      Rcpp::Named("lambda_str") = lambda_str,
                      Rcpp::Named("beta_str") = beta_str,
                      Rcpp::Named("alpha_str") = alpha_str,
                      Rcpp::Named("c_str") = c_str,
                      Rcpp::Named("at_u") = at_u
                      ,Rcpp::Named("at_c") = at_c
                      ,Rcpp::Named("at_alpha") = at_alpha
                      ,Rcpp::Named("at_beta") = at_beta
                        );
}