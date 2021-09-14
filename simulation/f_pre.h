#ifndef _F_PRE_H
#define _F_PRE_H

#include "f_inv.h"
#include "f_mmult.h" 
#include "f_madd.h"
#include "f_mvnorm.h"
#include "f_integral.h"

#include <Rcpp.h> 
using namespace Rcpp;

// [[Rcpp::export]]
List f_pre(int IF_QUADRATIC, NumericMatrix all_u_test, double c, NumericVector kns, NumericVector u_label_test, NumericMatrix xi_test, NumericVector alpha, NumericMatrix x_test, NumericMatrix y_test, NumericVector b, NumericVector latent_label, NumericVector psi_eps, NumericMatrix Fail_label_tdelta, NumericVector t_test, NumericVector delta_test, NumericVector time_test, NumericMatrix B_test, NumericMatrix u_test, NumericVector NT_test, NumericVector Queue_test, NumericVector DELTA_test, NumericMatrix z_test, NumericVector OT_test, NumericVector J, NumericMatrix at_u_test, NumericVector gamma, NumericVector beta, NumericVector lambda, double c_u_test, NumericMatrix G, int K, int S){
  NumericMatrix u_test_out = u_test;
  
    
    
  const int N_test = DELTA_test.size();
  const int NU = G.nrow();
  const int NX = x_test.nrow();
  const int NZ = z_test.nrow();
  const int NXI = xi_test.nrow();
  
  const int t_length = t_test.size();
  const int delta_length = delta_test.size();
  NumericVector AUC(t_length*delta_length);
  
  NumericMatrix pi_test(t_length*delta_length,N_test);
  NumericMatrix denominator_pair(t_length*delta_length,4);
  NumericMatrix numerator_pair(t_length*delta_length,4);
  
  
  NumericMatrix ratio_u(t_length*delta_length,N_test);
  
  NumericMatrix gamma1(1,NXI);
  NumericMatrix beta1(1,NZ);
  gamma1(0,_) = gamma;
  beta1(0,_) = beta;
  
  NumericMatrix sigma_u0 = f_inv(G);
  

  
  int dim_eta = latent_label[0];
  
  double eta_temp = 0;
  for(int k=0; k<dim_eta; ++k){
    eta_temp += pow(b[k],2)/psi_eps[k];
  }
  NumericMatrix alpha1(NU,NX);
  for(int k=0; k<NU; ++k){
    alpha1(k,_) = alpha[Range(k*NX,(k+1)*NX-1)];
  }
  NumericMatrix second_level = f_mmult(alpha1,x_test);
  
  NumericVector gammaxi = f_mmult(gamma1,xi_test)(0,_);
  NumericVector betaz = f_mmult(beta1,z_test)(0,_);
    for(int loop1=0; loop1<t_length; ++loop1){
      for(int loop2=0; loop2<delta_length; ++loop2){

        double t_temp1 = t_test[loop1];
        double delta_temp1 = delta_test[loop2];

        int location = (loop1*delta_length + loop2 ) *NU;

        NumericMatrix ft1(NU,1);
        NumericMatrix ft2(1,NU);
        NumericMatrix second_level1(NU,1);
        NumericMatrix xi1(NXI,1);
        NumericMatrix z1(NZ,1);
        for(int i=0; i<N_test; ++i){

          if(!((OT_test[i]<=t_temp1)&(DELTA_test[i]==1))){
            all_u_test(loop1*delta_length+loop2,i) += 1;
            
            NumericMatrix sigma_u = sigma_u0;
            NumericVector u_mean(NU);
            NumericVector u_mean_temp = u_test(_,i);
            u_mean = u_mean_temp[Range(location,location+NU-1)];
            NumericMatrix u_old1(NU,1);
            NumericMatrix u_old2(1,NU);
            u_old1(_,0) = u_mean;
            u_old2(0,_) = u_mean;

            int label_temp = 0;
            double t_temp = time_test[Queue_test[i]];


            while((t_temp<=t_temp1)&(label_temp<NT_test[i])){
              ft1(_,0) = B_test(_,Queue_test[i]+label_temp);
              ft2(0,_) = B_test(_,Queue_test[i]+label_temp);


              sigma_u = f_madd(sigma_u,f_mmult(ft1,ft2)*eta_temp);

              label_temp += 1;
              t_temp = time_test[Queue_test[i]+label_temp];

            }
            sigma_u = f_inv(sigma_u);

            NumericMatrix u_star1 = f_mvnorm(1,u_mean,sigma_u*c_u_test);

            int u_scale = u_label_test[loop1];
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
            while((t_temp<=t_temp1)&(label_temp<NT_test[i])){
              ft2(0,_) = B_test(_,Queue_test[i]+label_temp);
              for(int k=0; k<dim_eta; ++k){
                ratio_new += pow(y_test(k,Queue_test[i]+label_temp)-f_mmult(ft2,f_madd(u_star1,second_level1))(0,0)*b[k],2)/(-2*psi_eps[k]);
                ratio_old += pow(y_test(k,Queue_test[i]+label_temp)-f_mmult(ft2,f_madd(u_old1,second_level1))(0,0)*b[k],2)/(-2*psi_eps[k]);
              }

              label_temp += 1;
              t_temp = time_test[Queue_test[i]+label_temp];

            }

            ratio_new += (-1)*exp(gammaxi[i]+betaz[i])*f_integral(IF_QUADRATIC, NX-1, x_test(_,i), alpha, kns, u_star1(_,0), J, lambda, Fail_label_tdelta(0,loop1*delta_length+loop2), K, NU, t_temp1, c, 0, Fail_label_tdelta(0,loop1*delta_length+loop2));
            ratio_old += (-1)*exp(gammaxi[i]+betaz[i])*f_integral(IF_QUADRATIC, NX-1, x_test(_,i), alpha, kns, u_mean, J, lambda, Fail_label_tdelta(0,loop1*delta_length+loop2), K, NU, t_temp1, c, 0, Fail_label_tdelta(0,loop1*delta_length+loop2));

            ratio_u(loop1*delta_length+loop2,i) = exp(ratio_new-ratio_old);

            double randnum = runif(1)[0];
            if(randnum<ratio_u(loop1*delta_length+loop2,i)){
              NumericVector u_test_out_temp = u_test_out(_,i);
              u_test_out_temp[Range(location,location+NU-1)] = u_star1(_,0);
              u_test_out(_,i) = u_test_out_temp;
              at_u_test(loop1*delta_length+loop2,i) += 1;
            }
            
            NumericVector u_test_out_temp = u_test_out(_,i);

            pi_test(loop1*delta_length+loop2,i) = exp( exp(gammaxi[i]+betaz[i])*( (-1) *f_integral(IF_QUADRATIC, NX-1, x_test(_,i), alpha, kns, u_test_out_temp[Range(location,location+NU-1)], J, lambda, Fail_label_tdelta(1,loop1*delta_length+loop2), K, NU, t_temp1+delta_temp1, c, 0, Fail_label_tdelta(1,loop1*delta_length+loop2)) + f_integral(IF_QUADRATIC, NX-1, x_test(_,i), alpha, kns, u_test_out_temp[Range(location,location+NU-1)], J, lambda, Fail_label_tdelta(0,loop1*delta_length+loop2), K, NU, t_temp1, c, 0, Fail_label_tdelta(0,loop1*delta_length+loop2))   ) );



          }



        }

              // calculate AUC value
              for(int l1=0; l1<N_test; ++l1){
                for(int l2=0; l2<N_test; ++l2){

                  if(l1!=l2){

                    // four senarios for the pairs (l1,l2)
                    if((OT_test[l1]>t_temp1)&(OT_test[l1]<=(t_temp1+delta_temp1))&(DELTA_test[l1]==1)&(OT_test[l2]>(t_temp1+delta_temp1))){
                      denominator_pair(loop1*delta_length+loop2,0) += 1;
                      if(pi_test(loop1*delta_length+loop2,l1)<pi_test(loop1*delta_length+loop2,l2)){
                        numerator_pair(loop1*delta_length+loop2,0) += 1;
                      }


                    }else if((OT_test[l1]>t_temp1)&(OT_test[l1]<=(t_temp1+delta_temp1))&(DELTA_test[l1]==0)&(OT_test[l2]>(t_temp1+delta_temp1))){

                      denominator_pair(loop1*delta_length+loop2,1) += 1-pi_test(loop1*delta_length+loop2,l1);
                      if(pi_test(loop1*delta_length+loop2,l1)<pi_test(loop1*delta_length+loop2,l2)){
                        numerator_pair(loop1*delta_length+loop2,1) += 1-pi_test(loop1*delta_length+loop2,l1);
                      }

                    }else if((OT_test[l1]>t_temp1)&(OT_test[l1]<=(t_temp1+delta_temp1))&(DELTA_test[l1]==1)&(OT_test[l2]<=(t_temp1+delta_temp1))&(DELTA_test[l2]==0)&(OT_test[l1]<OT_test[l2])){

                      denominator_pair(loop1*delta_length+loop2,2) += pi_test(loop1*delta_length+loop2,l2);
                      if(pi_test(loop1*delta_length+loop2,l1)<pi_test(loop1*delta_length+loop2,l2)){
                        numerator_pair(loop1*delta_length+loop2,2) += pi_test(loop1*delta_length+loop2,l2);
                      }

                    }else if((OT_test[l1]>t_temp1)&(OT_test[l1]<=(t_temp1+delta_temp1))&(DELTA_test[l1]==0)&(OT_test[l2]<=(t_temp1+delta_temp1))&(DELTA_test[l2]==0)&(OT_test[l1]<OT_test[l2])){

                      denominator_pair(loop1*delta_length+loop2,3) += (1-pi_test(loop1*delta_length+loop2,l1))*pi_test(loop1*delta_length+loop2,l2);
                      if(pi_test(loop1*delta_length+loop2,l1)<pi_test(loop1*delta_length+loop2,l2)){
                        numerator_pair(loop1*delta_length+loop2,3) += (1-pi_test(loop1*delta_length+loop2,l1))*pi_test(loop1*delta_length+loop2,l2);
                      }

                    }


                  }

                }
              }

              // for(int k=0; k<4; ++k){
              //   if(denominator_pair(i*delta_length+j,k)!=0){
              //     AUC(i*delta_length+j) += numerator_pair(i*delta_length+j,k)/(denominator_pair(i*delta_length+j,k)*S_MCMC);
              //   }
              // }

              AUC[loop1*delta_length+loop2] = (sum(numerator_pair(loop1*delta_length+loop2,_))/sum(denominator_pair(loop1*delta_length+loop2,_)) );






      }
    }
  
  NumericVector u_test_out_vec(t_length*delta_length*NU*N_test);
  NumericVector at_u_test_vec(t_length*delta_length*N_test);
  NumericVector all_u_test_vec(t_length*delta_length*N_test);
  NumericVector numerator_pair_vec(t_length*delta_length*4);
  NumericVector denominator_pair_vec(t_length*delta_length*4);
  NumericVector pi_test_vec(t_length*delta_length*N_test);
  //NumericVector ratio_u_vec(t_length*delta_length*N_test);
  
  for(int i=0; i<t_length; ++i){
    for(int j=0; j<delta_length; ++j){
      for(int k=0; k<NU; ++k){
        u_test_out_vec[Range((i*delta_length+j)*NU*N_test+k*N_test,(i*delta_length+j)*NU*N_test+(k+1)*N_test-1)] = u_test_out((i*delta_length+j)*NU+k,_);
        
      }
      at_u_test_vec[Range((i*delta_length+j)*N_test,(i*delta_length+j+1)*N_test-1)] = at_u_test(i*delta_length+j,_);
      all_u_test_vec[Range((i*delta_length+j)*N_test,(i*delta_length+j+1)*N_test-1)] = all_u_test(i*delta_length+j,_);
      //ratio_u_vec[Range((i*delta_length+j)*N_test,(i*delta_length+j+1)*N_test-1)] = ratio_u(i*delta_length+j,_);
      
      pi_test_vec[Range((i*delta_length+j)*N_test,(i*delta_length+j+1)*N_test-1)] = pi_test(i*delta_length+j,_);
      numerator_pair_vec[Range((i*delta_length+j)*4,(i*delta_length+j+1)*4-1)] = numerator_pair(i*delta_length+j,_);
      denominator_pair_vec[Range((i*delta_length+j)*4,(i*delta_length+j+1)*4-1)] = denominator_pair(i*delta_length+j,_);
    }
  }
  
  
  return List::create(Rcpp::Named("u_test_vec") = u_test_out_vec,
                      Rcpp::Named("at_u_test_vec") = at_u_test_vec,
                      Rcpp::Named("all_u_test_vec") = all_u_test_vec,
                      Rcpp::Named("numerator_vec") = numerator_pair_vec,
                      Rcpp::Named("denominator_vec") = denominator_pair_vec,
                      Rcpp::Named("AUC") = AUC,
                      Rcpp::Named("pi_test_vec") = pi_test_vec
                      ,Rcpp::Named("ratio_u") = ratio_u  
                        
  );
}
#endif
