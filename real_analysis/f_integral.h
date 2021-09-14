#ifndef _F_INTEGRAL_H
#define _F_INTEGRAL_H


#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double f_integral(int IF_QUADRATIC, int NX, NumericVector x, NumericVector alpha, NumericVector kns, NumericVector u, NumericVector J, NumericVector lambda, int Fail_label, int K, int NU, double OT, double c, int start, int end) {
  double integral =0;
  NumericVector second_level(NU);
  for(int j=0; j<NU; ++j){
    for(int k=0; k<(NX+1); ++k){
    	second_level[j] += alpha[j*(NX+1)+k]*x[k];	
	}       	
  }
  
  for(int i=start; i<end; ++i){

    if(i==(Fail_label-1)){
      for(int j=0; j<K; ++j){

        double t = J[i]+(OT-J[i])*j/K;
        double temp = u[0]+second_level[0];
        if(NU>1){
          temp += (u[1]+second_level[1])*t;
          if(NU>2){
            if(IF_QUADRATIC==1){
              temp += (u[2]+second_level[2])*pow(t,2);
            }else{
              for(int k=2; k<NU; ++k){
                temp += (u[k]+second_level[k])*(-1)*((std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[k-2],3))+pow(t-kns[k-2],3)))/(2*(kns[NU-1]-kns[k-2]))-(std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[NU-2],3))+pow(t-kns[NU-2],3)))/(2*(kns[NU-1]-kns[NU-2])));
              }
            }
          }
        }
        
        integral += lambda[i]*exp(c*temp)*(OT-J[i])/K;
      }
    }else{
      for(int j=0; j<K; ++j){
        double t = J[i]+(J[i+1]-J[i])*j/K;
        double temp = u[0]+second_level[0];
        if(NU>1){
          temp += (u[1]+second_level[1])*t;
          if(NU>2){
            if(IF_QUADRATIC==1){
              temp += (u[2]+second_level[2])*pow(t,2);
            }else{
              for(int k=2; k<NU; ++k){
                temp += (u[k]+second_level[k])*(-1)*((std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[k-2],3))+pow(t-kns[k-2],3)))/(2*(kns[NU-1]-kns[k-2]))-(std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[NU-2],3))+pow(t-kns[NU-2],3)))/(2*(kns[NU-1]-kns[NU-2])));
              }
            }
          }
        }

        integral += lambda[i]*exp(c*temp)*(J[i+1]-J[i])/K;
      }
    }

  }
  
  return integral;  
}

#endif
