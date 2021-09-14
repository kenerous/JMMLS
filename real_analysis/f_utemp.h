#ifndef _F_UTEMP_H
#define _F_UTEMP_H

#include "f_mmult.h"
#include "f_madd.h" 

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix f_utemp(NumericVector x, NumericVector alpha, int NX, NumericVector kns, NumericVector J, NumericVector lambda, int Fail_label, int K, int NU, double OT, double c, int start, int end) {
  NumericMatrix integral(NU,NU);
  NumericMatrix spline1(NU,1);
  NumericMatrix spline2(1,NU);
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
        double temp = second_level[0] + second_level[1]*t;
        spline1(0,0) = 1;
        spline2(0,0) = 1;
        spline1(1,0) = t;
        spline2(0,1) = t;
        for(int k=2; k<NU; ++k){
          spline1(k,0) = (-1)*((std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[k-2],3))+pow(t-kns[k-2],3)))/(2*(kns[NU-1]-kns[k-2]))-(std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[NU-2],3))+pow(t-kns[NU-2],3)))/(2*(kns[NU-1]-kns[NU-2])));
          spline2(0,k) = spline1(k,0);
          temp += second_level[k]*spline1(k,0);
        }
        integral = f_madd(integral, f_mmult(spline1,spline2)*exp(c*temp)*lambda[i]*pow(c,2)*(OT-J[i])/K);
      }
    }else{
      for(int j=0; j<K; ++j){
        double t = J[i]+(J[i+1]-J[i])*j/K;
        double temp = second_level[0] + second_level[1]*t;
        spline1(0,0) = 1;
        spline2(0,0) = 1;
        spline1(1,0) = t;
        spline2(0,1) = t;
        for(int k=2; k<NU; ++k){
          spline1(k,0) = (-1)*((std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[k-2],3))+pow(t-kns[k-2],3)))/(2*(kns[NU-1]-kns[k-2]))-(std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[NU-2],3))+pow(t-kns[NU-2],3)))/(2*(kns[NU-1]-kns[NU-2])));
          spline2(0,k) = spline1(k,0);
          temp += second_level[k]*spline1(k,0);
        }
        integral = f_madd(integral, f_mmult(spline1,spline2)*exp(c*temp)*lambda[i]*pow(c,2)*(J[i+1]-J[i])/K);
      }
    }

  }
  
  return integral;  
}


#endif
