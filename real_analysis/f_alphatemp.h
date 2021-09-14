#ifndef _F_ALPHATEMP_H
#define _F_ALPHATEMP_H

#include "f_mmult.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix f_alphatemp(NumericVector x, NumericVector u, int NX, NumericVector kns, NumericVector J, NumericVector lambda, int Fail_label, int K, int NU, double OT, double c, int start, int end) {
  NumericMatrix integral(NU*(NX+1),NU*(NX+1));
  NumericMatrix spline1(NU*(NX+1),1);
  NumericMatrix spline2(1,NU*(NX+1));
  NumericVector second_level(NU*(NX+1));
  
  
  for(int i=start; i<end; ++i){

    if(i==(Fail_label-1)){
      for(int j=0; j<K; ++j){
        double t = J[i]+(OT-J[i])*j/K;
        double temp = u[0] + u[1]*t;
        for(int k=0; k<(NX+1); ++k){
        	spline1(k,0) = x[k];
        	spline2(0,k) = spline1(k,0);
        	spline1(NX+1+k,0) = t*x[NX+1+k];
        	spline2(0,NX+1+k) = spline1(0,NX+1+k);
		}
        for(int k=2; k<NU; ++k){
        	for(int l=0; l<(NX+1); ++l){
				spline1(k*(NX+1)+l,0) = x[k*(NX+1)+l]*(std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[k-2],3))+pow(t-kns[k-2],3)))/(2*(kns[NU-1]-kns[k-2]))-(std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[NU-2],3))+pow(t-kns[NU-2],3)))/(2*(kns[NU-1]-kns[NU-2]));
          		spline2(0,k*(NX+1)+l) = spline1(k*(NX+1)+l,0);        		
			}
            temp += u[k]*spline1(k,0);
        }
        integral += f_mmult(spline1,spline2)*exp(c*temp)*lambda[i]*pow(c,2)*(OT-J[i])/K;
      }
    }else{
      for(int j=0; j<K; ++j){
        double t = J[i]+(J[i+1]-J[i])*j/K;
        double temp = u[0] + u[1]*t;
        for(int k=0; k<(NX+1); ++k){
        	spline1(k,0) = x[k];
        	spline2(0,k) = spline1(k,0);
        	spline1(NX+1+k,0) = t*x[NX+1+k];
        	spline2(0,NX+1+k) = spline1(0,NX+1+k);
		}
        for(int k=2; k<NU; ++k){
        	for(int l=0; l<(NX+1); ++l){
				spline1(k*(NX+1)+l,0) = x[k*(NX+1)+l]*(std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[k-2],3))+pow(t-kns[k-2],3)))/(2*(kns[NU-1]-kns[k-2]))-(std::abs(pow(t-kns[NU-1],3))+pow(t-kns[NU-1],3)-(std::abs(pow(t-kns[NU-2],3))+pow(t-kns[NU-2],3)))/(2*(kns[NU-1]-kns[NU-2]));
          		spline2(0,k*(NX+1)+l) = spline1(k*(NX+1)+l,0);        		
			}
            temp += u[k]*spline1(k,0);
        }
        integral += f_mmult(spline1,spline2)*exp(c*temp)*lambda[i]*pow(c,2)*(J[i+1]-J[i])/K;
      }
    }

  }
  
  return integral;  
}


#endif
