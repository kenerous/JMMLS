#ifndef _F_INVGAUSS_H
#define _F_INVGAUSS_H

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector f_invgauss(int n, double mu, double lambda){

	NumericVector random_vector(n);
	double z,y,x,u;

	for(int i=0; i<n; ++i){
		z=rnorm(1,0,1)[0];
		y=z*z;
		x=mu+0.5*mu*mu*y/lambda - 0.5*(mu/lambda)*sqrt(4*mu*lambda*y+mu*mu*y*y);
		u=runif(1,0,1)[0];
		if(u <= mu/(mu+x)){
			random_vector(i)=x;
		}else{
			random_vector(i)=mu*mu/x;
		};
	}
	return(random_vector);
}

#endif
