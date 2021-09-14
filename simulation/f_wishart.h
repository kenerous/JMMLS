#ifndef _F_WISHART_H
#define _F_WISHART_H

#include <RcppArmadillo.h>
using namespace arma; 
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]



NumericMatrix f_wishart(int df, NumericMatrix sigma) {
  arma::mat sigma_temp= Rcpp::as<arma::mat>(sigma);
  int ncols = sigma_temp.n_cols;
  mat Y = randn(df, ncols);
  mat x = (Y * chol(sigma_temp)).t();
  return(wrap(x * x.t()));
}

#endif
