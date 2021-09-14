#ifndef _F_DET_H
#define _F_DET_H
#include <RcppEigen.h>

using namespace Eigen;
using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
double f_det(NumericMatrix AA){
	const Eigen::Map<Eigen::MatrixXd> AA1(as<Eigen::Map<Eigen::MatrixXd> >(AA));
     return  AA1.determinant();
}

#endif
