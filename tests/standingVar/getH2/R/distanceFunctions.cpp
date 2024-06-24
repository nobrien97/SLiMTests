#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]

// Frobenius norm
double EucNorm(const mat& M) {
  return norm(M, "fro");
}

// When alpha = 0, log euclidean distance
double distLogEuclidean(const mat& P1, const mat& P2) {
  vec eigval1;
  mat eigvec1;
  eig_sym(eigval1, eigvec1, P1);  // Eigen decomposition of P1
  
  mat logP1 = eigvec1 * diagmat(log(eigval1)) * eigvec1.t();
  
  vec eigval2;
  mat eigvec2;
  eig_sym(eigval2, eigvec2, P2);  // Eigen decomposition of P2
  
  mat logP2 = eigvec2 * diagmat(log(eigval2)) * eigvec2.t();
  
  return EucNorm(logP1 - logP2);
}


// Euclidean-Power distance
double distPowerEuclidean(const mat& P1, const mat& P2, double alpha = 0.5) {
  if (alpha != 0) {
    vec eigval1;
    mat eigvec1;
    eig_sym(eigval1, eigvec1, P1);  // Eigen decomposition of P1
    
    mat Q1 = eigvec1 * diagmat(pow(abs(eigval1), alpha)) * eigvec1.t();
    
    vec eigval2;
    mat eigvec2;
    eig_sym(eigval2, eigvec2, P2);  // Eigen decomposition of P2
    
    mat Q2 = eigvec2 * diagmat(pow(abs(eigval2), alpha)) * eigvec2.t();
    return EucNorm(Q1 - Q2) / alpha;
    }
  
  return distLogEuclidean(P1, P2);
}

// [[Rcpp::export]]
NumericMatrix distanceMatrix(List covarianceMatrices, double alpha = 0.5) {
  int n = covarianceMatrices.size();
  
  // make sure alpha > 0: for some reason abs() causes errors
  if (alpha < 0) {
    alpha = -alpha;
  }

  // result matrix
  NumericMatrix distances(n,n);

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      mat P1 = covarianceMatrices[i];
      mat P2 = covarianceMatrices[j];
      double dist = distPowerEuclidean(P1, P2, alpha);
      distances(i, j) = dist;
      distances(j, i) = dist;
    }
  }
  
  return distances;
}
