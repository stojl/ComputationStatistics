#include <Rcpp.h>
using namespace Rcpp;

double kernel(double x) {
  if(std::abs(x) > 1) {
    return 0;
  } else {
    return 0.75 * (1 - x * x);
  }
}

// [[Rcpp::export]]
NumericVector trisum(NumericVector x, double lambda) {
  int n = x.size();
  NumericVector K(n);
  
  for(int i = 1; i < n; ++i) {
    for(int j = 0; j <= i; ++j) {
      double tmp = kernel((x[i] - x[j]) / lambda);
      K[i] += tmp;
      K[j] += tmp;
    }
  }
  
  return K;
}

// [[Rcpp::export]]
double trisum_cv(NumericVector x, double lambda) {
  int n = x.size();
  NumericVector K(n);
  
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < i; ++j) {
      double tmp = kernel((x[i] - x[j]) / lambda);
      K[i] += tmp;
      K[j] += tmp;
    }
  }
  
  return n * log((n - 1) * lambda) - sum(log(K)); 
}