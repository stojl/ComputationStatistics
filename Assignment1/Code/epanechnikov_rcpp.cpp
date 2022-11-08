#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector e_kernel_cpp(NumericVector x) {
  int n = x.size();
  NumericVector result(n);
  
  for(int i = 0; i < n; ++i) {
    result[i] = std::abs(x[i]) <= 1 ? 0.75 * (1 - x[i] * x[i]) : 0;
  }
  
  return result;
}

