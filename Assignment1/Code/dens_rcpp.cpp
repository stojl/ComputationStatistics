#include <Rcpp.h>
#include <limits>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector e_kernel_cpp(NumericVector x) {
  int n = x.size();
  NumericVector result(n);
  for(int i = 0; i < n; ++i)
    result[i] = std::abs(x[i]) <= 1 ? 0.75 * (1 - x[i] * x[i]) : 0;
  return result;
}


// [[Rcpp::export]]
NumericVector dens_rcpp_partial(NumericVector x,
                                Function kernel,
                                double bandwidth,
                                NumericVector points) {
  int n = x.size(), m = points.size();
  NumericVector result(m);
  for(int i = 0; i < m; ++i) {
    NumericVector call = kernel((points[i] - x) / bandwidth);
    result[i] = sum(call) / (n * bandwidth);
  }
  return result;
}

// [[Rcpp::export]]
double bw_cv_rcpp_partial(NumericVector x,
                          Function kernel,
                          double bandwidth) {
  int n = x.size();
  NumericVector K(n);
  double result;
  for(int i = 1; i < n; ++i) {
    Range r = Range(0, i - 1);
    NumericVector tmp = kernel((x[r] - x[i]) / bandwidth);
    for(int j = 0; j < i; ++j)
      K[j] += tmp[j], K[i] += tmp[j];
      
  }
  for(int s = 0; s < n; ++s)
    if(K[s] > 0) result += std::log(K[s]);
  return n * log((n - 1) * bandwidth) - result;
}

// [[Rcpp::export]]
NumericVector dens_rcpp(NumericVector x,
                        double bandwidth,
                        NumericVector points) {
  int n = x.size(), m = points.size();
  NumericVector result(m);
  for(int i = 0; i < m; ++i) {
    NumericVector call = e_kernel_cpp((points[i] - x) / bandwidth);
    result[i] = sum(call) / (n * bandwidth);
  }
  return result;
}

// [[Rcpp::export]]
double bw_cv_rcpp(NumericVector x,
                  double bandwidth) {
  int n = x.size();
  NumericVector K(n);
  double result;
  for(int i = 1; i < n; ++i) {
    Range r = Range(0, i - 1);
    NumericVector tmp = e_kernel_cpp((x[i] - x[r]) / bandwidth);
    for(int j = 0; j < i; ++j)
      K[j] += tmp[j], K[i] += tmp[j];
  }
  for(int s = 0; s < n; ++s)
    if(K[s] > 0) result += std::log(K[s]);
  return n * std::log((n - 1) * bandwidth) - result;
}




