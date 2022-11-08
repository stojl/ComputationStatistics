#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector RCPP_poisdens(NumericVector y, NumericVector x, NumericVector z) {
  int n = y.size();
  NumericVector expyx(n);
  for(int i = 0; i < n; ++i)
    expyx[i] = sum(exp(y[i] * x));
  return exp(y * sum(x * z) - expyx);
}

// [[Rcpp::export]]
NumericVector RCPP_poisdens_derv(NumericVector y, 
                                 NumericVector x, 
                                 NumericVector z) {
  int n = y.size();
  NumericVector expyx(n);
  NumericVector x_expyx(n);
  for(int i = 0; i < n; ++i) {
    expyx[i] = sum(exp(y[i] * x));
    x_expyx[i] = sum(x * exp(y[i] * x));
  }
  double xz = sum(x * z);
  return exp(y * xz - expyx) * (xz - x_expyx);
}

// [[Rcpp::export]]
double RCPP_env_density(double x, 
                        NumericVector a, 
                        NumericVector b, 
                        NumericVector z) {
  int m = z.size();
  if(x < z[0] || x > z[m - 1]) return 0;
  int maxi;
  for(maxi = 1; maxi < m; ++maxi)
    if(x <= z[maxi]) break;
  maxi -= 1;
  return std::exp(a[maxi] * x + b[maxi]);
}

// [[Rcpp::export]]
double RCPP_env_quantile(double x, 
                         NumericVector a, 
                         NumericVector b, 
                         NumericVector z,
                         NumericVector az,
                         NumericVector Q) {
  int m = Q.size(), maxi;
  double c = Q[m - 1];
  for(maxi = 0; maxi < m; ++maxi)
    if(c * x <= Q[maxi]) break;
  maxi -= 1;
  double y = c * x - Q[maxi];
  return std::log(a[maxi] * y * 
                  std::exp(-b[maxi]) + std::exp(az[maxi])) / a[maxi];
}