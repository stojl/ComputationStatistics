#include <Rcpp.h>
using namespace Rcpp;

double dist(double x1, double x2, double y1, double y2) {
  double x_sq = x1 * x1 + x2 * x2;
  double y_sq = y1 * y1 + y2 * y2;
  double xy = x1 * y1 + x2 * y2;
  return std::sqrt(x_sq + y_sq - 2 * xy);
}

// [[Rcpp::export]]
List CPP_EM(NumericVector par0, 
                    NumericVector x, 
                    double nu,
                    int maxit,
                    double eps) {
  int n = x.size();
  NumericVector par = clone(par0);
  NumericVector par1(2);
  int i;
  for(i = 0; i < maxit; ++i) {
    NumericVector EW = (nu + 1) / 
      (1 + (x - par[0]) * (x - par[0]) / (nu * par[1]));
    par1[0] = sum(EW * x) / sum(EW);
    par1[1] = sum(EW * (x - par[0]) * (x - par[0])) / (n * nu);
    double norm_new = dist(par[0], par[1], par1[0], par1[1]);
    double norm_old = std::sqrt(par1[0] * par1[0] + par1[1] * par1[1]);
    if(norm_new < eps * (norm_old + eps)) break;
    par[0] = par1[0];
    par[1] = par1[1];
  }
  if(i != maxit) i += 1;
  return List::create(par1, i);
}

// [[Rcpp::export]]
NumericVector phi_cpp(NumericVector par0,
                      NumericVector x,
                      double nu)  {
  int n = x.size();
  NumericVector par(2);
  NumericVector EW = (nu + 1) / 
    (1 + (x - par0[0]) * (x - par0[0]) / (nu * par0[1]));
  
  par[0] = sum(EW * x) / sum(EW);
  par[1] = sum(EW * (x - par[0]) * (x - par[0])) / (n * nu);
  
  return par;
}

// [[Rcpp::export]]
double Q_cpp(NumericVector par,
              NumericVector par1,
              NumericVector x,
              double nu) {
  
  int n = x.size();
  double value = 0.5 * log(par[1]) * n;
  NumericVector C2 = (nu + 1) / 
    (1 + (x - par1[0]) * (x - par1[0]) / (nu * par1[1]));
  return value + sum(C2 * (x - par[0]) * (x - par[0]) / (2 * nu * par[1]));
}