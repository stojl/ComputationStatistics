#include <Rcpp.h>

using namespace Rcpp;

double inline norm_diff(NumericVector &x, NumericVector &y) {
  return std::sqrt(x[0] * x[0] + y[0] * y[0] - 2 * x[0] * y[0] +
                   x[1] * x[1] + y[1] * y[1] - 2 * x[1] * y[1]);
}

// [[Rcpp::export]]
int EM_rcpp(NumericVector x,
                      double nu,
                      NumericVector parm,
                      double epsilon = 1e-7,
                      int maxit = 15) {
  
  NumericVector par = clone(parm);
  NumericVector par1(2);
  
  NumericVector EW(x.size());
  
  int k = 0;
  
  for(int i = 0; i < maxit; ++i) {
    k = i + 1;
    EW = (nu + 1) / (1 + (x - par[0]) * (x - par[0]) / (nu * par[1]));
    
    par1[0] = sum(EW * x) / sum(EW);
    par1[1] = mean(EW * pow(x - par1[0], 2)) / nu;
    
    if(norm_diff(par, par1) < epsilon * (std::sqrt(sum(par * par)) + epsilon)) {
      break;
    }
    
    par[0] = par1[0];
    par[1] = par1[1];
    
    
    
  }
  
  NumericVector out(2);
  out[0] = par1[0];
  out[1] = par1[1];
  return k;
}
/*
NumericVector EM_rcpp(NumericVector x,
                      double nu,
                      NumericVector parm,
                      long double epsilon = 1e-7,
                      int maxit = 15) {
  
  double mu = parm[0];
  double sigma = parm[1];
  int k = -1;
  
  long double par[2];
  long double par1[2];
  long double diff[2];
  
  par[0] = parm[0];
  par[1] = parm[1];
  
  NumericVector EW(x.size());
  
  for(int i = 0; i < maxit; ++i) {
    EW = (nu + 1) / (1 + (x - par[0]) * (x - par[0]) / (nu * par[1]));
    
    par1[0] = sum(EW * x) / sum(EW);
    par1[1] = mean(EW * pow(x - par1[0], 2)) / nu;
    
    long double norm = std::sqrt((par[0] - par1[0]) * (par[0] - par1[0]) +
      (par[1] - par1[1]) * (par[1] - par1[1]));
    
    long double par_norm = std::sqrt(par[0] * par[0] + par[1] * par[1]);
    
    if(norm < epsilon * (par_norm + epsilon)) {
      break;
    }
    
    par[0] = par1[0];
    par[1] = par1[1];
    
  }
  
  NumericVector out(2);
  out[0] = par1[0];
  out[1] = par1[1];
  return out;
}
*/