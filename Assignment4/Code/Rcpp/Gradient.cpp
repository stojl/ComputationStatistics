#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
NumericVector gradient_rcpp(NumericVector par,
                            NumericVector x,
                            NumericVector y) {
  
  int N = x.size();
  NumericVector gr(4);
  
  for(int i = 0; i < N; ++i) {
    double elogx, da, db, dg, dr, yf, logx;
    logx = std::log(x[i]);
    elogx = std::exp(par[1] * logx - par[0]);
    dr = 1 / (1 + elogx);
    dg = 1 - dr;
    da = elogx * (par[3] - par[2]) * dr * dr;
    db = -da * logx;
    yf = y[i] - par[2] - (par[3] - par[2]) * dr;
    gr[0] -= da * yf / 2;
    gr[1] -= db * yf / 2;
    gr[2] -= dg * yf / 2;
    gr[3] -= dr * yf / 2;
  }
  
  return gr / N;
  
}

// [[Rcpp::export]]
NumericVector epoch_rcpp(NumericVector par0, 
                    NumericVector x, 
                    NumericVector y,
                    int minisize,
                    double gamma) {
  int n, MAX;
  Rcpp::Range r(0, minisize - 1);
  
  n = x.size();
  MAX = std::floor(n / minisize);
  NumericVector par = clone(par0);
  
  for(int i = 0; i < MAX; ++i, r + minisize) {
     par = par - gamma * gradient_rcpp(par, x[r], y[r]);
  }

  return par;
}

// [[Rcpp::export]]
NumericVector epoch_rcpp_partial(NumericVector par0, 
                         NumericVector x, 
                         NumericVector y,
                         Function gr,
                         int minisize,
                         double gamma) {
  
  int n = x.size();
  int MAX = std::floor(n / minisize);
  NumericVector par = clone(par0);
  
  Rcpp::Range r(0, minisize - 1);
  
  for(int i = 0; i < MAX; ++i, r + minisize) {
    par = par - gamma * as<NumericVector>(gr(par, x[r], y[r]));
  }
  
  return par;
}