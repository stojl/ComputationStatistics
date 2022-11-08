#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// [[Rcpp::export]]
double CPP_likelihood(NumericVector par, NumericVector x, double nu) {
  double K = sum(log(1 + (x - par[0]) * (x - par[0]) / (nu * par[1])));
  return log(par[1]) / 2 + (nu + 1) * K / (2 * x.size());
}

// [[Rcpp::export]]
NumericVector CPP_gradient(NumericVector par, NumericVector x, double nu) {
  int n = x.size();
  double K_mu = 0, K_sigma = 0;
  for(int i = 0; i < n; ++i) {
    double C1 = (x[i] - par[0]) / 
      (1 + (x[i] - par[0]) * (x[i] - par[0]) / (nu * par[1]));
    K_mu += C1;
    K_sigma += C1 * (x[i] - par[0]);
  }
  NumericVector grad(2);
  grad[0] = -(nu + 1) * K_mu / (n * nu * par[1]);
  grad[1] = 1 / (2 * par[1]) - 
    (nu + 1) * K_sigma / (2 * n * nu * par[1] * par[1]);
  return grad;
}

// [[Rcpp::export]]
NumericMatrix CPP_hessian(NumericVector par, NumericVector x, double nu) {
  int n = x.size();
  double K1 = 0, K2 = 0, K3 = 0, K4 = 0, K5 = 0, K6 = 0;
  for(int i = 0; i < n; ++i) {
    double C0 = 1 / (1 + (x[i] - par[0]) * (x[i] - par[0]) / (nu * par[1]));
    double C1 = C0 * (x[i] - par[0]);
    double C2 = C1 * (x[i] - par[0]);
    K1 += C0;
    K2 += C1;
    K3 += C1 * C1;
    K4 += C2;
    K5 += C2 * C2;
    K6 += C1 * C1 * (x[i] - par[0]);
  }
  NumericMatrix hess(2, 2);
  hess(0, 0) = (nu + 1) * K1 / (n * nu * par[1]) +
    2 * (nu + 1) * K3 / (n * nu * nu * par[1] * par[1]);
  hess(0, 1) = (nu + 1) * K2 / (n * nu * par[1] * par[1]) -
    (nu + 1) * K6 / (n * nu * par[1] * par[1] * par[1]);
  hess(1, 0) = hess(0, 1);
  hess(1, 1) = -1 / (2 * par[1] * par[1]) +
    (nu + 1) * K4 / (n * nu * par[1] * par[1] * par[1]) -
    (nu + 1) * K5 / (2 * n * nu * nu * par[1] * par[1] * par[1] * par[1]);
  return hess;
}
// [[Rcpp::export]]
List FitT(NumericVector par0, NumericVector x, double nu, double c,
          double d, double gamma0, int maxit, double eps) {
  
  arma::mat::fixed<2, 2> hess;
  arma::vec::fixed<2> grad;
  arma::vec::fixed<2> p1;
  arma::vec::fixed<2> p2;
  NumericVector par1(2), par = clone(par0);
  int n = x.size();
  int i;
  for(i = 0; i < maxit; ++i) {
    double value = CPP_likelihood(par, x, nu);
    NumericVector gr = CPP_gradient(par, x, nu);
    NumericMatrix hessian = CPP_hessian(par, x, nu);
    grad[0] = gr[0], grad[1] = gr[1];
    hess(0, 0) = hessian(0, 0), hess(0, 1) = hessian(0, 1);
    hess(1, 0) = hessian(1, 0), hess(1, 1) = hessian(1, 1);
    arma::vec::fixed<2> rho = -arma::solve(hess, grad);
    double gamma = gamma0;
    par1[0] = par[0] + gamma * rho[0];
    par1[1] = par[0] + gamma * rho[1];
    double h_prime = rho[0] * gr[0] + rho[1] * gr[1];
    while(CPP_likelihood(par1, x, nu) > value + c * gamma * h_prime) {
      gamma *= d;
      par1[0] = par[0] + gamma * rho[0];
      par1[1] = par[1] + gamma * rho[1];
    }
    p1[0] = par1[0], p1[1] = par1[1];
    p2[0] = par1[0] - par[0], p2[0] = par1[1] - par[1];
    double norm = arma::norm(p1, 2);
    double norm_new = arma::norm(p2, 2);
    if(norm_new < eps * (norm + eps)) break;
    par[0] = par1[0], par[1] = par1[1];
  }
  if(i != maxit) i -= 1;
  NumericVector iter(1);
  iter[0] = i;
  return List::create(par1, iter);
}
