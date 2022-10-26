#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

double loglike(NumericVector &x,
               double nu,
               arma::vec::fixed<2> &parms) {
  
  int n = x.size();
  
  double K = sum(log(1 + (x - parms[0]) * (x - parms[0]) / (nu * parms[1])));
  
  return std::log(parms[1]) / 2 + (nu + 1) * K / (2 * n);
}

arma::vec::fixed<2> gradient(NumericVector &x,
                      double nu,
                      arma::vec::fixed<2> &parms) {
  
  int n = x.size();
  
  NumericVector C1 = (x - parms[0]) / (1 + (x - parms[0]) * (x - parms[0]) / (nu * parms[1]));
    
  double K_mu = sum(C1);
  double K_sigma = sum(C1 * (x - parms[0]));
  
  double grad_mu = -(nu + 1) * K_mu / (n * nu * parms[1]);
  double grad_sigma = 1 / (2 * parms[1]) - 
    (nu + 1) * K_sigma / (2 * n * nu * parms[1] * parms[1]);
  
  arma::vec::fixed<2> out;
  
  out(0) = grad_mu;
  out(1) = grad_sigma;
  
  return out;

}

arma::mat::fixed<2, 2> hessian(NumericVector &x, 
                       double nu,
                       arma::vec::fixed<2> &parms) {
  
  int n = x.size();
  
  NumericVector C0 = 1 / (1 + (x - parms[0]) * (x - parms[0]) / (nu * parms[1]));
  NumericVector C1 = C0 * (x - parms[0]);
  NumericVector C2 = C1 * (x - parms[0]);
  
  double hess_mu = (nu + 1) * sum(C0) / (n * nu * parms[1]) + 
    2 * (nu + 1) * sum(C1 * C1) / (n * (nu * parms[1]) * (nu * parms[1]));
  
  double hess_sigma = -1 / (2 * parms[1] * parms[1]) + 
    (nu + 1) * sum(C2) / (n * nu * parms[1] * parms[1] * parms[1]) -
    (nu + 1) * sum(C2 * C2) / (2 * n * nu * nu * parms[1] * parms[1] * parms[1] * parms[1]);
    
  double hess_mu_sigma = (nu + 1) * sum(C1) / (n * nu * parms[1] * parms[1]) -
    (nu + 1) * sum(C1 * C1 * (x - parms[0])) / (n * nu * parms[1] * parms[1] * parms[1]);
  
  arma::mat::fixed<2, 2> out;
  
  out(0, 0) = hess_mu;
  out(0, 1) = hess_mu_sigma;
  out(1, 0) = hess_mu_sigma;
  out(1, 1) = hess_sigma;
  
  return out;
}

double norm_diff(arma::vec::fixed<2> &x, arma::vec::fixed<2> &y) {
  return std::sqrt(x[0] * x[0] + y[0] * y[0] - 2 * x[0] * y[0] +
    x[1] * x[1] + y[1] * y[1] - 2 * x[1] * y[1]);
}

// [[Rcpp::export]]
NumericVector Newton_rcpp(NumericVector x,
                          NumericVector par0,
                          double nu,
                          double d = 0.8,
                          double c = 0.2,
                          double gamma0 = 1,
                          double epsilon = 1e-7,
                          int maxit = 500
                          ) {
  
  arma::vec::fixed<2> par = as<arma::vec>(par0);
  
  for(int i = 0; i < maxit; ++i) {
    double value = loglike(x, nu, par);
    auto grad = gradient(x, nu, par);
    auto hess = hessian(x, nu, par);
    
    arma::vec::fixed<2> rho = arma::solve(hess, grad);
    
    double gamma = gamma0;
    
    arma::vec::fixed<2> par1 = par + gamma * rho;
    
    double h_prime = arma::dot(grad, rho);
    
    while(loglike(x, nu, par1) > value + c * gamma * h_prime) {
      gamma *= d;
      par1 = par + gamma * rho;
    }
    
    if(norm_diff(par, par1) < epsilon * (std::sqrt(arma::sum(par % par)) + epsilon)) {
      par = par1;
      break;
    }
    
    par = par1;
  }
  
  return wrap(par);
}