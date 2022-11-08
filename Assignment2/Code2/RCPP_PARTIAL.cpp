#include <Rcpp.h>
#include <vector>
using namespace Rcpp;

double env_quantile(double x,
                    NumericVector &a,
                    NumericVector &b,
                    std::vector<double> &az,
                    std::vector<double> &Q) {
  int n = a.size();
  double c = Q[n];
  int maxi;
  for(maxi = 0; maxi < n + 1; ++maxi)
    if(c * x <= Q[maxi]) break;
  maxi -= 1;
  double y = c * x - Q[maxi];
  return std::log(a[maxi] * y * std::exp(-b[maxi]) + 
                  std::exp(az[maxi])) / a[maxi];
}

double env_density(double x,
                   NumericVector &a,
                   NumericVector &b,
                   NumericVector &z) {
  int n = a.size();
  if(x > z[n] || x < z[0])
    return 0;
  int maxi;
  for(maxi = 0; maxi < n + 1; ++maxi)
    if(x <= z[maxi]) break;
  maxi -= 1;
  return std::exp(a[maxi] * x + b[maxi]);
}

// [[Rcpp::export]]
List RCPP_adap_samp_partial(int n,
                                     Function density,
                                     NumericVector a,
                                     NumericVector b,
                                     NumericVector z) {
  int m = a.size();
  std::vector<double> Q(m + 1);
  std::vector<double> az(m);
  Q[0] = 0;
  for(int i = 1; i < m + 1; ++i) {
    az[i - 1] = a[i - 1] * z[i - 1];
    Q[i] = Q[i - 1] + 
      std::exp(b[i - 1]) * 
      (std::exp(a[i - 1] * z[i]) - std::exp(az[i - 1])) / a[i - 1];
  }
  NumericVector samples(n);
  int accepts = 0, tries = 0;
  for(int i = 0; i < n; ++i) {
    int reject = 1;
    while(reject == 1) {
      ++tries;
      double u0 = R::runif(0, 1);
      double u1 = R::runif(0, 1);
      double y0 = env_quantile(u0, a, b, az, Q);
      double env_y0 = env_density(y0, a, b, z);
      NumericVector dens_y0 = density(y0);
      if(u1 <= dens_y0[0] / env_y0) {
        reject = 0;
        samples[i] = y0;
        ++accepts;
      }
    }
  }
  NumericVector rate(1);
  rate[0] = ((double) tries - (double) accepts) / (double) tries;
  return List::create(samples, rate);
}