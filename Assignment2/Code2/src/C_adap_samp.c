#include <R.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Random.h> // ACCESS TO UNIF NUMBER GENERATOR

double env_quantile(double x,
                    double *a,
                    double *b,
                    double *az,
                    double *Q,
                    int n) {
  double c = Q[n];
  int maxi;
  for(maxi = 0; maxi < n + 1; ++maxi)
    if(c * x <= Q[maxi]) break;
  maxi -= 1;
  double y = c * x - Q[maxi];
  return log(a[maxi] * y * exp(-b[maxi]) + exp(az[maxi])) / a[maxi];
}

double env_density(double x,
                   double *a,
                   double *b,
                   double *z,
                   int n) {
  if(x > z[n] || x < z[0])
    return 0;
  int maxi;
  for(maxi = 0; maxi < n + 1; ++maxi)
    if(x <= z[maxi]) break;
  maxi -= 1;
  return exp(a[maxi] * x + b[maxi]);
}

SEXP C_adap_samp(SEXP n,
                 SEXP density,
                 SEXP a,
                 SEXP b,
                 SEXP z,
                 SEXP rho) {
  int m = length(a);
  int N = INTEGER(n)[0];
  double *Q = (double *)malloc(sizeof(double) * (m + 1));
  double *az = (double *)malloc(sizeof(double) * m);
  double *a_ = REAL(a), *b_ = REAL(b), *z_ = REAL(z);

  Q[0] = 0;
  for(int i = 1; i < m + 1; ++i) {
    az[i - 1] = a_[i - 1] * z_[i - 1];
    Q[i] = Q[i - 1] + 
      exp(b_[i - 1]) * (exp(a_[i - 1] * z_[i]) - exp(az[i - 1])) / a_[i - 1];
  }
  
  SEXP density_call = PROTECT(lang2(density, R_NilValue));
  SEXP samples = PROTECT(allocVector(REALSXP, N));
  double *samples_ = REAL(samples);
  int accepts = 0, tries = 0;
  GetRNGstate();
  for(int i = 0; i < N; ++i) {
    int reject = 1;
    while(reject == 1) {
      ++tries;
      double u0 = unif_rand();
      double u1 = unif_rand();
      double y0 = env_quantile(u0, a_, b_, az, Q, m);
      double env_y0 = env_density(y0, a_, b_, z_, m);
      SETCADR(density_call, ScalarReal(y0));
      SEXP dens_y0 = eval(density_call, rho);
      if(u1 <= REAL(dens_y0)[0] / env_y0) {
        reject = 0;
        samples_[i] = y0;
        ++accepts;
      }
    }
  }
  PutRNGstate();
  SEXP values = PROTECT(allocVector(VECSXP, 2));
  double rate = ((double) tries - (double) accepts) / (double) tries;
  SET_VECTOR_ELT(values, 0, samples);
  SET_VECTOR_ELT(values, 1, ScalarReal(rate));
  UNPROTECT(3);
  free(Q);
  free(az);
  return values;
}