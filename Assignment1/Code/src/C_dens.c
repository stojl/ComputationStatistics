#include <Rinternals.h>
#include <R.h>

SEXP C_dens(SEXP x, SEXP p, SEXP kernel, SEXP bw, SEXP rho) {
  int n = length(x), m = length(p); 
  SEXP dens = PROTECT(allocVector(REALSXP, m));
  SEXP tmp = PROTECT(allocVector(REALSXP, n));
  SEXP K_Call = PROTECT(lang2(kernel, R_NilValue));
  double *x_ = REAL(x), *p_ = REAL(p), *tmp_ = REAL(tmp), *dens_ = REAL(dens);
  double bw_ = REAL(bw)[0];
  
  memset(dens_, 0, sizeof(double) * m);
  
  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j)
      tmp_[j] = (p_[i] - x_[j]) / bw_;
    SETCADR(K_Call, tmp);
    SEXP result = eval(K_Call, rho);
    double *result_ = REAL(result);
    for(int j = 0; j < n; ++j)
      dens_[i] += result_[j];
    dens_[i] /= n * bw_;
  }
  UNPROTECT(3);
  return dens;
}
/*
SEXP C_dens(SEXP x, SEXP p, SEXP kernel, SEXP bw, SEXP rho) {
  int n = length(x);
  int m = length(p);
  SEXP dens = PROTECT(allocVector(REALSXP, m));
  SEXP K_Call = PROTECT(lang2(kernel, R_NilValue));
  double *x_ = REAL(x);
  double *p_ = REAL(p);
  double *dens_ = REAL(dens);
  double bw_ = REAL(bw)[0];
  
  memset(dens_, 0, sizeof(double) * m);
  
  for(int i = 0; i < m; ++i) {
    for(int j = 0; j < n; ++j) {
      SETCADR(K_Call, ScalarReal((p_[i] - x_[j]) / bw_));
      SEXP result = eval(K_Call, rho);
      dens_[i] += REAL(result)[0];
    }
    dens_[i] /= n * bw_;
  }
  UNPROTECT(2);
  return dens;
}*/