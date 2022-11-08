#include <Rinternals.h>
#include <R.h>
#include <float.h>

SEXP C_cv(SEXP x, SEXP fn, SEXP lambda, SEXP rho) {
  if(REAL(lambda)[0] < DBL_EPSILON) return ScalarReal(INFINITY);
  int n = length(x);
  SEXP K = PROTECT(allocVector(REALSXP, n));
  SEXP fn_call = PROTECT(lang2(fn, R_NilValue));
  SEXP out = PROTECT(ScalarReal(0.0));
  double *x_ = REAL(x), *K_ = REAL(K), *out_ = REAL(out);
  double h_ = REAL(lambda)[0];
  memset(K_, 0, sizeof(double) * n);
  for(int i = 1; i < n; ++i) {
    SEXP tmp = PROTECT(allocVector(REALSXP, i));
    double *tmp_ = REAL(tmp);
    for(int k = 0; k < i; ++k)
      tmp_[k] = (x_[i] - x_[k]) / h_;
    SETCADR(fn_call, tmp);
    SEXP s = eval(fn_call, rho);
    double *s_ = REAL(s);
    UNPROTECT(1);
    for(int j = 0; j < i; ++j) {
      K_[i] += s_[j];
      K_[j] += s_[j];
    }
  }
  for(int i = 0; i < n; ++i)
    if(K_[i] > DBL_EPSILON) *out_ += log(K_[i]);
  *out_ = n * log((n - 1) * h_) - (*out_);
  UNPROTECT(3);
  return out;
}
/*
SEXP C_cv(SEXP x, SEXP fn, SEXP lambda, SEXP rho) {
  int n = length(x);
  SEXP K = PROTECT(allocVector(REALSXP, n));
  SEXP fn_call = PROTECT(lang2(fn, R_NilValue));
  SEXP out = PROTECT(ScalarReal(0.0));
  double *x_ = REAL(x), *K_ = REAL(K), *out_ = REAL(out);
  double h_ = REAL(lambda)[0];
  
  memset(K_, 0, sizeof(double) * n);
  
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < i; ++j) {
      SETCADR(fn_call, ScalarReal((x_[i] - x_[j]) / h_));
      SEXP s = eval(fn_call, rho);
      double tmp = REAL(s)[0];
      K_[i] += tmp;
      K_[j] += tmp;
    }
  }
  for(int i = 0; i < n; ++i)
    if(K_[i] > DBL_EPSILON) *out_ += log(K_[i]);
    
    *out_ = n * log((n - 1) * h_) - (*out_);
    UNPROTECT(3);
    return out; 
}*/
