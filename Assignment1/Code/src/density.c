#include <Rinternals.h>
#include <R.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>
#include <math.h>
#include <float.h>

SEXP trisum_cv(SEXP x, SEXP fn, SEXP rho, SEXP lambda) {
  
  int n = length(x);
  
  double *xx = REAL(x);
  double *h = REAL(lambda);
  
  SEXP K = PROTECT(allocVector(REALSXP, n));
  double *KK = REAL(K);
  memset(KK, 0, sizeof(double) * n);
  
  SEXP fn_call = PROTECT(lang2(fn, R_NilValue));
  SEXP s = PROTECT(ScalarReal(0.0));

  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < i; ++j) {
      SETCADR(fn_call, ScalarReal((xx[i] - xx[j]) / (*h)));
      s = eval(fn_call, rho);
      double tmp = REAL(s)[0];
      
      KK[i] += tmp;
      KK[j] += tmp;
    }
  }
  
  SEXP out = ScalarReal(0.0);
  double *oo = REAL(out);
  *oo = 0;
  
  for(int i = 0; i < n; ++i)
    if(KK[i] > DBL_EPSILON) {
      *oo += log(KK[i]);
    }
  
  *oo = n * log((n - 1) * (*h)) - (*oo);

  UNPROTECT(3);
  return out; 
}

