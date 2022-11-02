#include <Rinternals.h>
#include <R.h>
#include <R_ext/Rdynload.h>
#include <math.h>

SEXP epoch_batch(SEXP par0, SEXP index, SEXP loss_gr, SEXP gamma, SEXP mbs, SEXP rho) {
  int mbs2, M, n;
  n = length(par0);
  if(INTEGER(mbs)[0] > length(index)) {
    mbs2 = length(index);
  } else {
    mbs2 = INTEGER(mbs)[0];
  }
  M = floor(length(index) / mbs2);
  SEXP gr_call = PROTECT(lang3(loss_gr, R_NilValue, R_NilValue));
  SEXP par = PROTECT(allocVector(REALSXP, n));
  SEXP m_index = PROTECT(allocVector(INTSXP, mbs2));
  int *iptr = INTEGER(index);
  int *miptr = INTEGER(m_index);
  SEXP gr;
  for(int i = 0; i < n; ++i)
    REAL(par)[i] = REAL(par0)[i];
  for(int i = 0; i < M; ++i) {
    for(int j = 0; j < mbs2; ++j) {
      miptr[j] = iptr[i * mbs2 + j];
    }
    SETCADR(gr_call, par);
    SETCADDR(gr_call, m_index);
    gr = eval(gr_call, rho);
    for(int k = 0; k < n; ++k)
      REAL(par)[k] = REAL(par)[k] - REAL(gamma)[0] * REAL(gr)[k];
  }
  UNPROTECT(3);
  return par;
}