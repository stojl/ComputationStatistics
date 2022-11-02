#include <Rinternals.h>
#include <R.h>

SEXP sgd(SEXP par0, SEXP loss_gr, SEXP N, SEXP batch, SEXP epoch, SEXP gamma0, SEXP maxit, SEXP rho) {
  int n_maxit = asInteger(maxit);
  double *gam = REAL(gamma0);
  SEXP par = PROTECT(allocVector(REALSXP, length(par0)));
  for(int i = 0; i < length(par0); ++i)
    REAL(par)[i] = REAL(par0)[i];
  SEXP batch_call = PROTECT(lang2(batch, R_NilValue));
  SEXP s, t;
  for(int i = 0; i < n_maxit; ++i) {
    SETCADR(batch_call, N);
    SEXP index = eval(batch_call, rho);
    SEXP gami = PROTECT(ScalarReal(gam[i]));
    t = s = PROTECT(allocList(5));
    SET_TYPEOF(s, LANGSXP);
    SETCAR(t, epoch); t = CDR(t);
    SETCAR(t, par); t = CDR(t);
    SETCAR(t, index); t = CDR(t);
    SETCAR(t, loss_gr); t = CDR(t);
    SETCAR(t, gami);
    par = eval(s, rho);
    UNPROTECT(2);
  }
  UNPROTECT(2);
  return par;
}