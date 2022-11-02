#include <Rinternals.h>
#include <R.h>
#include <stdlib.h>
#include <R_ext/Rdynload.h>

extern SEXP trisum_CV(SEXP x, SEXP fn, SEXP rho, SEXP lambda)

static const R_CallMethodDef CallEntries[] = {
  {"trisum_cv", (DL_FUNC) &trisum_CV, 4},
  {NULL, NULL, NULL, 0}
};

void R_init_Density(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}