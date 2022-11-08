#include <R.h>
#include <math.h>
#include <Rinternals.h>

extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                   double* b, int* ldb, int* info );

extern void dgesvd_(char *jobu, char *jobvt, int *m, int *n, double *a,
                    int *lda, double *s, double *u, int *ldu, double *vt,
                    int *ldvt, double *work, int *lwork, int *info);

double norm2(double *x, int m) {
  char jobu = 'N', jobvt = 'N';
  int n = 1, info, lwork;
  if(m > 3) {
    lwork = m + 2;
  } else {
    lwork = 5;
  }
  double *U, *VT, s;
  double *A = (double*)malloc(sizeof(double) * m);
  double *work = (double*)malloc(sizeof(double) * lwork);
  for(int i = 0; i < m; ++i)
    A[i] = x[i];
  dgesvd_(&jobu, &jobvt, &m, &n, A, &m, &s, 
          U, &m, VT, &m, work, &lwork, &info);
  free(A);
  free(work);
  return s;
}

void solve(int n, double *a, double *b, double *result, int *info) {
  for(int i = 0; i < n; ++i)
    result[i] = b[i];
  int nrhs = 1;
  int *ipiv = (int*)malloc(sizeof(int) * n);
  dgesv_(&n, &nrhs, a, &n, ipiv, result, &n, info);
  free(ipiv);
}

SEXP C_Newton(SEXP par0, SEXP H, SEXP gr, SEXP hess,
              SEXP d, SEXP c, SEXP gamma0, 
              SEXP eps, SEXP maxit,
              SEXP env) {
  int n = length(par0), info;
  SEXP H_call = PROTECT(lang2(H, R_NilValue));
  SEXP gr_call = PROTECT(lang2(gr, R_NilValue));
  SEXP hess_call = PROTECT(lang2(hess, R_NilValue));
  SEXP value, grad, hessian, Hpar;
  SEXP par = PROTECT(allocVector(REALSXP, n));
  SEXP par1 = PROTECT(allocVector(REALSXP, n));
  double value_, *grad_, *hessian_, *par1_ = REAL(par1), *par_ = REAL(par);
  double gamma0_ = REAL(gamma0)[0], d_ = REAL(d)[0], c_ = REAL(c)[0];
  double eps_ = REAL(eps)[0], maxit_ = INTEGER(maxit)[0];
  double *rho = (double*)malloc(sizeof(double) * n);
  
  for(int i = 0; i < n; ++i)
    REAL(par)[i] = REAL(par0)[i];
  int k;
  
  for(k = 0; k < maxit_; ++k) {
    double gamma = gamma0_;
    double h_prime;
    SETCADR(H_call, par);
    value = PROTECT(eval(H_call, env));
    SETCADR(gr_call, par);
    grad = PROTECT(eval(gr_call, env));
    SETCADR(hess_call, par);
    hessian = PROTECT(eval(hess_call, env));
    hessian_ = REAL(hessian);
    grad_ = REAL(grad);
    for(int j = 0; j < n; ++j)
      grad_[j] = -grad_[j];
    value_ = REAL(value)[0];
    solve(n, hessian_, grad_, rho, &info);
    if(info != 0) break;
    for(int j = 0; j < n; ++j) {
      par1_[j] = par_[j] + gamma * rho[j];
      h_prime += grad_[j] * rho [j];
    }
    while(1 == 1) {
      SETCADR(H_call, par1);
      Hpar = eval(H_call, env);
      if(ISNAN(REAL(Hpar)[0]) || REAL(Hpar)[0] > value_ + c_ * gamma * h_prime) {
        gamma = d_ * gamma;
        for(int j = 0; j < n; ++j)
          par1_[j] = par_[j] + gamma * rho[j];
      } else break;
    }
    UNPROTECT(3);
    double norm = norm2(par1_, n);
    for(int j = 0; j < n; ++j)
      par_[j] -= par1_[j];
    double norm_new = norm2(par_, n);
    if(norm_new < eps_ * (norm + eps_)) break;
    for(int j = 0; j < n; ++j) {
      par_[j] = par1_[j];
    }
  }
  SEXP result = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result, 0, par1);
  if(k != maxit_) k += 1;
  if(info != 0) {
    SET_VECTOR_ELT(result, 1, ScalarInteger(-1));
  } else {
    SET_VECTOR_ELT(result, 1, ScalarInteger(k));
  }
  UNPROTECT(6);
  free(rho);
  return result;
}