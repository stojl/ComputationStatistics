Rcpp::sourceCpp("gradients.cpp")
source("hess_grad.R")
source("Newton.R")
source("EM-algorithm.R")
Rcpp::sourceCpp("EM_rcpp.cpp")
set.seed(3939392)
nu1 <- 0.2
X <- extraDistr::rlst(100000, df = nu1, mu = 5, sigma = sqrt(2))

Newton_rcpp(X, c(median(X), IQR(X)), nu1, epsilon = 1e-7, maxit = 9L)

grad_func <- gradl(X, nu1)
func <- logl(X, nu1)
hess_func <- hessl(X, nu1)

Newton(c(median(X), IQR(X)), func, grad_func, hess_func, gamma0 = 1, d = 0.8, c = 0.2)


GradientTest(X, nu1, c(1,2))
grad_func(c(1,2))

LoglikeTest(X, nu1, c(1, 2))
func(c(1,2))

HessianTest(X, nu1, c(1, 2))
hess_func(c(1, 2))

EM_rcpp(X, nu1, c(median(X), IQR(X)))
EM(c(median(X), IQR(X)), x = X, nu = nu1)
                               
