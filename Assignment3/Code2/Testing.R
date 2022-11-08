library(dplyr)
library(magrittr)
library(tidyr)
library(ggplot2)
library(extraDistr)
library(bench)
library(extraDistr)
library(numDeriv)
library(CSwR)
theme_set(theme_bw())
source("EM.R")
source("RCPP_EM.R")
source("Newton.R")
source("Newton_C.R")
source("Newton_RCPP.R")
source("GradientAndHessian.R")
source("RCPP_GradientAndHessian.R")
set.seed(12345)
nu <- 2.5
nu2 <- 0.5
X1 <- extraDistr::rlst(10000, nu, 5, sqrt(3))
X2 <- extraDistr::rlst(10000, nu2, 5, sqrt(3))

H_opt <- get_like(X1, nu)(optim(c(5, 3), get_like(X1, nu), method = "BFGS", control = list(reltol = 1e-30))$par)
H_opt2 <- get_like(X2, nu2)(optim(c(5, 3), get_like(X2, nu2), method = "BFGS", control = list(reltol = 1e-30))$par)

get_sum <- function(x, name, opt, x1, nu) {
  summary(x) %>% 
    rowwise() %>% 
    mutate(diff = get_like(x1, nu)(c(par1.1, par1.2)) - opt,
           alg = as.factor(name))
}
etracer1 <- tracer("par1", N = 0)
EM(x = X1, 
    nu = nu,
    maxit = 1000,
    min.eps = 1e-7,
    cb = etracer1$tracer)

ntracer1 <- tracer("par1", N = 0)
Newton(
  par = c(median(X1), IQR(X1)),
  H = get_like(X1, nu),
  gr = get_grad(X1, nu),
  hess = get_hess(X1, nu),
  maxit = 500,
  min.eps = 1e-7,
  cb = ntracer1$tracer
)

etracer2 <- tracer("par1", N = 0)
EM(x = X2, 
   nu = nu2,
   maxit = 1000,
   min.eps = 1e-7,
   cb = etracer2$tracer)

ntracer2 <- tracer("par1", N = 0)
Newton(
  par = c(median(X2), IQR(X2)),
  H = get_like(X2, nu2),
  gr = get_grad(X2, nu2),
  hess = get_hess(X2, nu2),
  maxit = 500,
  min.eps = 1e-7,
  cb = ntracer2$tracer
)

bind_rows(get_sum(ntracer1, "Newton", H_opt, X1, nu),
          get_sum(etracer1, "EM", H_opt, X1, nu),
          get_sum(ntracer2, "Newton0.5", H_opt2, X2, nu2),
          get_sum(etracer2, "EM0.5", H_opt2, X2, nu2)) %>% 
  ggplot(aes(x = .time, y = diff, col = alg)) +
  geom_line() +
  scale_y_log10()

bm0 <- mark(
  EM = EM(x = X1, 
          nu = nu,
          maxit = 500,
          min.eps = 1e-7),
  NEWTON = Newton(
    par = c(median(X1), IQR(X1)),
    H = get_like(X1, nu),
    gr = get_grad(X1, nu),
    hess = get_hess(X1, nu),
    maxit = 500,
    min.eps = 1e-7
  ),
  EM0.5 = EM(x = X2, 
          nu = nu2,
          maxit = 500,
          min.eps = 1e-7),
  NEWTON0.5 = Newton(
    par = c(median(X2), IQR(X2)),
    H = get_like(X2, nu2),
    gr = get_grad(X2, nu2),
    hess = get_hess(X2, nu2),
    maxit = 500,
    min.eps = 1e-7
  ),
  
  check = FALSE,
  min_iterations = 100
)

bm0

set.seed(12345)
nu <- 2.5
X1 <- extraDistr::rlst(100000, nu, 5, sqrt(3))

profvis::profvis(
  Newton(
    par = c(median(X1), IQR(X1)),
    H = get_like(X1, nu),
    gr = get_grad(X1, nu),
    hess = get_hess(X1, nu),
    maxit = 500,
    min.eps = 1e-7
  )
)

set.seed(12345)
nu <- 2.5
X1 <- extraDistr::rlst(10000, nu, 5, sqrt(3))


benchmark_cpp_grad_hess<- function() {
  test_like <- get_like(X1, nu)
  test_grad <- get_grad(X1, nu)
  test_hess <- get_hess(X1, nu)
  
  test_like_cpp <- get_like_cpp(X1, nu)
  test_grad_cpp <- get_grad_cpp(X1, nu)
  test_hess_cpp <- get_hess_cpp(X1, nu)
  
  mark(
    LIKE = test_like(c(1, 2)),
    LIKE_CPP = test_like_cpp(c(1, 2)),
    
    GRAD = test_grad(c(1, 2)),
    GRAD_CPP = test_grad_cpp(c(1, 2)),
    
    HESS = test_hess(c(1, 2)),
    HESS_CPP = test_hess_cpp(c(1, 2)),
    check = FALSE,
    min_iterations = 100
  )
}

bm1 <- benchmark_cpp_grad_hess()

ntracer3 <- tracer("par1", N = 0)
Newton(
  par = c(median(X1), IQR(X1)),
  H = get_like_cpp(X1, nu),
  gr = get_grad_cpp(X1, nu),
  hess = get_hess_cpp(X1, nu),
  maxit = 500,
  min.eps = 1e-7,
  cb = ntracer3$tracer
)

ntracer4 <- tracer("par1", N = 0)
Newton(
  par = c(median(X2), IQR(X2)),
  H = get_like_cpp(X2, nu2),
  gr = get_grad_cpp(X2, nu2),
  hess = get_hess_cpp(X2, nu2),
  maxit = 500,
  min.eps = 1e-7,
  cb = ntracer4$tracer
)

bind_rows(get_sum(ntracer1, "Newton", H_opt, X1, nu),
          get_sum(ntracer3, "Newton_CPP", H_opt, X1, nu),
          get_sum(etracer1, "EM", H_opt, X1, nu),
          get_sum(ntracer2, "Newton0.5", H_opt2, X2, nu2),
          get_sum(ntracer4, "Newton0.5_CPP", H_opt2, X2, nu2),
          get_sum(etracer2, "EM0.5", H_opt2, X2, nu2)) %>% 
  ggplot(aes(x = .time, y = diff, col = alg)) +
  geom_line() +
  scale_y_log10()


bm2 <- mark(
  EM = EM(x = X1, 
          nu = nu,
          maxit = 500,
          min.eps = 1e-7),
  NEWTON = Newton(
    par = c(median(X1), IQR(X1)),
    H = get_like(X1, nu),
    gr = get_grad(X1, nu),
    hess = get_hess(X1, nu),
    maxit = 500,
    min.eps = 1e-7
  ),
  NEWTON_CPP = Newton(
    par = c(median(X1), IQR(X1)),
    H = get_like_cpp(X1, nu),
    gr = get_grad_cpp(X1, nu),
    hess = get_hess_cpp(X1, nu),
    maxit = 500,
    min.eps = 1e-7
  ),
  EM0.5 = EM(x = X2, 
          nu = nu2,
          maxit = 500,
          min.eps = 1e-7),
  NEWTON0.5 = Newton(
    par = c(median(X2), IQR(X2)),
    H = get_like(X2, nu2),
    gr = get_grad(X2, nu2),
    hess = get_hess(X2, nu2),
    maxit = 500,
    min.eps = 1e-7
  ),
  NEWTON0.5_CPP = Newton(
    par = c(median(X2), IQR(X2)),
    H = get_like_cpp(X2, nu2),
    gr = get_grad_cpp(X2, nu2),
    hess = get_hess_cpp(X2, nu2),
    maxit = 500,
    min.eps = 1e-7
  ),
  check = FALSE,
  min_iterations = 100
)

bm2


bm3 <- mark(
  EM = EM(x = X1, 
          nu = nu,
          maxit = 500,
          min.eps = 1e-7),
  EM_CPP = EM_cpp(x = X1,
              nu = nu,
              maxit = 500,
              min.eps = 1e-7),
  NEWTON = Newton(
    par = c(median(X1), IQR(X1)),
    H = get_like(X1, nu),
    gr = get_grad(X1, nu),
    hess = get_hess(X1, nu),
    maxit = 500,
    min.eps = 1e-7
  ),
  NEWTON_CPP = Newton(
    par = c(median(X1), IQR(X1)),
    H = get_like_cpp(X1, nu),
    gr = get_grad_cpp(X1, nu),
    hess = get_hess_cpp(X1, nu),
    maxit = 500,
    min.eps = 1e-7
  ),
  NEWTON_C = Newton_C(
    par = c(median(X1), IQR(X1)),
    H = get_like(X1, nu),
    gr = get_grad(X1, nu),
    hess = get_hess(X1, nu),
    maxit = 500,
    min.eps = 1e-7),
  NEWTON_C_CPP = Newton_C(
    par = c(median(X1), IQR(X1)),
    H = get_like_cpp(X1, nu),
    gr = get_grad_cpp(X1, nu),
    hess = get_hess_cpp(X1, nu),
    maxit = 500,
    min.eps = 1e-7),
  EM0.5 = EM(x = X2, 
             nu = nu2,
             maxit = 500,
             min.eps = 1e-7),
  EM0.5_CPP = EM_cpp(x = X2,
                  nu = nu2,
                  maxit = 500,
                  min.eps = 1e-7),
  NEWTON0.5 = Newton(
    par = c(median(X2), IQR(X2)),
    H = get_like(X2, nu2),
    gr = get_grad(X2, nu2),
    hess = get_hess(X2, nu2),
    maxit = 500,
    min.eps = 1e-7
  ),
  NEWTON0.5_CPP = Newton(
    par = c(median(X2), IQR(X2)),
    H = get_like_cpp(X2, nu2),
    gr = get_grad_cpp(X2, nu2),
    hess = get_hess_cpp(X2, nu2),
    maxit = 500,
    min.eps = 1e-7
  ),
  NEWTON0.5_C = Newton_C(
    par = c(median(X2), IQR(X2)),
    H = get_like(X2, nu2),
    gr = get_grad(X2, nu2),
    hess = get_hess(X2, nu2),
    maxit = 500,
    min.eps = 1e-7),
  NEWTON0.5_C_CPP = Newton_C(
    par = c(median(X2), IQR(X2)),
    H = get_like_cpp(X2, nu2),
    gr = get_grad_cpp(X2, nu2),
    hess = get_hess_cpp(X2, nu2),
    maxit = 500,
    min.eps = 1e-7),
  check = FALSE,
  min_iterations = 100
)
bm3

bm4 <- mark(
  EM = EM(x = X1, 
          nu = nu,
          maxit = 500,
          min.eps = 1e-7),
  EM_CPP = EM_cpp(x = X1,
                  nu = nu,
                  maxit = 500,
                  min.eps = 1e-7),
  NEWTON = Newton(
    par = c(median(X1), IQR(X1)),
    H = get_like(X1, nu),
    gr = get_grad(X1, nu),
    hess = get_hess(X1, nu),
    maxit = 500,
    min.eps = 1e-7
  ),
  NEWTON_CPP = Newton(
    par = c(median(X1), IQR(X1)),
    H = get_like_cpp(X1, nu),
    gr = get_grad_cpp(X1, nu),
    hess = get_hess_cpp(X1, nu),
    maxit = 500,
    min.eps = 1e-7
  ),
  NEWTON_C = Newton_C(
    par = c(median(X1), IQR(X1)),
    H = get_like(X1, nu),
    gr = get_grad(X1, nu),
    hess = get_hess(X1, nu),
    maxit = 500,
    min.eps = 1e-7),
  NEWTON_C_CPP = Newton_C(
    par = c(median(X1), IQR(X1)),
    H = get_like_cpp(X1, nu),
    gr = get_grad_cpp(X1, nu),
    hess = get_hess_cpp(X1, nu),
    maxit = 500,
    min.eps = 1e-7),
  NEWTON_ARMA = Newton_cpp(
    par = c(median(X1), IQR(X1)),
    x = X1,
    nu = nu,
    maxit = 500L,
    min.eps = 1e-9
  ),
  EM0.5 = EM(x = X2, 
             nu = nu2,
             maxit = 500,
             min.eps = 1e-7),
  EM0.5_CPP = EM_cpp(x = X2,
                     nu = nu2,
                     maxit = 500,
                     min.eps = 1e-7),
  NEWTON0.5 = Newton(
    par = c(median(X2), IQR(X2)),
    H = get_like(X2, nu2),
    gr = get_grad(X2, nu2),
    hess = get_hess(X2, nu2),
    maxit = 500,
    min.eps = 1e-7
  ),
  NEWTON0.5_CPP = Newton(
    par = c(median(X2), IQR(X2)),
    H = get_like_cpp(X2, nu2),
    gr = get_grad_cpp(X2, nu2),
    hess = get_hess_cpp(X2, nu2),
    maxit = 500,
    min.eps = 1e-7
  ),
  NEWTON0.5_C = Newton_C(
    par = c(median(X2), IQR(X2)),
    H = get_like(X2, nu2),
    gr = get_grad(X2, nu2),
    hess = get_hess(X2, nu2),
    maxit = 500,
    min.eps = 1e-7),
  NEWTON0.5_C_CPP = Newton_C(
    par = c(median(X2), IQR(X2)),
    H = get_like_cpp(X2, nu2),
    gr = get_grad_cpp(X2, nu2),
    hess = get_hess_cpp(X2, nu2),
    maxit = 500,
    min.eps = 1e-7),
  NEWTON0.5_ARMA = Newton_cpp(
    par = c(median(X2), IQR(X2)),
    x = X2,
    nu = nu2,
    maxit = 500L,
    min.eps = 1e-1
  ),
  check = FALSE,
  min_iterations = 100
)
bm4
