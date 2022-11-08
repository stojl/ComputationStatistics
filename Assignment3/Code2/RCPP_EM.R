Rcpp::sourceCpp("src/CPP_EM.cpp")
EM_cpp <- function(par = NULL, x, nu, cb = NULL, maxit = 500, min.eps = 1e-7) {
  if(is.null(par)) par <- c(median(x), IQR(x))
  par1 <- numeric(2)
  n <- length(x)
  EW <- numeric(n)
  result <- CPP_EM(par, x, nu, maxit, min.eps)
  if(result[[2]] == maxit) warning("Maximum number of itertaions ", maxit, " reached.")
  names(result[[1]]) <- c("mu", "sigma_sq")
  structure(
    list(par = c(result[[1]]), 
         iterations = result[[2]], 
         nu = nu, 
         x = x),
    class = "em_estimate"
    )
}

EM_fisher <- function(mle, x, nu) {
  Q <- Q_func(x, nu, mle)
  phi <- phi_func(x, nu)
  IY <- numDeriv::hessian(Q, mle)
  IX <- (diag(1, 2) - t(numDeriv::jacobian(phi, mle))) %*% IY
  solve(IX)
}

Q_func <- function(x, nu, mle) {
  force(x); force(nu); force(par);
  function(par) Q_cpp(par, mle, x, nu)
}

phi_func <- function(x, nu) {
  force(x); force(nu)
  function(par) phi_cpp(par, x, nu)
}

confint.em_estimate <- function(x, level = 0.95) {
  qq <- level + (1 - level) / 2
  qq <- qnorm(qq)
  invf <- EM_fisher(x$par, x$x, x$nu)
  mu <- x$par[1] + c(-1, 1) * qq * sqrt(invf[1, 1])
  sigma <- x$par[2] + c(-1, 1) * qq * sqrt(invf[2, 2])
  names(mu) <- c("lwr", "upr")
  names(sigma) <- c("lwr", "upr")
  list(mu = mu, sigma_sq = sigma)
}

print.em_estimate <- function(x) {
  print(x[1:3])
}
