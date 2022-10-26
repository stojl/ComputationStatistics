E_step <- function(x, nu) {
  function(par) {
    mu <- par[1]
    sigma <- par[2]
    test <- (nu + 1) / (1 + ((x - mu)^2) / (nu * sigma))
    test
  }
}

M_step <- function(x, nu) {
  function(EW) {
    mu <- sum(EW * x) / sum(EW)
    sigma <- mean(EW * (x - mu)^2) / nu
    c(mu, sigma)
  }
}

EM <- function(par = NULL, x, nu, cb = NULL, maxit = 500, min.eps = 1e-7) {
  E <- E_step(x, nu)
  M <- M_step(x, nu) 
  if(is.null(par)) {
    par <- c(median(x), IQR(x))
  }
  for(i in 1:maxit) {
    EW <- E(par)
    par1 <- M(EW)
    if(!is.null(cb)) cb()
    if(norm(par - par1, "2") < min.eps * (norm(par, "2") + min.eps)) break
    par <- par1
  }
  if(i == maxit) warning("Maximum number of itertaions ", maxit, " reached.")
  names(par1) <- c("mu", "sigma")
  list(par1 = c(par, nu = nu), iterations = i)
}

EM2 <- function(par = NULL, x, nu, cb = NULL, maxit = 500, min.eps = 1e-7) {
  if(is.null(par)) {
    par <- c(median(x), IQR(x))
  }
  par1 <- numeric(2)
  n <- length(x)
  EW <- numeric(n)
  for(i in 1:maxit) {
    EW <- (nu + 1) / (1 + ((x - par[1])^2) / (nu * par[2]))
    par1[1] <- sum(EW * x) / sum(EW)
    par1[2] <- sum(EW * (x - par[1])^2) / (n * nu)
    if(!is.null(cb)) cb()
    if(norm(par - par1, "2") < min.eps * (norm(par, "2") + min.eps)) break
    par <- par1
  }
  if(i == maxit) warning("Maximum number of itertaions ", maxit, " reached.")
  names(par1) <- c("mu", "sigma")
  list(par = c(par1, nu = nu), iterations = i)
}
