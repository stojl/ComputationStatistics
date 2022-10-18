logl <- function(x, nu) {
  n <- length(x)
  force(nu)
  
  function(par) {
    mu <- par[1]
    sigma <- par[2]
    
    K <- sum(log(1 + (x - mu)^2 / (nu * sigma)))
    
    log(sigma) / 2 + (nu + 1) * K / (2 * n)
  }
}

gradl <- function(x, nu) {
  n <- length(x)
  force(nu)
  
  function(par) {
    mu <- par[1]
    sigma <- par[2]
    
    C1 <- (x - mu) / (1 + (x - mu)^2 / (nu * sigma))
    
    K_mu <- sum(C1)
    K_sigma <- sum(C1 * (x - mu))
    
    grad_mu <- -(nu + 1) * K_mu / (n * nu * sigma)
    grad_sigma <- 1 / (2 * sigma) - (nu + 1) * K_sigma / (2 * n * nu * sigma^2)
    
    c(grad_mu, grad_sigma)
  }
}

hessl <- function(x, nu) {
  n <- length(x)
  force(nu)
  
  function(par){
    mu <- par[1]
    sigma <- par[2]
    
    C0 <- 1 / (1 + (x - mu)^2 / (nu * sigma))
    C1 <- C0 * (x - mu)
    C2 <- C1 * (x - mu)
    
    hess_mu <- -(nu + 1) * sum(C0) / (2 * n * nu * sigma) + 
      (nu + 1) * sum(C1^2) / (n * (nu * sigma)^2)
    
    hess_sigma <- -1 / (2 * sigma^2) + 
      (nu + 1) * sum(C2) / (n * nu * sigma^3) -
      (nu + 1) * sum(C2^2) / (2 * n * nu^2 * sigma^4)
    
    hess_mu_sigma <- (nu + 1) * sum(C1) / (n * nu * sigma^2) -
      (nu + 1) * sum(C1^2 * (x - mu)) / (n * nu * sigma^4)
    
    hess <- c(hess_mu, hess_mu_sigma, hess_mu_sigma, hess_sigma)
    dim(hess) <- c(2, 2)
    hess
  }
}

nu1 <- 1.5
X <- extraDistr::rlst(100000, df = nu1, mu = 30, sigma = sqrt(2))

loglike <- logl(X, nu1)
loggrad <- gradl(X, nu1)

marginal_MLE_nu <- function(x, nu, ...) {
  loglike <- logl(x, nu)
  loggrad <- gradl(x, nu)
  
  mu0 <- mean(x)
  sigma0 <- var(x) * nu / (nu + 2)
  
  optim(c(mu0, sigma0), loglike, loggrad, ...)
}

marginal_MLE_nu(X, nu1, method = "BFGS")
EM_1(c(mean(X), var(X) * nu1 / (nu1 + 2)), X, nu1, maxit = 1000)
