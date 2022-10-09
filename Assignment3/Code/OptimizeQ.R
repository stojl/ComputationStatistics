like_Q <- function(x, mu_p, sigma_p, nu_p) {
  n <- length(x)
  C1 <- n * digamma((nu_p + 1) / nu_p) + 
    n * log(2) -
    sum(log(1 + (x - mu_p)^2 / (nu_p * sigma_p)))
  
  function(nu) {
    -C1 / nu -
      n * log(2) * (nu + 1) / 2 -
      n * log(gamma(nu / 2))
  }
}

grad_Q <- function(x, mu_p, sigma_p, nu_p) {
  n <- length(x)
  
  C1 <- n * digamma((nu_p + 1) / nu_p) + 
    n * log(2) -
    sum(log(1 + (x - mu_p)^2 / (nu_p * sigma_p)))
  
  function(nu) {
    out <- C1 / nu^2 -
      n * log(2) / 2 - 
      n * digamma(nu / 2) / 2
    
    out / n
  }
}

f <-like_Q(X$x, 5, 1.95, 3)
gradf <- grad_Q(X$x, 4.95, 1.95, 2.22)
gradg <- function(x) -gradf(x) / length(X$x)
g <- function(x) -f(x) / length(X$x)

newton <- function(x0, gr, gam, a, N) {
  for(i in 1:N) {
    x0 <- x0 + gam * gr(x0)
  }
  x0
}

newton(2, gradg, 0.1, 1, 50)

estimate_full <- function(par, x, N, k) {
  browser()
  for(i in 1:N) {
    w <- 1 + (x - par[1])^2 / (par[3] * par[2])
    mu <- sum(w * x) / sum(w)
    gradf <- grad_Q2(x, par[1], par[2], par[3])
    nu <- newton(par[3], gradf, 0.1, 1, k)
    sigma <- mean(w * (x - mu)^2) * gamma(1 + 1 / par[3]) * 2^(-1 - 1 / par[3]) / nu
    par[1] <- mu
    par[2] <- sigma
    par[3] <- nu
  }
  par
}

grad_Q2 <- function(x, mu_p, sigma_p, nu_p) {
  n <- length(x)
  
  C1 <- digamma(1 + nu_p) - log(1 + (x - mu_p)^2 / (nu_p * sigma_p))
  
  K1 <- sum(C1)
  
  function(nu) {
    out <- digamma(nu / 2) / 2 - K1 / (2 * n)
    out
  }
}

grad_Q3 <- function(x, mu_p, sigma_p, nu_p) {
  n <- length(x)
  
  C1 <- digamma(1 + nu_p) - log(1 + (x - mu_p)^2 / (nu_p * sigma_p))
  C2 <- (nu_p + 1) / (1 + (x - mu_p)^2 / (nu_p * sigma_p))
  
  K1 <- sum(C1)
  
  
  function(par) {
    K2 <- sum(C2 * (x - par[1])^2)
    
    mu_out <- K2 / (par[2]* par[3])
    
    nu_out <- 1 / par[3] + 
      digamma(par[3] / 2) -
      K1 / n +
      K2 / (n * par[3]^2 * par[2])
    
    sigma_out <- 1 / par[1] -
      K2 / (n * par[3] * par[2]^2)
    
    c(mu_out, sigma_out, nu_out)
  }
  
}

like_Q3 <- function(x, mu_p, sigma_p, nu_p) {
  n <- length(x)
  
  C1 <- digamma(1 + nu_p) - log(1 + (x - mu_p)^2 / (nu_p * sigma_p))
  C2 <- (nu_p + 1) / (1 + (x - mu_p)^2 / (nu_p * sigma_p))
  
  K1 <- sum(C1)
  
  function(par) {
    K2 <- sum(C2 * (x - par[1])^2)
    
    log(par[3]) +
      log(par[2]) +
      log(gamma(par[3] / 2)) / 2 -
      K1 * par[3] / n +
      K2 / (n * par[3] * par[2])
  }
}

EM <- function(par, x, N) {
  for(i in 1:N) {
browser()
    Q <- like_Q3(x, par[1], par[2], par[3])
    Qgr <- grad_Q3(x, par[1], par[2], par[3])
    
    sigma_nu <- optim(c(par[1], par[2], par[3]), Q, Qgr, method = "BFGS")
    
    par[1] <- sigma_nu$par[1]
    par[2] <- sigma_nu$par[2]
    par[3] <- sigma_nu$par[3]
  }
  
  par
}
