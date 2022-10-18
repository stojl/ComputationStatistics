funcl <- function(x) {
  force(x)
  n <- length(x)
  
  function(par) {
    mu <- par[1]
    sigma <- par[2]
    nu <- par[3]
    
    K1 <- sum(log(1 + (x - mu)^2 / (nu * sigma)))
    
    -log(gamma((nu + 1) / 2)) +
      log(nu) / 2 +
      log(sigma) / 2 +
      log(gamma(nu / 2)) +
      (nu + 1) * K1 / (2 * n)
  }
}

gradl <- function(x) {
  force(x)
  n <- length(x)
  function(par) {
    mu <- par[1]
    sigma <- par[2]
    nu <- par[3]
    
    C1 <- (x - mu) / (1 + (x - mu)^2 / (nu * sigma))
    
    K1 <- sum(C1)
    K2 <- sum(C1 * (x - mu))
    K3 <- sum(log(C1 * (x - mu)))
    
    mu_out <- -K1 * (nu + 1) / (n * nu * sigma)
    sigma_out <- 1 / (2 * sigma) - (nu + 1) * K2 / (2 * n * nu * sigma^2)
    nu_out <- -digamma((nu + 1) / 2) + 
      1 / (2 * nu) +
      digamma(nu) / 2 +
      K3 / (2 * n) -
      (nu + 1) * K2 / (2 * nu^2 * sigma * n)
    
    c(mu_out, sigma_out, nu_out)
  }
  
}

estimate_marginal_likelihood <- function(x) {
  like <- funcl(x)
  
  out <- suppressWarnings(optim(c(mean(x), var(x), 1), like)$par)
  
  names(out) <- c("mu", "sigma", "nu")
  
  out
}

X <- extraDistr::rlst(100000, df = 1.5, mu = 30, sigma = sqrt(2))

estimate_marginal_likelihood(X)

densl <- function(mu, sigma, nu) {
  force(mu)
  force(sigma)
  force(nu)
  function(x) {
    gamma((nu + 1) / 2) * (1 + (x - mu)^2 / (nu * sigma))^(-(nu + 1) / 2) / 
      (gamma(nu / 2) * sqrt(pi * nu * sigma))
  }
}

dens1 <- densl(30, 2, 100)
dens2 <- densl(30, 2, 30)

plot(seq(20, 40, 0.05), dens1(seq(20, 40, 0.05)))
points(seq(20, 40, 0.05), dens2(seq(20, 40, 0.05)))
points(seq(20, 40, 0.05), dnorm(seq(20, 40, 0.05), 30, sqrt(2)))


kullback_div <- function(mu, sigma, nu) {
  force(mu)
  force(sigma)
  force(nu)
  
  function(nup) {
    log(gamma((nup + 1) / 2)) - log(gamma((nu + 1) / 2)) -
      log(nup) / 2 + log(nu) / 2 -
      log(gamma(nup / 2)) + log(gamma(nu / 2)) -
      (nup + 1) / 2 + (nu + 1) / 2 +
      (nu + 1) * nup / (2 * nu * (nup - 2)) -
      (nup + 1) / (2 * (nup - 2))
  }
}

test <- kullback_div(30, 2, 100)
plot(seq(20, 200, 1), test(seq(20, 200, 1)))
