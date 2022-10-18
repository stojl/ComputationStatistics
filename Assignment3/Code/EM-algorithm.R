E_step <- function(x, nu) {
  force(x)
  force(nu)
  function(par) {
    # mu = par[1]
    # sigma^2 = par[2]
    test <- (nu + 1) / (1 + ((x - par[1])^2) / (nu * par[2]))
    test
  }
}

M_step <- function(x, nu) {
  force(x)
  force(nu)
  function(EW) {
    mu <- sum(EW * x) / sum(EW)
    sigma <- mean(EW * (x - mu)^2) / nu
    c(mu, sigma)
  }
}

EM <- function(par, x, nu, maxit = 500, min.eps = 1e-7) {
  E <- E_step(x, nu)
  M <- M_step(x, nu)
  for(i in 1:maxit) {
    EW <- E(par)
    new_par <- M(EW)
    if(sum((new_par - par)^2) < min.eps * (sum(par^2) + min.eps)) {
      par <- new_par
      break
    }
    par <- new_par
    if(i == maxit) warning("Maximum number of itertaions reached.")
  }
  names(par) <- c("mu", "sigma")
  list(par = c(par, nu = nu), iterations = i)
}

EM_2 <- function(par, x, nu, maxit = 500, min.eps = 1e-7) {
  new_par <- par
  for(i in 1:maxit) {
    w <- (nu + 1) / (1 + ((x - par[1])^2) / (nu * par[2]))
    new_par[1] <- sum(w * x) / sum(w)
    new_par[2] <- mean(w * (x - par[1])^2) / nu
    if(sum((new_par - par)^2) < min.eps * (sum(par^2) + min.eps)) {
      par <- new_par
      break
    }
    par <- new_par
    if(i == maxit) warning("Maximum number of itertaions reached.")
  }
  par
}

X <- extraDistr::rlst(100000, df = 6, mu = 5, sigma = sqrt(2))

set.seed(3939392)
simulate_X <- function(N, mu, sigma, nu) {
  W <- rchisq(N, nu)
  X <- rnorm(N, mu, sqrt(nu * sigma / W))
  
  list(x = X, w = W)
}

full_mle <- function(X, W, nu) {
  stopifnot(length(X) == length(W))
  mu <- sum(X * W) / sum(W)
  sigma <- sum(W * (X - mu)^2) / (nu * (length(X)))
  list(mu = mu, sigma = sigma)
}

samples <- simulate_X(100000, 5, 1.5, 3)
full_mle(samples$x, samples$w, 3)


EM(c(1,2), samples$x, 3)



