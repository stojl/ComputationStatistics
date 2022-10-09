E_step <- function(x, nu) {
  force(x)
  force(nu)
  
  function(p_theta) {
    mu <- mean(x)
    w <- (nu + 1) / (1 + ((x - p_theta[1])^2) / (nu * p_theta[2]))
    mu <- sum(w * x) / sum(w)
    sigma <- mean(w * ((x - mu)^2)) / nu
    
    c(mu, sigma)
  }
}

simulate_X <- function(N, mu, sigma, nu) {
  W <- rchisq(N, nu)
  X <- rnorm(N, mu, sqrt(nu * sigma / W))
  
  list(x = X, w = W)
}



estimate <- function(x, theta, N) {
  E <- E_step(x, 2)
  for(i in 1:N) {
    theta <- E(theta)
  }
  theta
}



full_mle <- function(X, W, nu, df) {
  stopifnot(length(X) == length(W))
  mu <- sum(X * W) / sum(W)
  sigma <- sum(W * (X - mu)^2) / (nu * (length(X) - df))
  list(mu = mu, sigma = sigma)
}


full_mle(X$x, X$w, 3)

X <- simulate_X(100000, 5, 1.5, 3)
X <- extraDistr::rlst(100000, 6, 5, sqrt(1.5))
result <- estimate(X, c(3,6), 1000)
result

estimate_df(X, result[2])

estimate_df <- function(x, sigma) {
  v <- var(x)
  2 / (1 - sigma / v)
}


