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
EM_tracer <- CSwR::tracer(c("new_par"))
EM(x = X, nu = nu1, cb = EM_tracer$tracer)

EM_info <- summary(EM_tracer)
EM_info
