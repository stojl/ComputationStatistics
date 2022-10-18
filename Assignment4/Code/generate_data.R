create_regfunc <- function(alpha, beta, gamma, rho) {
  force(alpha)
  force(beta)
  force(gamma)
  force(rho)
  
  function(x) {
    gamma + (rho - gamma) / (1 + exp(beta * log(x) - alpha))
  }
}

generate_x <- function(n, sd) {
  exp(rnorm(n, mean = 0, sd = sd))
}

generate_y <- function(x, sigma, alpha = 0, beta = 1, gamma = 1, rho = 1) {
  regfunc <- create_regfunc(alpha, beta, gamma, rho)
  regfunc(x) + rnorm(length(x), mean = 0, sd = sigma)
}

X <- generate_x(1000, 2)
Y <- generate_y(X, 1.5, alpha = 1, beta = 3, gamma = 5, rho = 0.2)
