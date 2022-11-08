poisdata <- read.csv("Poisson.csv", header = TRUE)
Rcpp::sourceCpp("RCPP_SAMP.cpp")
poisreg <- function(x, z) {
  force(x)
  force(z)
  function(y) {
    expyx <- sapply(y, function(s) sum(exp(s * x)))
    exp(y * sum(x * z) - expyx)
  }
}

poisreg_cpp <- function(x, z) {
  force(x)
  force(z)
  function(y) {
    RCPP_poisdens(y, x, z)
  }
}

poisreg_derv <- function(x, z) {
  force(x)
  force(z)
  function(y) {
    expyx <- sapply(y, function(s) sum(exp(s * x)))
    x_expyx <- sapply(y, function(s) sum(x * exp(s * x)))
    xz <- sum(x * z) 
    exp(y * xz - expyx) * (xz - x_expyx)
  }
}

poisreg_derv_cpp <- function(x, z) {
  force(x)
  force(z)
  function(y) {
    RCPP_poisdens_derv(y, x, z)
  }
}

gausdens <- function(mu, sigma) {
  force(mu)
  force(sigma)
  function(x) {
    exp(-(x - mu)^2 / (2 * sigma))
  }
}

poisdens <- poisreg(poisdata$x, poisdata$z)
poisdens_derv <- poisreg_derv(poisdata$x, poisdata$z)
poisdens_cpp <- poisreg(poisdata$x, poisdata$z)
poisdens_derv_cpp <- poisreg_derv(poisdata$x, poisdata$z)
