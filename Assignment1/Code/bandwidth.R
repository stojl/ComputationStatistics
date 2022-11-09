if(!is.loaded("src/dens.o")) {
  dyn.load("src/dens.o")
}
Rcpp::sourceCpp("dens_rcpp.cpp")
bw_oracle <- function(x, kernel) {
  n <- length(x)
  K <- integrate(function(x) kernel(x)^2, -Inf, Inf)$value
  sigma2 <- integrate(function(x) kernel(x) * x^2, -Inf, Inf)$value
  sigma <- min(sd(x), IQR(x) / 1.34)
  (8 * sqrt(pi) * K / (3 * sigma2^2))^(1/5) * sigma * n^(-1/5)
}

bw_cv <- function(x, kernel, max_bw = 2) {
  cv_func <- function(l) .Call("C_cv", x, kernel, l, environment())
  optimize(cv_func, c(.Machine$double.eps, max_bw))$minimum
}

bw_cv_cpp_partial <- function(x, kernel, max_bw = 2) {
  cv_func <- function(l) bw_cv_rcpp_partial(x, kernel, l)
  optimize(cv_func, c(.Machine$double.eps, max_bw))$minimum
}

bw_cv_cpp <- function(x, kernel, max_bw = 2) {
  cv_func <- function(l) bw_cv_rcpp(x, l)
  optimize(cv_func, c(.Machine$double.eps, max_bw))$minimum
}

bw_cv_R <- function(x, kernel, max_bw = 2) {
  cv_func <- function(l) {
    n <- length(x)
    K <- numeric(n)
    for(i in 2:n) {
      index <- 1:(i - 1)
      tmp <- kernel((x[i] - x[index]) / l)
      K[i] <- K[i] + sum(tmp)
      K[index] <- K[index] + tmp[index]
    }
    cv <- sum(log(K[K > 0]))
    n * log((n - 1) * l) - cv
  }
  optimize(cv_func, c(.Machine$double.eps, max_bw))$minimum
}

bw_cv_R2 <- function(x, kernel, max_bw = 2) {
  cv_func <- function(l) {
    n <- length(x)
    K <- numeric(n)
    for(i in 1:n) {
      K[i] <- sum(kernel((x[i] - x[-i]) / l))
    }
    cv <- sum(log(K[K > 0]))
    n * log((n - 1) * l) - cv
  }
  optimize(cv_func, c(.Machine$double.eps, max_bw))$minimum
}
