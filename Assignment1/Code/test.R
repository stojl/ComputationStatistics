dyn.load("src/test.o")

kernel <- function(x) {
  if(abs(x) > 1) {
    0
  } else {
    0.75 * (1 - x^2)
  }
}

x <- c(rnorm(512, -1, 0.1), rnorm(512, 1, 0.1))
.Call("trisum_cv", x, kernel, environment(), 0.2)
objfunc <- function(lambda) {
  .Call("trisum_cv", x, kernel, environment(), lambda)
}

plot(seq(0.01, 0.1, 0.01), sapply(seq(0.01, 0.1, 0.01), objfunc))
plot(seq(0.1, 1, 0.05), sapply(seq(0.1, 1, 0.05), function(lambda) trisum_cv(x, lambda)))

CV <- function(x, kernel) {
  objfunc <- function(lambda) {
    .Call("trisum_cv", x, kernel, environment(), lambda)
  }
  
  optimise(objfunc, c(0.01, 1))$minimum
}

plot(x, trisum(x, CV(x, kernel)))

trisum_cv_R <- function(x, kernel, lambda) {
  K <- numeric(length(x))
  n <- length(x)
  
  for(i in 2:n) {
    for(j in 1:(i - 1)) {
      tmp <- kernel((x[i] - x[j]) / lambda)
      K[i] <- K[i] + tmp
      K[j] <- K[j] + tmp
    }
  }
  
  n * log((n - 1) * lambda) - sum(log(K))
}

microbenchmark::microbenchmark(
  R = trisum_cv_R(x, kernel, 0.5),
  RAPI = .Call("trisum_cv", x, kernel, environment(), 0.5)
)
