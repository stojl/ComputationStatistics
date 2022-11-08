R_dens <- function(x, p, kernel, bandwidth) {
  m <- length(p)
  n <- length(x)
  result <- numeric(m)
  for(i in 1:m) {
    result[i] <- sum(kernel((p[i] - x) / bandwidth))
  }
  result / (n * bandwidth)
}

R_dens_for <- function(x, p, kernel, bandwidth) {
  m <- length(p)
  n <- length(x)
  result <- numeric(m)
  for(i in 1:m) {
    for(j in 1:n) {
      result[i] <- result[i] + kernel((p[i] - x[j]) / bandwidth)
    }
  }
  result / (n * bandwidth)
}

R_dens_apply <- function(x, p, kernel, bandwidth) {
  app_func <- function(t) sum(kernel((t - x) / bandwidth))
  result <- sapply(p, app_func)
  result / (length(x) * bandwidth)
}