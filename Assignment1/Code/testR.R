cv_func <- function(x, kernel, l) {
  if(l < .Machine$double.eps) Inf
  n <- length(x)
  K <- numeric(n)
  for(i in 2:n) {
    index <- 1:(i - 1)
    tmp <- kernel(x[i] - x[index])
    K[i] <- K[i] + sum(tmp)
    K[index] <- K[index] + tmp[index]
  }
  cv <- sum(log(K[K > .Machine$double.eps]))
  n * log((n - 1) * l) - cv
}