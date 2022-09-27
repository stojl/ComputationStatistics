invFunc <- function(theta) {
  force(theta)
  
  c <- (exp(theta * 2) - exp(theta * (-1.9))) / theta
  
  function(x) {
    log(c * theta * x + exp(theta * (-1.9))) / theta
  }
  
}

ruin_prob_3 <- function(n, m, theta) {
  
  quantFunc <- invFunc(theta)
  
  samples <- quantFunc(runif(n * m))
  
  dim(samples) <- c(n, m)
  
  w <- numeric(m)
  
  for(i in seq_along(w)) {
    w[i] <- exp(-theta * sum(samples[, i]))
  }
  
  w <- w / sum(w)
  
  results <- logical(m)
  
  for(i in 1:m) {
    S <- 30 + cumsum(samples[, i])
    results[i] <- any(S <= 0)
  }
  
  sum(results * w)
}
