ruin_prob_1 <- function(n, m) {

  samples <- runif(n * m, min = -1.9, max = 2)
  
  dim(samples) <- c(n, m)
  
  results <- logical(m)
  
  for(i in 1:m) {
    S <- 30 + cumsum(samples[, i])
    results[i] <- any(S <= 0)
  }
  
  mean(results)
}
