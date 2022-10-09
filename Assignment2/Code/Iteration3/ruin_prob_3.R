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
  
  w_star <- numeric(m)
  
  for(i in seq_along(w_star)) {
    w_star[i] <- exp(-theta * sum(samples[, i]))
  }
  
  w <- w_star / sum(w_star)
  
  
  results <- logical(m)
  
  for(i in 1:m) {
    S <- 30 + cumsum(samples[, i])
    results[i] <- any(S <= 0)
  }
  
  sigma_IS <- var(results * w_star)
  gamma <- cov(results * w_star, w_star)
  sigma_w <- var(w_star)
  c <- mean(w_star)
  mu <- sum(results * w)
  
  sd <- sqrt((sigma_IS + mu^2 * sigma_w - 2 * mu * gamma) / (m * c^2))
  
  error <- mu + c(-1, 1) * sd * 1.96
  
  list(mu = mu, sd = sd, conf = error)
}
