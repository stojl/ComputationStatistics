source("rng_stream.R")

ruin_prob_4_MC <- function(n, m, theta, batch_size, sd_threshold) {
  
  get_samples <- rng_stream(batch_size, runif, -1.9, 2)
  
  quantile_function <- invFunc(theta)  
  current_sd <- sd_threshold
  
  w_star <- numeric(m)
  
  while(current_sd >= sd_threshold) {
    
    samples <- quantile_function(get_samples(batch_size))
    
    dim(samples) <- c(n, m)
    
    
    
    for(i in seq_along(w_star)) {
      w_star[i] <- exp(-theta * sum(samples[, i]))
    }
    
    w <- w_star / sum(w_star)
    
  }
    
}