source("rng_stream.R")
source("./Iteration3/ruin_prob_3.R")
source("running_statistics.R")

ruin_prob_5 <- function(n, batch_size, sd_threshold, min_samples = 2000) {
  
  current_sd <- sd_threshold
  
  results <- logical(batch_size)
  
  sigma_MC <- running_variance()
  mu <- running_mean()
  m_current <- 0
  
  while(current_sd >= sd_threshold | m_current <= min_samples) {
    
    samples <- runif(batch_size * n, -1.9, 2)
    
    m_current <- m_current + batch_size
    
    dim(samples) <- c(n, batch_size)
    
    results <- ruin_vector(samples, n, batch_size)
    
    sigma_MC(results)
    mu(results)
    
    estimate <- mu()
    
    sd <- sqrt(sigma_MC() / m_current)
    current_sd <- sd
    error <- estimate + c(-1, 1) * sd * 1.96
  }
  
  list(mu = estimate, 
       sd = sd, 
       conf = error, 
       samples = m_current, 
       epochs = m_current / batch_size, 
       batch_size = batch_size)
}
