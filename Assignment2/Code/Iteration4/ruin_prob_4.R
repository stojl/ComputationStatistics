source("rng_stream.R")
source("./Iteration3/ruin_prob_3.R")
source("running_statistics.R")

ruin_prob_4 <- function(n, theta, batch_size, sd_threshold, min_samples = 2000) {
  
  quantile_function <- invFunc(theta)  
  current_sd <- sd_threshold
  
  w_star <- numeric(batch_size)
  results <- logical(batch_size)
  
  sigma_IS <- running_variance()
  gamma <- running_covariance()
  sigma_w <- running_variance()
  c_mean <- running_mean()
  norm_sum <- 0
  mu <- 0
  m_current <- 0
  
  while(current_sd >= sd_threshold | m_current <= min_samples) {
    
    samples <- runif(batch_size * n)
    
    samples <- quantile_function(samples)
    
    m_current <- m_current + batch_size
    
    dim(samples) <- c(n, batch_size)
    
    w_star <- exp(-theta * colSums(samples))
    
    results <- ruin_vector(samples, n, batch_size)
    
    sigma_IS(results * w_star)
    gamma(results * w_star, w_star)
    sigma_w(w_star)
    c_mean(w_star)
    mu <- mu + sum(results * w_star)
    norm_sum <- norm_sum + sum(w_star)

    estimate <- mu / norm_sum
    sd <- sqrt((sigma_IS() + estimate^2 * sigma_w() - 2 * estimate * gamma()) / (m_current * c_mean()^2))
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
