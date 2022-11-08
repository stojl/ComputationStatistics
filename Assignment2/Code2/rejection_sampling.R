rejection_sampling <- function(n, 
                               density, 
                               env_density, 
                               env_sampler, 
                               alpha,
                               seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  samples <- numeric(n)
  succes <- tries <- 0
  for(s in 1:n) {
    reject <- TRUE
    while(reject) {
      tries <- tries + 1
      u0 <- runif(1)
      y0 <- env_sampler()
      env_y0 <- env_density(y0)
      dens_y0 <- density(y0)
      if(u0 <= alpha * dens_y0 / env_y0) {
        reject <- FALSE
        samples[s] <- y0
        succes <- succes + 1
      }
    }
  }
  list(samples, (tries - succes) / tries)
}
