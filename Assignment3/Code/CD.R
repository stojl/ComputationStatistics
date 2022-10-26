CG <- function(
    par, 
    H,
    gr,
    d = 0.8, 
    c = 0.1, 
    gamma0 = 1, 
    epsilon = 1e-7,
    maxiter = 500,
    cb = NULL
) {
  p <- length(par)
  m <- 1
  rho0 <- numeric(p)
  for(i in 1:maxiter) {
    value <- H(par)
    grad <- gr(par)
    grad_norm_sq <- sum(grad^2)
    gamma <- gamma0
    rho <- - grad + grad_norm_sq * rho0
    h_prime <- drop(t(grad) %*% rho)
    if(m > p || h_prime >= 0) {
      rho <- - grad
      h_prime <- - grad_norm_sq 
      m <- 1
    }
    par1 <- par + gamma * rho
    while(H(par1) > value + c * gamma * h_prime) {
      gamma <- d * gamma
      par1 <- par + gamma * rho
    }
    if(!is.null(cb)) cb()
    rho0 <- rho / grad_norm_sq
    if(norm(par - par1, "2") < epsilon * (norm(par, "2") + epsilon)) {
      par <- par1
      break
    } 
    par <- par1
    m <- m + 1
  }
  
  if(i == maxiter)
    warning("Maximal number, ", maxiter, ", of iterations reached")
  par
}