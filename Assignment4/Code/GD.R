GD <- function(par, 
               H, 
               gr, 
               d = 0.8, 
               c = 0.1, 
               gamma0 = 1, 
               epsilon = 1e-7, 
               maxiter = 500,
               backtrack = TRUE,
               cb = NULL
) {
  for(i in 1:maxiter) {
    if(backtrack) value <- H(par)
    grad <- gr(par)
    h_prime <- sum(grad^2)
    gamma <- gamma0
    par1 <- par - gamma * grad
    if(!is.null(cb)) cb()
    if(backtrack) {
      while(H(par1) > value - c * gamma * h_prime) {
        gamma <- d * gamma
        par1 <- par - gamma * grad
      }
    }
    if(norm(par - par1, "2") < epsilon * (norm(par, "2") + epsilon)) {
      break
    } 
    par <- par1
  }
  if(i == maxiter)
    warning("Maximal number, ", maxiter, ", of iterations reached")
  par1
}
