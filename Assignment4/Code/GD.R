GD <- function(par0, 
               H, 
               gr, 
               d = 0.8, 
               c = 0.1, 
               gamma0 = 1, 
               epsilon = 1e-7, 
               maxit = 500,
               backtrack = TRUE,
               cb = NULL
) {
  for(i in 1:maxit) {
    if(backtrack) value <- H(par0)
    grad <- gr(par0)
    h_prime <- sum(grad^2)
    gamma <- gamma0
    par <- par0 - gamma * grad
    if(!is.null(cb)) cb()
    if(backtrack) {
      while(H(par) > value - c * gamma * h_prime) {
        gamma <- d * gamma
        par <- par0 - gamma * grad
      }
    }
    if(norm(par0 - par, "2") < epsilon * (norm(par0, "2") + epsilon)) {
      break
    } 
    par0 <- par
  }
  if(i == maxit)
    warning("Maximal number, ", maxit, ", of iterations reached")
  par
}
