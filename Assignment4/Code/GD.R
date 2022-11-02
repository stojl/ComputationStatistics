GD <- function(par0, 
               loss, 
               loss_gr, 
               d = 0.8, 
               c = 0.1, 
               gamma0 = 1, 
               maxit = 500,
               stop_criteria = NULL,
               backtrack = TRUE,
               cb = NULL
) {
  for(i in 1:maxit) {
    if(backtrack) value <- loss(par0)
    grad <- loss_gr(par0)
    h_prime <- sum(grad^2)
    gamma <- gamma0
    par <- par0 - gamma * grad
    if(!is.null(cb)) cb()
    if(backtrack) {
      while(loss(par) > value - c * gamma * h_prime) {
        gamma <- d * gamma
        par <- par0 - gamma * grad
      }
    }
    if(!is.null(stop_criteria)) {
      if(stop_criteria(par, par0, loss_gr, loss)) break
    }
    par0 <- par
  }
  if(i == maxit)
    warning("Maximal number, ", maxit, ", of iterations reached")
  par
}
