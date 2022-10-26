Newton <- function(par, H, gr, hess, d = 0.8, c = 0.2, gamma0 = 1, epsilon = 1e-7, maxiter = 500, cb = NULL) {
  for(i in 1:maxiter) {
    value <- H(par)
    grad <- gr(par)
    Hessian <- hess(par) 
    rho <- - drop(solve(Hessian, grad)) 
    gamma <- gamma0
    par1 <- par + gamma * rho
    h_prime <- t(grad) %*% rho 
    while(min(H(par1), Inf, na.rm = TRUE) > value +  c * gamma * h_prime) { 
      gamma <- d * gamma 
      par1 <- par + gamma * rho
    }
    if(!is.null(cb)) cb()
    if(norm(par - par1, "2") < epsilon * (norm(par, "2") + epsilon)) break 
    par <- par1 
  }
  if(i == maxiter)
    warning("Maximal number, ", maxiter, ", of iterations reached")
  par1
}
