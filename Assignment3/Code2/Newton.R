Newton <- function(par, H, gr, hess, d = 0.8, c = 0.2, gamma0 = 1, min.eps = 1e-7, maxit = 500, cb = NULL) {
  par1 <- par
  if(!is.null(cb)) cb()
  for(i in 1:maxit) {
    value <- H(par)
    grad <- gr(par)
    Hessian <- hess(par) 
    rho <- - drop(solve(Hessian, grad)) 
    gamma <- gamma0
    par1 <- par + gamma * rho
    h_prime <- crossprod(grad, rho)
    while(min(H(par1), Inf, na.rm = TRUE) > value +  c * gamma * h_prime) { 
      gamma <- d * gamma 
      par1 <- par + gamma * rho
    }
    if(!is.null(cb)) cb()
    if(norm(par - par1, "2") < min.eps * (norm(par1, "2") + min.eps)) break 
    par <- par1 
  }
  if(i == maxit) warning("Maximal number, ", maxit, ", of iterations reached")
  list(par = par1, iterations = i)
}
