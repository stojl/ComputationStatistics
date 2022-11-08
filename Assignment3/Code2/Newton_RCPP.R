Newton_cpp <- function(par, x, nu, d = 0.8, c = 0.2, gamma0 = 1, min.eps = 1e-7, maxit = 500) {
  result <- FitT(par, x, nu, c, d, gamma0, maxit, min.eps)
  if(result[[2]] == maxit)
    warning("Maximal number, ", maxit, ", of iterations reached")
  list(par = result[[1]], iterations = result[[2]])
}
