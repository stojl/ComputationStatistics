if(!is.loaded("src/test.o")) {
  dyn.load("src/test.o")
}
Newton_C <- function(par, H, gr, hess, d = 0.8, c = 0.2, gamma0 = 1, min.eps = 1e-7, maxit = 500) {
  result <- .Call("C_Newton",
                  par,
                  H,
                  gr,
                  hess,
                  d,
                  c,
                  gamma0,
                  min.eps,
                  as.integer(maxit),
                  environment())
  if(result[[2]] == -1L)
    warning("Matrix solve went wrong.")
  if(result[[2]] == maxit)
    warning("Maximal number, ", maxit, ", of iterations reached")
  list(par = result[[1]], iterations = result[[2]])
}