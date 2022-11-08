EM <- function(x, nu, cb = NULL, maxit = 500, min.eps = 1e-7, par = NULL) {
  if(is.null(par)) par <- c(median(x), IQR(x))
  par1 <- par
  if(!is.null(cb)) cb()
  n <- length(x)
  EW <- numeric(n)
  for(i in 1:maxit) {
    EW <- (nu + 1) / (1 + ((x - par[1])^2) / (nu * par[2]))
    par1[1] <- sum(EW * x) / sum(EW)
    par1[2] <- sum(EW * (x - par[1])^2) / (n * nu)
    if(!is.null(cb)) cb()
    if(norm(par - par1, "2") < min.eps * (norm(par, "2") + min.eps)) break
    par <- par1
  }
  if(i == maxit) warning("Maximum number of itertaions ", maxit, " reached.")
  names(par1) <- c("mu", "sigma")
  list(par = c(par1, nu = nu), iterations = i)
}
