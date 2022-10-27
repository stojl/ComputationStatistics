SGD <- function(x,
                y,
                batch = NULL,
                par0,
                loss,
                loss_gr,
                gamma0 = 1,
                epoch = NULL,
                ...,
                maxit = 15,
                cb = NULL) {
  if(length(x) != length(y)) {
    stop("Length of x: ", length(x), 
         " is not equal to length of y: ", 
         length(y), ".")
  }
  if(is.numeric(gamma0)) {
    if(length(gamma0) == 1) {
      gamma <- rep(gamma0, maxit)
    } else {
      gamma <- c(gamma0, rep(gamma0[length(gamma0)], maxit - length(gamma0)))
    }
  } else if (is.function(gamma0)) {
    gamma <- gamma0(1:maxit)
  } else {
    stop("gamma0 must be a numeric or a function.")
  }
  for(i in 1:maxit) {
    xy <- batch(x, y)
    par <- epoch(par0, xy$x, xy$y, loss_gr, gamma[i])
    if(!is.null(cb)) cb()
    par0 <- par
  }
  par
}
