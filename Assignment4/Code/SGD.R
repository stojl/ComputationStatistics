SGD <- function(par0,
                loss_gr,
                N,
                batch,
                epoch,
                gamma0 = 1,
                maxit = 15,
                stop_criteria = NULL,
                loss = NULL,
                cb = NULL) {
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
    index <- batch(N)
    par <- epoch(par0, index, loss_gr, gamma[i])
    if(!is.null(cb)) cb()
    if(!is.null(stop_criteria)) {
      if(stop_criteria(par, par0, loss_gr, loss)) break
    }
    par0 <- par
  }
  if(!is.null(stop_criteria) & maxit == i) {
    warning("Maximum number of iterations ", i, " reached.")
  }
  par
}
