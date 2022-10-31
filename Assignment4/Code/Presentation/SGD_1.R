SGD_1 <- function(par0,
                  loss_gr,
                  N,
                  gamma0 = 1,
                  maxit = 15,
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
  par <- par0
  for(i in 1:maxit) {
    index <- sample(N)
    for(j in 1:N) {
      gr <- loss_gr(par, index[j])
      par <- par - gamma[i] * gr
    }
    if(!is.null(cb)) cb()
    par0 <- par
  }
  par
}
