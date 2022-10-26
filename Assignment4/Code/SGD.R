SGD <- function(x,
                y,
                N,
                par0,
                loss,
                loss_gr,
                gamma0 = 1,
                epoch = NULL,
                ...,
                maxit = 15,
                cb = NULL) {
  
  if(is.numeric(gamma0)) {
    if(length(gamma0) == 1) {
      gamma <- function(x) gamma0
    } else {
      gamma <- function(x) gamma0[1]
      warning("gamma0 has length > 1. First value is chosen.")
    }
  }
  
  par <- par0
  
  for(i in 1:maxit) {
    rand_indicies <- sample(seq_along(x), N)
    x_rand <- x[rand_indicies]
    y_rand <- y[rand_indicies]
    par <- epoch(par0, x_rand, y_rand, loss_gr, gamma)
    if(!is.null(cb)) cb()
    par0 <- par
  }
  par0
}
