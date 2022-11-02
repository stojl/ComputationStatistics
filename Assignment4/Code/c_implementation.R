dyn.load("c/src/sgd_o.o")

SGD2 <- function(par0,
                 loss_gr,
                 N,
                 batch,
                 epoch,
                 gamma0 = 1,
                 maxit = 15,
                 loss) {
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
  
  .Call("sgd", 
        par0, 
        loss_gr, 
        N, 
        batch, 
        epoch,
        gamma,
        maxit,
        environment())
}

epoch_test <- function(mini_batch_size = 1) {
  function(par0, index, loss_gr, gamma0) {
    .Call("epoch_batch",
          par0,
          index,
          loss_gr,
          gamma0,
          as.integer(mini_batch_size),
          environment())
  }
}
