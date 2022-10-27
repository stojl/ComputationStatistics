epoch_adam <- function(batch_size = 1, alpha = 0.95, beta = 0.9) {
  rho <- v <- 0
  function(par0, x, y, loss_gr, gamma) {
    M <- floor(length(x) / batch_size)
    par <- par0
    for(i in 1:M) {
      batch_indicies <- ((i - 1) * batch_size + 1):(i * batch_size)
      gr <- loss_gr(par, x[batch_indicies], y[batch_indicies])
      rho <<- alpha * rho + (1 - alpha) * gr
      v <<- beta * v + (1 - beta) * gr^2
      par <- par - gamma * (rho / (sqrt(v) + 1e-8))
    }
    par
  }
}

epoch_momentum <- function(batch_size = 1, beta = 0.95) {
  rho <- 0
  function(par0, x, y, loss_gr, gamma) {
    M <- floor(length(x) / batch_size)
    par <- par0
    for(i in 1:M) {
      batch_indicies <- ((i - 1) * batch_size + 1):(i * batch_size)
      gr <- loss_gr(par, x[batch_indicies], y[batch_indicies])
      rho <<- beta * rho + (1 - beta) * gr
      par <- par - gamma * rho
    }
    par
  }
}

epoch_batch <- function(batch_size = 1) {
  function(par0, x, y, loss_gr, gamma) {
    M <- floor(length(x) / batch_size)
    par <- par0
    for(i in 1:M) {
      batch_indicies <- ((i - 1) * batch_size + 1):(i * batch_size)
      gr <- loss_gr(par, x[batch_indicies], y[batch_indicies])
      par <- par - gamma * gr
    }
    par
  }
}

if(!exists("epoch_rcpp")) Rcpp::sourceCpp("Rcpp/Gradient.cpp")
epoch_batch_rcpp <- function(batch_size) {
  function(par0, x, y, loss_gr, gamma) {
    epoch_rcpp(par0, x, y, batch_size, gamma)
  }
}

epoch_batch_rcpp_partial <- function(batch_size) {
  function(par0, x, y, loss_gr, gamma) {
    epoch_rcpp_partial(par0, x, y, loss_gr, batch_size, gamma)
  }
}

epoch_full_batch <- function(iterations = 1) {
  function(par0, x, y, loss_gr, gamma) {
    for(i in 1:iterations) {
      gr <- loss_gr(par0, x, y)
      par0 <- par0 - gamma * gr
    }
    par0
  }
}

epoch_full_momentum <- function(iterations = 1, beta = 0.9) {
  rho <- 0
  function(par0, x, y, loss_gr, gamma) {
    for(i in 1:iterations) {
      gr <- loss_gr(par0, x, y)
      rho <<- beta * rho + (1 - beta) * gr
      par0 <- par0 - gamma * rho
    }
    par0
  }
}

epoch_full_adam <- function(iterations = 1, alpha = 0.95, beta = 0.9) {
  rho <- v <- 0
  function(par0, x, y, loss_gr, gamma) {
    for(i in 1:iterations) {
      gr <- loss_gr(par0, x, y)
      rho <<- alpha * rho + (1 - alpha) * gr
      v <<- beta * v + (1 - beta) * gr^2
      par0 <- par0 - gamma * (rho / (sqrt(v) + 1e-8))
    }
    par0
  }
}