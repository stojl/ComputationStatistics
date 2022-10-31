epoch_adam <- function(mini_batch_size = 1, alpha = 0.95, beta = 0.9) {
  rho <- v <- 0
  function(par0, index, loss_gr, gamma) {
    mini_batch_size <- min(length(index), mini_batch_size)
    M <- floor(length(index) / mini_batch_size)
    par <- par0
    for(i in 1:M) {
      mini_batch_index <- ((i - 1) * mini_batch_size + 1):(i * mini_batch_size)
      gr <- loss_gr(par, index[mini_batch_index])
      rho <<- alpha * rho + (1 - alpha) * gr
      v <<- beta * v + (1 - beta) * gr^2
      par <- par - gamma * (rho / (sqrt(v) + 1e-8))
    }
    par
  }
}

epoch_momentum <- function(mini_batch_size = 1, beta = 0.95) {
  rho <- 0
  function(par0, index, loss_gr, gamma) {
    mini_batch_size <- min(length(index), mini_batch_size)
    M <- floor(length(index) / mini_batch_size)
    par <- par0
    for(i in 1:M) {
      mini_batch_index <- ((i - 1) * mini_batch_size + 1):(i * mini_batch_size)
      gr <- loss_gr(par, index[mini_batch_index])
      rho <<- beta * rho + (1 - beta) * gr
      par <- par - gamma * rho
    }
    par
  }
}

epoch_batch <- function(mini_batch_size = 1) {
  function(par0, index, loss_gr, gamma) {
    mini_batch_size <- min(length(index), mini_batch_size)
    M <- floor(length(index) / mini_batch_size)
    par <- par0
    for(i in 1:M) {
      mini_batch_index <- ((i - 1) * mini_batch_size + 1):(i * mini_batch_size)
      gr <- loss_gr(par, index[mini_batch_index])
      par <- par - gamma * gr
    }
    par
  }
}

if(!exists("epoch_rcpp")) Rcpp::sourceCpp("Rcpp/Gradient.cpp")
epoch_batch_rcpp <- function(batch_size, x, y) {
  function(par0, index, loss_gr, gamma) {
    epoch_rcpp(par0, x[index], y[index], batch_size, gamma)
  }
}

epoch_batch_rcpp_partial <- function(batch_size, x, y) {
  function(par0, index, loss_gr, gamma) {
    epoch_rcpp_partial(par0, x[index], y[index], loss_gr, batch_size, gamma)
  }
}

epoch_full <- function() {
  function(par0, index, loss_gr, gamma) {
    for(i in 1:length(index)) {
      gr <- loss_gr(par0, index[i])
      par0 <- par0 - gamma * gr
    }
    par0
  }
}

epoch_full_momentum <- function(beta = 0.9) {
  rho <- 0
  function(par0, index, loss_gr, gamma) {
    for(i in 1:length(index)) {
      gr <- loss_gr(par0, index[i])
      rho <<- beta * rho + (1 - beta) * gr
      par0 <- par0 - gamma * rho
    }
    par0
  }
}

epoch_full_adam <- function(alpha = 0.95, beta = 0.9) {
  rho <- v <- 0
  function(par0, index, loss_gr, gamma) {
    for(i in 1:length(i)) {
      gr <- loss_gr(par0, index[i])
      rho <<- alpha * rho + (1 - alpha) * gr
      v <<- beta * v + (1 - beta) * gr^2
      par0 <- par0 - gamma * (rho / (sqrt(v) + 1e-8))
    }
    par0
  }
}