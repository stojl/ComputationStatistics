# Code for the log dose response curve function
source("generate_data.R")
source("LDR.R")

# Implementation of Stochastic Gradient Descent
source("SGD.R")
source("epochs.R")
source("batch.R")
source("decay.R")

library(CSwR)
library(profvis)
library(ggplot2)

X <- generate_x(10000000, 2)
Y <- generate_y(X, 0.5, alpha = 1, beta = 3, gamma = 5, rho = 1)

SG_tracer <- tracer("")

SGD(
  par0 = c(3,2,6,3),
  loss_gr = log_gradient_rcpp(X, Y),
  N = length(X),
  batch = batch_random_chunk(100000, 2),
  epoch = epoch_momentum(10000),
  gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.5, n1 = 1000),
  maxit = 30
)

bm <- mark(
  RCPP = SGD(
    par0 = c(3,2,6,3),
    loss_gr = log_gradient_rcpp(X, Y),
    N = length(X),
    batch = batch_random_chunk(500),
    epoch = epoch_batch(100),
    gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 1000),
    maxit = 500
  ),
  R = SGD(
    par0 = c(3,2,6,3),
    loss_gr = log_gradient(X, Y),
    N = length(X),
    batch = batch_random_chunk(500),
    epoch = epoch_batch(100),
    gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 1000),
    maxit = 500
  ),
  min_time = 20,
  check = FALSE
)
autoplot(bm)


profvis::profvis(
  SGD(
    par0 = c(3,2,6,3),
    loss_gr = log_gradient_rcpp(X, Y),
    N = length(X),
    batch = batch_random_chunk(100000, Inf),
    epoch = epoch_momentum(10000),
    gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.5, n1 = 1000),
    maxit = 10
  )
)
