source("SGD.R")
source("epochs.R")
source("func_grad/loss_grad.R")
source("Newton_SGD.R")
source("generate_data.R")
source("batch.R")
source("decay.R")
library(CSwR)
library(ggplot2)

X <- generate_x(1000000, 2)
Y <- generate_y(X, 0.5, alpha = 1, beta = 3, gamma = 5, rho = 1)

X1 <- X[order(X)]
Y1 <- Y[order(X)]

nt <- get_hessian_grad(X, Y)
nt_loss <- get_loss(X, Y)
opt_gr <- get_gradient(X, Y)

optim(par = c(3, 1, 8, 4), fn = nt_loss, gr = opt_gr, method = "BFGS", control = list(maxit = 500))
nt_tracer <- CSwR::tracer("diff", expr = quote(diff <- H(par)))
Newton(c(3, 1, 6, 4), nt_loss, nt, cb = nt_tracer$tracer)


GD(par = c(3, 1, 8, 4), nt_loss, opt_gr, gamma0 = 1, maxiter = 500, backtrack = FALSE)

SG_trace <- tracer("diff", N = 1, expr = quote(diff <- loss(par0, x, y)))
SGD(X, Y, batch = batch_random_chunk(500, shuffle = FALSE), c(3, 1, 8, 4), 
    loss_function, loss_gradient, 
    gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 2000), 
    epoch = epoch_batch_rcpp_partial(100), 
    cb = NULL, 
    maxit = 500)
SG_sum <- summary(SG_trace)
autoplot(SG_sum, diff, log = TRUE)


microbenchmark::microbenchmark(
  SGD(X, Y, batch = batch_random_chunk(500, shuffle = FALSE), c(3, 1, 8, 4), 
      loss_function, loss_gradient, 
      gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 2000), 
      epoch = epoch_batch_rcpp(100), 
      cb = NULL, 
      maxit = 500),
  
  SGD(X, Y, batch = batch_random_chunk(500, shuffle = FALSE), c(3, 1, 8, 4), 
      loss_function, loss_gradient, 
      gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 2000), 
      epoch = epoch_batch(100), 
      cb = NULL, 
      maxit = 500)
)

result <- mark(
  CPP = SGD(X, Y, batch = batch_random_chunk(500, shuffle = FALSE), c(3, 1, 8, 4), 
      loss_function, loss_gradient, 
      gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 2000), 
      epoch = epoch_batch_rcpp(100), 
      cb = NULL, 
      maxit = 500),
  
  R_GR = SGD(X, Y, batch = batch_random_chunk(500, shuffle = FALSE), c(3, 1, 8, 4), 
            loss_function, loss_gradient, 
            gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 2000), 
            epoch = epoch_batch_rcpp_partial(100), 
            cb = NULL, 
            maxit = 500),
  
  CPP_GR = SGD(X, Y, batch = batch_random_chunk(500, shuffle = FALSE), c(3, 1, 8, 4), 
            loss_function, gradient, 
            gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 2000), 
            epoch = epoch_batch(100), 
            cb = NULL, 
            maxit = 500),
  
  R = SGD(X, Y, batch = batch_random_chunk(500, shuffle = FALSE), c(3, 1, 8, 4), 
      loss_function, loss_gradient, 
      gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 2000), 
      epoch = epoch_batch(100), 
      cb = NULL, 
      maxit = 500),
  min_time = 20,
  check = FALSE
)
autoplot(result)
