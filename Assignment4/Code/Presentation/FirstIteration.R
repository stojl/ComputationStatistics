source("generate_data.R")
source("Presentation/SGD_1.R")
library(ggplot2)
library(magrittr)

X <- generate_x(1000, 2)
Y <- generate_y(X, 0.5, alpha = 1, beta = 3, gamma = 5, rho = 1)

data.frame(x = X, y = Y) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_point(shape = 1, alpha = 0.5) +
  geom_smooth(se = FALSE) +
  scale_x_log10() +
  theme_bw()

source("LDR.R")
source("decay.R")
library(CSwR)

opt_value <- optim(c(1, 3, 5, 1), log_loss(X, Y))$value

SGD_tracer <- tracer("diff", expr = quote(diff <- loss(par) - opt_value))
SGD_1(c(3, 6, 3, 3), 
      log_gradient(X, Y),
      N = length(X),
      gamma0 = decay_scheduler(gamma0 = 0.1, gamma1 = 0.01, n1 = 200),
      maxit = 200,
      loss = log_loss(X, Y),
      cb = SGD_tracer$tracer)

summary(SGD_tracer) %>% 
  ggplot(aes(x = .time, y = diff)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_bw()

source("SGD.R")
source("epochs.R")
source("batch.R")

SGD_tracer2 <- tracer("diff", expr = quote(diff <- loss(par) - opt_value))
SGD(par0 = c(3, 6, 3, 3),
    loss_gr = log_gradient(X, Y),
    N = length(X),
    batch = batch_random(replace = FALSE),
    epoch = epoch_batch(500),
    gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 200),
    loss = log_loss(X, Y),
    maxit = 200,
    cb = SGD_tracer2$tracer)

summary(SGD_tracer2) %>% 
  ggplot(aes(x = .time, y = diff)) +
  geom_point() +
  geom_line() +
  geom_point(aes(x = .time, y = diff), data = summary(SGD_tracer), col = "red") +
  geom_line(aes(x = .time, y = diff), data = summary(SGD_tracer), col = "red") +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw()

SGD_tracer3 <- tracer("diff", expr = quote(diff <- loss(par) - opt_value))
SGD(par0 = c(3, 6, 3, 3),
    loss_gr = log_gradient(X, Y),
    N = length(X),
    batch = batch_random(replace = FALSE),
    epoch = epoch_full(),
    gamma0 = decay_scheduler(gamma0 = 0.1, gamma1 = 0.01, n1 = 200),
    maxit = 200,
    loss = log_loss(X, Y),
    cb = SGD_tracer3$tracer)

summary(SGD_tracer3) %>% 
  ggplot(aes(x = .time, y = diff)) +
  geom_point() +
  geom_line() +
  geom_point(aes(x = .time, y = diff), data = summary(SGD_tracer), col = "red") +
  geom_line(aes(x = .time, y = diff), data = summary(SGD_tracer), col = "red") +
  scale_y_log10() +
  theme_bw()

library(tidyr)
library(bench)

bm <- mark(SGD = SGD(par0 = c(3, 6, 3, 3),
               loss_gr = log_gradient(X, Y),
               N = length(X),
               batch = batch_random(replace = FALSE),
               epoch = epoch_full(),
               gamma0 = decay_scheduler(gamma0 = 0.1, gamma1 = 0.01, n1 = 200),
               maxit = 10),
           SGD_1 = SGD_1(c(3, 6, 3, 3), 
                         log_gradient(X, Y),
                         N = length(X),
                         gamma0 = decay_scheduler(gamma0 = 0.1, gamma1 = 0.01, n1 = 200),
                         maxit = 10),
           min_time = 60,
           check = FALSE)

autoplot(bm)

source("GD.R")
library(dplyr)

get_summary <- function(x, loss, value) {
  summary(x) %>% 
    rowwise() %>% 
    mutate(diff = loss(c(par.1, par.2, par.3, par.4)) - value)
}

GD_tracer <- tracer("par", N = 0)
GD(par0 = c(3, 6, 3, 3),
   H = log_loss(X, Y),
   gr = log_gradient(X, Y),
   gamma0 = 1,
   maxit = 1500,
   backtrack = FALSE,
   cb = GD_tracer$tracer)

GD_sum <- get_summary(GD_tracer, log_loss(X, Y), opt_value)

GD_sum %>% 
  ggplot(aes(x = .time, y = diff)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_bw()

SGD_tracer4 <- tracer("par", N = 0)
SGD(par0 = c(3, 6, 3, 3),
    loss_gr = log_gradient(X, Y),
    N = length(X),
    batch = batch_random(replace = FALSE),
    epoch = epoch_batch(200),
    gamma0 = decay_scheduler(gamma0 = 1, gamma1 = 0.2, n1 = 500),
    loss = log_loss(X, Y),
    maxit = 1000,
    cb = SGD_tracer4$tracer)


SGD_sum <- get_summary(SGD_tracer4, log_loss(X, Y), opt_value)

SGD_sum %>% 
  ggplot(aes(x = .time, y = diff)) +
  geom_point() +
  geom_line() +
  geom_point(aes(x = .time, y = diff), data = GD_sum, col = "red") +
  geom_line(aes(x = .time, y = diff), data = GD_sum, col = "red") +
  scale_y_log10() +
  theme_bw()

X_big <- generate_x(1000000, 2)
Y_big <- generate_y(X_big, 0.5, alpha = 1, beta = 3, gamma = 5, rho = 1)

opt_parm <- optim(c(1, 3, 5, 1), log_loss(X_big, Y_big))
opt_value_big <- opt_parm$value

GD_tracer2 <- tracer("par", N = 0)
GD(par0 = c(3, 6, 3, 2),
   H = log_loss(X_big, Y_big),
   gr = log_gradient(X_big, Y_big),
   gamma0 = 3,
   maxit = 1000,
   backtrack = FALSE,
   cb = GD_tracer2$tracer)

GD_sum2 <- get_summary(GD_tracer2, log_loss(X_big, Y_big), opt_value_big)

GD_sum2 %>% 
  ggplot(aes(x = .time, y = diff)) +
  geom_point() +
  geom_line() +
  scale_y_log10() +
  theme_bw()

SGD_tracer5 <- tracer("par", N = 0)
SGD(par0 = c(3, 6, 3, 3),
    loss_gr = log_gradient(X_big, Y_big),
    N = length(X_big),
    batch = batch_random(10000, replace = TRUE),
    epoch = epoch_batch(200),
    gamma0 = decay_scheduler(gamma0 = 3, gamma1 = 0.04, n1 = 1000),
    maxit = 1000,
    cb = SGD_tracer5$tracer)

SGD_sum2 <- get_summary(SGD_tracer5, log_loss(X_big, Y_big), opt_value_big)
SGD_sum2 %>% 
  ggplot(aes(x = .time, y = diff)) +
  geom_point() +
  geom_line() +
  geom_point(aes(x = .time, y = diff), data = GD_sum2, col = "red") +
  geom_line(aes(x = .time, y = diff), data = GD_sum2, col = "red") +
  scale_y_log10() +
  scale_x_log10() +
  theme_bw()


source("RCPP_gradient.R")
X_CPP <- generate_x(10000, 2)
Y_CPP <- generate_y(X_CPP, 0.5, alpha = 1, beta = 3, gamma = 5, rho = 1)

opt_parm_cpp <- GD(par0 = c(1, 3, 5, 1),
                    H = log_loss(X_CPP, Y_CPP),
                    gr = log_gradient_rcpp(X_CPP, Y_CPP),
                    gamma0 = 3,
                    maxit = 10000,
                    backtrack = TRUE)

opt_value_cpp <- log_loss(X_CPP, Y_CPP)(opt_parm_cpp)

SGD_tracer6 <- tracer("par", N = 0)
SGD(par0 = c(3, 6, 3, 3),
    loss_gr = log_gradient_rcpp(X_CPP, Y_CPP),
    N = length(X_CPP),
    batch = batch_random(replace = TRUE),
    epoch = epoch_batch(500),
    gamma0 = decay_scheduler(gamma0 = 3, gamma1 = 0.04, n1 = 500),
    maxit = 500,
    cb = SGD_tracer6$tracer)

SGD_sum3 <- get_summary(SGD_tracer6, log_loss(X_CPP, Y_CPP), opt_value_cpp)

SGD_tracer7 <- tracer("par", N = 0)
SGD(par0 = c(3, 6, 3, 3),
    loss_gr = log_gradient(X_CPP, Y_CPP),
    N = length(X_CPP),
    batch = batch_random(replace = TRUE),
    epoch = epoch_batch(500),
    gamma0 = decay_scheduler(gamma0 = 3, gamma1 = 0.04, n1 = 500),
    maxit = 500,
    cb = SGD_tracer7$tracer)

SGD_sum4 <- get_summary(SGD_tracer7, log_loss(X_CPP, Y_CPP), opt_value_cpp)

SGD_sum4 %>% 
  ggplot(aes(x = .time, y = diff)) +
  geom_point() +
  geom_line() +
  geom_point(aes(x = .time, y = diff), data = SGD_sum3, col = "red") +
  geom_line(aes(x = .time, y = diff), data = SGD_sum3, col = "red") +
  scale_y_log10() +
  theme_bw()

bm2 <- mark(
  R = SGD(par0 = c(3, 6, 3, 3),
          loss_gr = log_gradient(X_CPP, Y_CPP),
          N = length(X_CPP),
          batch = batch_random(replace = TRUE),
          epoch = epoch_batch(500),
          gamma0 = decay_scheduler(gamma0 = 3, gamma1 = 0.04, n1 = 500),
          maxit = 500),
  CPP_GR = SGD(par0 = c(3, 6, 3, 3),
               loss_gr = log_gradient(X_CPP, Y_CPP),
               N = length(X_CPP),
               batch = batch_random(replace = TRUE),
               epoch = epoch_batch(500),
               gamma0 = decay_scheduler(gamma0 = 3, gamma1 = 0.04, n1 = 500),
               maxit = 500),
  CPP_EPOCH = SGD(par0 = c(3, 6, 3, 3),
                  loss_gr = log_gradient(X_CPP, Y_CPP),
                  N = length(X_CPP),
                  batch = batch_random(replace = TRUE),
                  epoch = epoch_batch_rcpp(500, X_CPP, Y_CPP),
                  gamma0 = decay_scheduler(gamma0 = 3, gamma1 = 0.04, n1 = 500),
                  maxit = 500),
  CPP_PARTIAL = SGD(par0 = c(3, 6, 3, 3),
                    loss_gr = log_gradient(X_CPP, Y_CPP),
                    N = length(X_CPP),
                    batch = batch_random(replace = TRUE),
                    epoch = epoch_batch_rcpp_partial(500, X_CPP, Y_CPP),
                    gamma0 = decay_scheduler(gamma0 = 3, gamma1 = 0.04, n1 = 500),
                    maxit = 500),
  CPP_FULL = ...,
  check = FALSE,
  min_time = 30
)

autoplot(bm)
