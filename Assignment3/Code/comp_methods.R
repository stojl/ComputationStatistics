source("EM-algorithm.R")
source("hess_grad.R")
source("GD.R")
source("CD.R")
source("Newton.R")
library(tibble)
library(dplyr)
library(ggplot2)
set.seed(3939392)
nu1 <- 0.5
X <- extraDistr::rlst(10000000, df = nu1, mu = 5, sigma = sqrt(2))

grad_func <- gradl(X, nu1)
func <- logl(X, nu1)
hess_func <- hessl(X, nu1)

GD_tracer <- CSwR::tracer(c("par", "conv"), N = 0, expr = quote(conv <- norm(par - par1, "2")))
GD_result <- 
  GD(c(median(X), IQR(X)), func, grad_func, gamma0 = 1, d = 0.8, c = 0.1, cb = GD_tracer$tracer)
GD_summary <- summary(GD_tracer)

EM_tracer <- CSwR::tracer(c("par", "conv"), N = 0, expr = quote(conv <- norm(par - par1, "2")))
EM_result <-
  EM(par = c(median(X), IQR(X)), x = X, nu = nu1, cb = EM_tracer$tracer)
EM_summary <- summary(EM_tracer)

CG_tracer <- CSwR::tracer(c("par", "conv"), N = 0, expr = quote(conv <- norm(par - par1, "2")))
CG_result <-
  CG(c(median(X), IQR(X)), func, grad_func, gamma0 = 1, d = 0.8, c = 0.1, cb = CG_tracer$tracer)
CG_summary <- summary(CG_tracer)

NT_tracer <- CSwR::tracer(c("par", "conv"), N = 0, expr = quote(conv <- norm(par - par1, "2")))
NT_result <-
  Newton(c(median(X), IQR(X)), func, grad_func, hess_func, gamma0 = 1, d = 0.8, c = 0.1, cb = NT_tracer$tracer)
NT_summary <- summary(NT_tracer)

results <- bind_rows(EM_summary, GD_summary, CG_summary, NT_summary, .id = "Alg") %>% 
  mutate(Algorithm = factor(Alg, labels = c("EM", "GD", "CGD", "NT")))

results %>% 
  ggplot(aes(x = .time, y = par.1, col = Algorithm)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  xlab("Time [s]") +
  ylab("mu")

results %>% 
  ggplot(aes(x = .time, y = par.2, col = Algorithm)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  xlab("Time [s]") +
  ylab("sigma^2")

results %>% 
  ggplot(aes(x = .time, y = conv, col = Algorithm)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  scale_y_log10() +
  xlab("Time [s]") +
  ylab("norm(par - par1, \"2\")")

library(profvis)
init <- c(median(X), IQR(X))
profvis(Newton(init, func, grad_func, hess_func, gamma0 = 1, d = 0.8, c = 0.1), interval = 0.005)
profvis(EM(par = init, x = X, nu = nu1), interval = 0.005)
