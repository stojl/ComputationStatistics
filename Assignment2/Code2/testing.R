library(ggplot2)
library(magrittr)
library(dplyr)
library(tidyr)
library(bench)
library(profvis)
theme_set(theme_bw())
library(gridExtra)
source("likelihood.R")


source("rejection_sampling.R")
source("adap_samp.R")
source("C_adap_samp.R")
source("RCPP_adaptive_sampling.R")
source("adap_samp_stream.R")
source("RCPP_adap_samp_partial.R")
norm_constant <- integrate(poisdens, 0, Inf)$value
mupois <- optimize(function(x) -poisdens(x), c(0, 0.5))$minimum
sigmapois <- integrate(function(x) (x - mupois)^2 * poisdens(x) / norm_constant, 0, Inf)$value
gaus1 <- gausdens(mupois, sqrt(sigmapois) / 14)

# PLOT OF UNNORMALIZED DENSITY
p1 <- data.frame(x = seq(0, 0.5, 0.001), y = poisdens(seq(0, 0.5, 0.001)) / (0.74 * 1e-40)) %>% 
  ggplot(aes(x = x, y = y)) +
  geom_line(size = 1) +
  geom_function(fun = gaus1, col = "red") +
  labs(title = "Unnormalized Density of Poisson Regression")

p2 <- data.frame(x = seq(0, 2.6, 0.001),
           y = gaus1(seq(0, 2.6, 0.001)) - 
             poisdens(seq(0, 2.6, 0.001)) / (0.74 * 1e-40)) %>% 
  ggplot(aes(x = x, y = y)) + 
  ylab("Difference") +
  geom_line(size = 1)

grid.arrange(p1, p2)

# GAUSSIAN ENVELOPE
envelope_gauss <- function(t) 0.74 * 1e-40 * gaus1(t)
all(envelope_gauss(seq(0, 2, 0.001)) - poisdens(seq(0, 2, 0.001)) >= 0)

adap_test <- adap_samp(100000, 
                  poisdens, 
                  poisdens_derv, c(0.1, 0.2, 0.3))

# IT WORKS!
data.frame(x = adap_test) %>% 
  ggplot(aes(x = x)) + 
  geom_density(aes(col = "Sample Density"), 
               size = 1) + 
  geom_function(aes(col = "True Density"), 
                fun = function(x) poisdens(x) / norm_constant, 
                size = 0.5)

gauss_test <- rejection_sampling(100000,
                   density = poisdens,
                   env_density = gaus1,
                   env_sampler = function() rnorm(1, mupois, sqrt(sqrt(sigmapois) / 14)),
                   alpha = 1 / (0.74 * 1e-40),
                   seed = 12348)

test_c <- adap_samp(100000,
                      density = poisdens,
                      density_deriv = poisdens_derv,
                      p = c(0.15, 0.2, 0.28, 0.32),
                      seed = 12348) 

data.frame(x = test_c[[1]]) %>% 
  ggplot(aes(x = x)) + 
  geom_density(aes(col = "Adap Sample Density"), 
               size = 1) + 
  geom_density(aes(x = x1, col = "Rej. Sample Density"), 
               data = data.frame(x1 = gauss_test[[1]])) +
  geom_function(aes(col = "True Density"), 
                fun = function(x) poisdens(x) / norm_constant, 
                size = 0.5)

bm0 <- mark(
  GAUSS = rejection_sampling(100,
                             density = poisdens,
                             env_density = gaus1,
                             env_sampler = function() rnorm(1, mupois, sqrt(sqrt(sigmapois) / 14)),
                             alpha = 1 / (0.74 * 1e-40),
                             seed = 12348),
  ADAP = adap_samp(100, 
                  poisdens, 
                  poisdens_derv, 
                  p = c(0.15, 0.2, 0.28, 0.32),
                  seed = 12348),
  check = FALSE,
  min_time = 10
)
autoplot(bm0)

profvis(adap_samp(100000, 
                  poisdens, 
                  poisdens_derv, 
                  p = c(0.15, 0.2, 0.28, 0.32),
                  seed = 12348))

bm1 <- mark(
  GAUSS = rejection_sampling(100,
                             density = poisdens_cpp,
                             env_density = gaus1,
                             env_sampler = function() rnorm(1, mupois, sqrt(sqrt(sigmapois) / 14)),
                             alpha = 1 / (0.74 * 1e-40),
                             seed = 12348),
  ADAP= adap_samp(100, 
                  poisdens, 
                  poisdens_derv, 
                  c(0.1, 0.18, 0.25, 0.3),
                  seed = 12348),
  ADAP_STR = adap_samp_stream(100, 
                              poisdens, 
                              poisdens_derv, 
                              c(0.1, 0.18, 0.25, 0.3),
                              seed = 12348),
  ADAP_STR_CPP = adap_samp_cpp(100, 
                              poisdens, 
                              poisdens_derv, 
                              c(0.1, 0.18, 0.25, 0.3),
                              seed = 12348),
  check = FALSE,
  min_time = 10)

bm1

bm2 <- mark(
  GAUSS = rejection_sampling(100,
                             density = poisdens_cpp,
                             env_density = gaus1,
                             env_sampler = function() rnorm(1, mupois, sqrt(sqrt(sigmapois) / 14)),
                             alpha = 1 / (0.74 * 1e-40),
                             seed = 12348),
  ADAP= adap_samp(100, 
                  poisdens, 
                  poisdens_derv, 
                  c(0.1, 0.18, 0.25, 0.3),
                  seed = 12348),
  ADAP_STR = adap_samp_stream(100, 
                              poisdens, 
                              poisdens_derv, 
                              c(0.1, 0.18, 0.25, 0.3),
                              seed = 12348),
  ADAP_STR_CPP = adap_samp_cpp(100, 
                               poisdens, 
                               poisdens_derv, 
                               c(0.1, 0.18, 0.25, 0.3),
                               seed = 12348),
  ADAP_PARTIAL_CPP = adap_samp_cpp_partial(100, 
                               poisdens, 
                               poisdens_derv, 
                               c(0.1, 0.18, 0.25, 0.3),
                               seed = 12348),
  check = FALSE,
  min_time = 10)

bm2

bm3 <- mark(
  GAUSS = rejection_sampling(100,
                             density = poisdens_cpp,
                             env_density = gaus1,
                             env_sampler = function() rnorm(1, mupois, sqrt(sqrt(sigmapois) / 14)),
                             alpha = 1 / (0.74 * 1e-40),
                             seed = 12348),
  ADAP= adap_samp(100, 
                  poisdens, 
                  poisdens_derv, 
                  c(0.1, 0.18, 0.25, 0.3),
                  seed = 12348),
  ADAP_STR = adap_samp_stream(100, 
                              poisdens, 
                              poisdens_derv, 
                              c(0.1, 0.18, 0.25, 0.3),
                              seed = 12348),
  ADAP_STR_CPP = adap_samp_cpp(100, 
                               poisdens, 
                               poisdens_derv, 
                               c(0.1, 0.18, 0.25, 0.3),
                               seed = 12348),
  ADAP_PARTIAL_CPP = adap_samp_cpp_partial(100, 
                                           poisdens, 
                                           poisdens_derv, 
                                           c(0.1, 0.18, 0.25, 0.3),
                                           seed = 12348),
  ADAP_C = adap_samp_c(100, 
                       poisdens, 
                       poisdens_derv, 
                       c(0.1, 0.18, 0.25, 0.3),
                       seed = 12348),
  check = FALSE,
  min_time = 10)

bm3


create_benchmark <- function(n_samples) {
  bm_func <- function(n) {
    mark(
      GAUSS = rejection_sampling(n,
                                 density = poisdens_cpp,
                                 env_density = gaus1,
                                 env_sampler = function() rnorm(1, mupois, sqrt(sqrt(sigmapois) / 14)),
                                 alpha = 1 / (0.74 * 1e-40),
                                 seed = 12348),
      ADAP= adap_samp(n, 
                      poisdens, 
                      poisdens_derv, 
                      c(0.1, 0.2, 0.3),
                      seed = 12348),
      ADAP_STR = adap_samp_stream(n, 
                                  poisdens, 
                                  poisdens_derv, 
                                  c(0.1, 0.2, 0.3),
                                  seed = 12348),
      ADAP_STR_CPP = adap_samp_cpp(n, 
                                   poisdens, 
                                   poisdens_derv, 
                                   c(0.1, 0.2, 0.3),
                                   seed = 12348),
      ADAP_PARTIAL_CPP = adap_samp_cpp_partial(n, 
                                               poisdens, 
                                               poisdens_derv, 
                                               c(0.1, 0.2, 0.3),
                                               seed = 12348),
      ADAP_C = adap_samp_c(n, 
                           poisdens, 
                           poisdens_derv, 
                           c(0.1, 0.2, 0.3),
                           seed = 12348),
      check = FALSE,
      time_unit = 'ms',
      min_iterations = 100) %>% 
      select(expression, median, mem_alloc) %>% 
      mutate(samples = n)
  }
  bms <- lapply(n_samples, bm_func)
  do.call(bind_rows, bms)
}

many_benchmarks <- create_benchmark(c(64, 128, 256, 512, 1024, 2048, 4096, 8192))

fac_name <- attributes(many_benchmarks$expression)$description

many_benchmarks %>% 
  mutate(name = as.factor(fac_name)) %>% 
  ggplot(aes(x = samples, y = median, col = name)) +
  geom_line() +
  ylab("Median Runtime [milliseconds]")


test_c <- adap_samp(100000,
                    density = poisdens,
                    density_deriv = poisdens_derv,
                    p = c(0.15, 0.2, 0.28, 0.32)) 

data.frame(x = test_c[[1]]) %>% 
  ggplot(aes(x = x)) + 
  geom_density(aes(col = "C Adap Sample Density"), 
               size = 1) + 
  geom_function(aes(col = "True Density"), 
                fun = function(x) poisdens(x) / norm_constant, 
                size = 0.5)
