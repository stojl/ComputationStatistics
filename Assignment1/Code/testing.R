source("dens.R")
source("bandwidth.R")
source("R_dens.R")
library(ggplot2)
# INITIALIZATION
e_kernel <- function(x) {
  0.75 * (1 - x^2) * (abs(x) <= 1)
}

g_kernel <- function(x) {
  exp(-x^2 / 2) * 0.39894228040143
}

x <- c(rnorm(200, -1, 0.1), 
       rnorm(200, -0.5, 0.1),
       rnorm(200, 0, 0.3),
       rnorm(200, 0.5, 0.7),
       rnorm(200, 1, 1))

funky_density <- function(x) {
  dnorm(x, -1, 0.1) / 5 +
    dnorm(x, -0.5, 0.1) / 5 +
    dnorm(x, 0, 0.3) / 5 +
    dnorm(x, 0.5, 0.7) / 5 +
    dnorm(x, 1, 1) / 5
}
# FIRST TEST OF IMPLEMENTATION
p2 <- seq(max(x) - 0.1, min(x) + 0.1, length.out = 1024L)
y1 <- R_dens_for(x, p2, e_kernel, 0.1)
test3 <- density(x, kernel = "epanechnikov", bw = 0.1 / sqrt(5))
theme_set(theme_bw())
ggplot(data.frame(x = p2, y = y1, col = "Density Estimate"), aes(x = x, y = y1)) +
  geom_line(size = 1) +
  geom_line(aes(x = x3, y = y3, col = "True density"), data = data.frame(x3 = seq(min(x), max(x), 0.01),
                                                                       y3 = funky_density(seq(min(x), max(x), 0.01))), size = 1) +
  geom_line(aes(x = x5, y = y5, col = "R density"), data = data.frame(x5 = test3$x, y5 = test3$y), size = 0.5)

# TEST OF IMPLEMENTATION
test <- dens_R1(x, e_kernel, bw_cv_R2, points = 1024L)
test1 <- dens_R1(x, e_kernel, bw_oracle, points = 1024L)
test2 <- density(x, kernel = "epanechnikov", bw = "ucv")
theme_set(theme_bw())
ggplot(data.frame(x = test$x, y = test$y), aes(x = x, y = y, col = "CV")) +
  geom_line(size = 1) +
  geom_line(aes(x = z, y = w, col = "density (ucv)"), data = data.frame(z = test2$x, w = test2$y), size = 0.5) +
  geom_line(aes(x = x, y = y, col = "True density"), data = data.frame(x = seq(min(x), max(x), 0.01),
                                                                       y = funky_density(seq(min(x), max(x), 0.01))), size = 1)  +
  geom_line(aes(x = s, y = t, col = "Oracle"), data = data.frame(s = test1$x, t = test1$y), size = 1)
# PROFILING OF IMPLEMENTATION
library(profvis)
profvis(dens_R1(x, e_kernel, bw_cv_R2, points = 1024L))

# DIFFERENT IMPLEMENTATIONS OF DENSITY CALCULATION IN R
p1 <- seq(min(x), max(x), length.out = 1024L)
library(bench)
library(dplyr)
library(magrittr)
bm <- mark(
  ONE_FOR_LOOP = R_dens(x, p1, e_kernel, 0.5),
  TWO_FOR_LOOP = R_dens_for(x, p1, e_kernel, 0.5),
  APPLY = R_dens_apply(x, p1, e_kernel, 0.5),
  check = FALSE,
  min_time = 20
)
autoplot(bm)
select(bm, expression, median, mem_alloc) %>% 
  arrange(median)

# DIFFERENT IMPLEMENTATIONS OF CV CALCULATION IN R
bm2 <- mark(
  TRIANGLE = bw_cv_R(x, e_kernel),
  RECTANGLE = bw_cv_R2(x, e_kernel),
  check = FALSE,
  min_time = 20
)
autoplot(bm2)
select(bm2, expression, median, mem_alloc) %>% 
  arrange(median)

# RCPP IMPLEMENTATION OF CV
bm3 <- mark(
  R = bw_cv_R(x, e_kernel),
  CPP = bw_cv_cpp_partial(x, e_kernel),
  R_CPP_KERNEL = bw_cv_R(x, e_kernel_cpp),
  CPP_CPP_KERNEL = bw_cv_cpp_partial(x, e_kernel_cpp),
  CPP_FULL = bw_cv_cpp(x, e_kernel_cpp),
  check = FALSE,
  min_time = 20
)
autoplot(bm3)
select(bm3, expression, median, mem_alloc) %>% 
  arrange(median)

# RCPP IMPLEMENTATION OF DENSITY CALCULATION
bm4 <- mark(
  R = R_dens(x, p1, e_kernel, 0.5),
  CPP = dens_rcpp_partial(x, e_kernel, 0.5, p1),
  R_CPP_KERNEL = R_dens(x, p1, e_kernel_cpp, 0.5),
  CPP_CPP_KERNEL = dens_rcpp_partial(x, e_kernel_cpp, 0.5, p1),
  CPP_FULL = dens_rcpp(x, 0.5, p1),
  check = FALSE,
  min_time = 20
)
autoplot(bm4)
select(bm4, expression, median, mem_alloc) %>% 
  arrange(median)

# C IMPLEMENTATION THROUGH R'S C API BANDWIDTH
bm5 <- mark(
  C = bw_cv(x, e_kernel),
  R = bw_cv_R(x, e_kernel),
  CPP = bw_cv_cpp_partial(x, e_kernel),
  R_CPP_KERNEL = bw_cv_R(x, e_kernel_cpp),
  C_CPP_KERNEL = bw_cv(x, e_kernel_cpp),
  CPP_CPP_KERNEL = bw_cv_cpp_partial(x, e_kernel_cpp),
  CPP_FULL = bw_cv_cpp(x, e_kernel_cpp),
  check = FALSE,
  min_time = 20
)
autoplot(bm5)
select(bm5, expression, median, mem_alloc) %>% 
  arrange(median)

# C IMPLEMENTATION THROUGH R'S C API DENSITY CALCULATION
bm6 <- mark(
  C = .Call("C_dens", x, p1, e_kernel, 0.5, environment()),
  R = R_dens(x, p1, e_kernel, 0.5),
  CPP = dens_rcpp_partial(x, e_kernel, 0.5, p1),
  C_CPP_KERNEL = .Call("C_dens", x, p1, e_kernel_cpp, 0.5, environment()),
  R_CPP_KERNEL = R_dens(x, p1, e_kernel_cpp, 0.5),
  CPP_CPP_KERNEL = dens_rcpp_partial(x, e_kernel_cpp, 0.5, p1),
  CPP_FULL = dens_rcpp(x, 0.5, p1),
  check = FALSE,
  min_time = 20
)
autoplot(bm6)
select(bm6, expression, median, mem_alloc) %>% 
  arrange(median)
