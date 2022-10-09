source("./Iteration1/ruin_prob_1.R")
source("./Iteration2/ruin_prob_2.R")
source("./Iteration3/ruin_prob_3.R")
source("./Iteration4/ruin_prob_4.R")
source("./Iteration5/ruin_prob_5.R")

library(magrittr)
library(bench)
library(ggplot2)
library(profvis)

ggplot2::theme_set(ggplot2::theme_bw())

profvis(ruin_prob_1(100, 1000000))
profvis(ruin_prob_2(100, 1000000))

microbenchmark::microbenchmark(
  ruin_prob_1(100, 100000),
  ruin_prob_2(100, 100000)
)

ruin <- ruin_prob_sum_1(100, 500000)

profvis::profvis()

microbenchmark::microbenchmark(ruin_prob_4(100, -0.13, 500, 0.00005),
                               ruin_prob_5(100, 500, 0.00005))

profvis(ruin_prob_4(100, -0.05, 500, 0.00005))
profvis(ruin_prob_5(100, 500, 0.00005))

dist <- numeric(1000L)

for(i in seq_along(dist)) {
  dist[i] <- ruin_prob_5(100, 500, 0.00005)$mu
}


result <- mark(IS = ruin_prob_4(100, -0.13, 500, 0.00005),
     MC = ruin_prob_5(100, 500, 0.00005),
     min_time = 60, check = FALSE)
autoplot(result)

