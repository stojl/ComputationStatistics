source("./Iteration1/ruin_prob_1.R")
source("./Iteration2/ruin_prob_2.R")
source("./Iteration3/ruin_prob_3.R")
ggplot2::theme_set(ggplot2::theme_bw())
library(magrittr)

profvis(ruin_prob_1(100, 1000000))
profvis(ruin_prob_2(100, 1000000))

microbenchmark::microbenchmark(
  ruin_prob_1(100, 100000),
  ruin_prob_2(100, 100000)
)

ruin <- ruin_prob_sum_1(100, 500000)



data.frame(x = 1:length(ruin), y = ruin) %>% 
  ggplot(aes(x = x, y = y)) + 
  geom_line()
