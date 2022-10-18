library(ggplot2)
library(dplyr)
library(tibble)
library(extraDistr)
library(tidyr)
library(magrittr)
plotframe <- 
  tibble(x = seq(30, 60, 0.5)) %>% 
  mutate(df2 = dlst(x, 2, 45, sqrt(3)),
         df30 = dlst(x, 30, 45, sqrt(3)),
         df100 = dlst(x, 100, 45, sqrt(3)),
         norm = dnorm(x, 45, sqrt(3))) %>% 
  pivot_longer(-x, names_to = "dist", values_to = "value")

plotframe %>% 
  ggplot(aes(x = x, y = value, col = dist)) +
  geom_line(size = 1) +
  theme_bw() +
  xlim(c(37, 53))
  