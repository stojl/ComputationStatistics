library(ggplot2)
library(dplyr)
library(tibble)
library(extraDistr)
library(tidyr)
library(magrittr)
plotframe <- 
  tibble(x = seq(30, 60, 0.5)) %>% 
  mutate(df1 = dlst(x, 1, 45, sqrt(3)),
         df0.2 = dlst(x, 0.2, 45, sqrt(3)),
         df2 = dlst(x, 2, 45, sqrt(3)),
         df30 = dlst(x, 30, 45, sqrt(3)),
         df100 = dlst(x, 100, 45, sqrt(3)),
         norm = dnorm(x, 45, sqrt(3))) %>% 
  pivot_longer(-x, names_to = "dist", 
               values_to = "value", 
               names_transform = list(dist = ~factor(.x, levels = c("df0.2",
                                                            "df1",
                                                            "df2",
                                                            "df30",
                                                            "df100",
                                                            "norm"))))

plotframe %>% 
  ggplot(aes(x = x, y = value, col = as.factor(dist))) +
  geom_line(size = 1) +
  theme_bw() +
  xlim(c(37, 53))
  