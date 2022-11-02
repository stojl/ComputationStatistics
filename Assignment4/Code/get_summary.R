library(dplyr)

get_summary <- function(x, loss, value) {
  summary(x) %>% 
    rowwise() %>% 
    mutate(diff = loss(c(par.1, par.2, par.3, par.4)) - value)
}