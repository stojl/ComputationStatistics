set.seed(34929392)
library(ggplot2)
library(dplyr)
library(magrittr)

test_robustness <- function(mu, sigma, nu, n) {
  result <- vector("list", n)
  initial_par <- numeric(2 * n)
  dim(initial_par) <- c(n, 2)
  X <- extraDistr::rlst(2000, df = nu, mu = mu, sigma = sqrt(sigma))
  for(i in 1:n) {
    mu_r <- rcauchy(1, location = mu, scale = 100)
    sigma_r <- extraDistr::rpareto(1, a = 0.15, b = sigma)
    initial_par[i, 1] <- mu_r
    initial_par[i, 2] <- sigma_r
    result[[i]] <- EM(c(mu_r, sigma_r), X, nu)
  }
  list(results = result, initial_par = initial_par)
}

convert_to_df <- function(df) {
  do.call(rbind, lapply(df$results, 
                        function(x) c(x$par, iterations = x$iterations))) %>% 
    as.data.frame()
}

set.seed(34929392)
test_2.5 <- test_robustness(0, 3, 2.5, 250)
test_0.5 <- test_robustness(0, 3, 0.5, 250)
test_30 <- test_robustness(0, 3, 30, 250)

result_df_2.5 <- convert_to_df(test_2.5)
result_df_0.5 <- convert_to_df(test_0.5) 
result_df_30 <- convert_to_df(test_30)

massive_result <- bind_rows(
  result_df_0.5,
  result_df_2.5,
  result_df_30
) %>% 
  mutate(nu = factor(nu, levels = c("0.5", "2.5", "30")))

massive_result %>% 
  ggplot(aes(x = mu, y = sigma, col = nu)) +
  geom_point() +
  geom_point(aes(x = 0, y = 3), size = 2, col = "red") +
  theme_bw()

massive_result %>% 
  ggplot(aes(y = nu, x = iterations)) +
  geom_boxplot() +
  xlab("Number of iterations to reach convergence criteria") +
  ylab("Degrees of freedom") +
  theme_bw()
