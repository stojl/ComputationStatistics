Rcpp::sourceCpp("./Iteration2/ruin_prob_2.cpp")

ruin_prob_2 <- function(n, m) {
  
  samples <- runif(n * m, min = -1.9, 2)
  
  ruin_prob_cpp(samples, n, m) / m
}
