if(!exists("gradient_rcpp")) Rcpp::sourceCpp("Rcpp/Gradient.cpp")

log_gradient_rcpp <- function(x, y) {
  force(x); force(y)
  
  function(par, index) {
    gradient_rcpp(par, x[index], y[index])
  }
}
