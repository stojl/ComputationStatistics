Rcpp::sourceCpp("src/CPP_GradientAndHessian.cpp")

get_like_cpp <- function(x, nu) {
  force(x); force(nu)
  function(par) CPP_likelihood(par, x, nu)
}

get_grad_cpp <- function(x, nu) {
  force(x); force(nu)
  function(par) CPP_gradient(par, x, nu)
}

get_hess_cpp <- function(x, nu) {
  force(x); force(nu)
  function(par) CPP_hessian(par, x, nu)
}