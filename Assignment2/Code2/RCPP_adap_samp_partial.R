Rcpp::sourceCpp("RCPP_PARTIAL.cpp")

adap_samp_cpp_partial <- function(n, 
                                  density, 
                                  density_deriv, 
                                  p, 
                                  zb = c(-Inf, Inf),
                                  seed = NULL) {
  if(!is.null(seed)) set.seed(seed)
  p <- sort(unique(p))
  densp <- density(p)
  a <- density_deriv(p) / densp
  b <- log(densp) - a * p
  a_diff <- a[-length(a)] - a[-1]
  check1 <- a[1] < 0 & zb[1] == -Inf
  check2 <- a[length(a)] > 0 & zb[2] == Inf
  if(check1 | check2)
    stop("Envelope is not integrable. Choose different points.")
  if(any(a == 0) | any(a_diff == 0))
    stop("Divison by zero. Choose different points.")
  z <- c(zb[1], (b[-1] - b[-length(b)]) / a_diff, zb[2])
  
  RCPP_adap_samp_partial(n, density, a, b ,z)
}
