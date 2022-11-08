test1 <- adap_samp(100, 
                   poisdens, 
                   poisdens_derv, 
                   c(0.1, 0.2, 0.3),
                   seed = 1234)
test2 <- adap_samp_cpp(100, 
                       poisdens, 
                       poisdens_derv, 
                       c(0.1, 0.2, 0.3),
                       seed = 1234)

bm1 <- mark(
  QUANT_R = test1$quant(seq(0.01, 1, 0.01)),
  QUANT_RCPP = test2$quant(seq(0.01, 1, 0.01)),
  DENS_R = test1$dens(seq(0, 0.5, 0.01)),
  DENS_RCPP = test2$dens(seq(0, 0.5, 0.01)),
  check = FALSE,
  min_time = 10
)

autoplot(bm1)

bm2 <- mark(
  DENS_RCPP = poisdens_cpp(seq(0, 1, 0.01)),
  DENS_R = poisdens(seq(0, 1, 0.01)),
  DENS_DERIV_CPP = poisdens_derv_cpp(seq(0, 1, 0.01)),
  DENS_DERIV_R = poisdens_derv(seq(0, 1, 0.01)),
  check = FALSE,
  min_time = 10
)
autoplot(bm2)