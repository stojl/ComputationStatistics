get_grad <- function(par, x, y) {
  alpha <- par[1]
  beta <- par[2]
  gamma <- par[3]
  rho <- par[4]
  
  n <- length(x)
  
  lx <- log(x)
  ex <- exp(beta * lx - alpha)
  
  drho <- 1 / (1 + ex)
  dgamma <- 1 - drho
  dalpha <- (rho - gamma) * ex / drho^2
  dbeta <- -dalpha * lx
  
  ss <- y - gamma - (rho - gamma) * drho
  
  drho <- sum(ss * drho)
  dgamma <- sum(ss * dgamma)
  dalpha <- sum(ss * dalpha)
  dbeta <- sum(ss * dbeta)
  
  -2 * c(dalpha, dbeta, dgamma, drho) / n
}


