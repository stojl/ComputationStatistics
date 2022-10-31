log_loss <- function(x, y) {
  force(x); force(y)
  
  function(par) {
    mean((y - par[3] - (par[4] - par[3]) / 
            (1 + exp(par[2] * log(x) - par[1])))^2, na.rm = TRUE)
  }
}

log_gradient <- function(x, y) {
  force(x); force(y)
  function(par, index) {
    elogx <- exp(par[2] * log(x[index]) - par[1])
    
    drho <- 1 / (1 + elogx)
    dgamma <- 1 - drho
    dalpha <- elogx * (par[4] - par[3]) * drho^2
    dbeta <- -dalpha * log(x[index])
    
    diff <- y[index] - par[3] - (par[4] - par[3]) * drho
    
    c(-mean(diff * dalpha, na.rm = TRUE) / 2,
      -mean(diff * dbeta, na.rm = TRUE) / 2,
      -mean(diff * dgamma, na.rm = TRUE) / 2,
      -mean(diff * drho, na.rm = TRUE) / 2)
  }
}
