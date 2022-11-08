get_env_quantile<- function(a, b, z) {
  force(a); force(b); force(z)
  az <- a * z[-length(z)]
  R <- exp(b) * (exp(a * z[-1]) - exp(az)) / a
  Q1 <- numeric(length(a) + 1)
  Q1[2:length(Q1)] <- cumsum(R)
  c <- Q1[length(Q1)]
  function(q) {
    ind <- c * q <= Q1
    maxi <- which.max(ind) - 1
    y <- c * q - Q1[maxi]
    log(a[maxi] * y * exp(-b[maxi]) + exp(az[maxi])) / a[maxi]
  }
}

get_env_density <- function(a, b, z) {
  force(a); force(b); force(z)
  function(x) {
    if(x > z[length(z)] | x < z[1]) return(0)
    maxi <- which.max(x <= z) - 1
    exp(a[maxi] * x + b[maxi])
  }
}

adap_samp <- function(n, density, density_deriv, p, zb = c(-Inf, Inf), seed = NULL) {
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
  env_quantile <- get_env_quantile(a, b, z)
  env_density <- get_env_density(a, b, z)
  samples <- numeric(n)
  succes <- tries <- 0
  for(s in 1:n) {
    reject <- TRUE
    while(reject) {
      tries <- tries + 1
      u0 <- runif(2)
      y0 <- env_quantile(u0[1])
      env_y0 <- env_density(y0)
      dens_y0 <- density(y0)
      if(u0[2] <= dens_y0 / env_y0) {
        reject <- FALSE
        succes <- succes + 1
        samples[s] <- y0
      }
    }
  }
  list(samples, (tries - succes) / tries)
}
