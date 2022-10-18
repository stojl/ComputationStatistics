inv_digamma <- function(x) {
  sel_x <- x >= -2.22
  y <- ifelse(x >= -2.22, exp(x + 0.5), -1 / (x - digamma(1)))
  for(i in 1:3) {
    y <- y - (digamma(y) - x) / trigamma(y)
  }
  y
}

E_step <- function(x) {
  force(x)
  function(par) {
    C1 <- digamma((par[3] + 1) / 2) - log(1/2 + ((x - par[1])^2) / (2 * par[3] * par[2]))
    C2 <- (par[3] + 1) / (1 + ((x - par[1])^2) / (par[3] * par[2]))
    list(C1, C2)
  }
}

M_step <- function(x) {
  force(x)
  function(C) {
    mu <- sum(C[[2]] * x) / sum(C[[2]])
    nu <- 2 * inv_digamma(mean(C[[1]]) - log(2))
    sigma <- mean(C[[2]] * (x - mu)^2) / nu
    c(mu, sigma, nu)
  }
}

EM_3 <- function(par0, x, maxit = 15, min.eps = 1e-7) {
  E <- E_step(x)
  M <- M_step(x)
  par <- par0
  for(i in 1:maxit) {
    C <- E(par)
    new_par <- M(C)
    if(norm(new_par - par, "2") < min.eps * (norm(new_par, "2") + min.eps)) {
      par <- par
      break
    }
    par <- new_par
    if(i == maxit) warning("Maximum itertaions reached.")
  }
  names(par) <- c("mu", "sigma", "nu")
  list(par = par, iterations = i)
}

X <- extraDistr::rlst(100000, df = 100, mu = 3000, sigma = sqrt(15))
EM_3(c(500,10,80), X, maxit = 10000)
