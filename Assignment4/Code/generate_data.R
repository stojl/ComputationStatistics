generate_x <- function(n, sd) {
  exp(rnorm(n, mean = 0, sd = sd))
}

generate_y <- function(x, sigma, alpha = 0, beta = 1, gamma = 1, rho = 1) {
  create_regfunc <- function(alpha, beta, gamma, rho) {
    force(alpha)
    force(beta)
    force(gamma)
    force(rho)
    
    function(x) {
      gamma + (rho - gamma) / (1 + exp(beta * log(x) - alpha))
    }
  }
  
  regfunc <- create_regfunc(alpha, beta, gamma, rho)
  regfunc(x) + rnorm(length(x), mean = 0, sd = sigma)
}

create_regfunc <- function(alpha, beta, gamma, rho) {
  force(alpha)
  force(beta)
  force(gamma)
  force(rho)
  
  function(x) {
    gamma + (rho - gamma) / (1 + exp(beta * log(x) - alpha))
  }
}
