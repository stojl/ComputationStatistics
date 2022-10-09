running_variance <- function(df = 1) {
  force(df)
  
  m <- 0
  SSQ <- 0
  S <- 0
  
  function(x = NULL) {
    
    if(!is.null(x)) {
      m <<- m + length(x)
      S <<- S + sum(x)
      SSQ <<- SSQ + sum(x^2)
    }
    
    V1 <- SSQ / (m - df)
    V2 <- S^2 / (m * (m - df))
    
    V1 - V2
  }
}

running_mean <- function() {
  m <- 0
  S <- 0
  
  function(x = NULL) {
    
    if(!is.null(x)) {
      m <<- m + length(x)
      S <<- S + sum(x)
    }
    
    S / m
  }
}

running_covariance <- function(df = 1) {
  
  force(df)
  m <- 0
  SP <- 0
  S1 <- 0
  S2 <- 0
  
  function(x = NULL, y = NULL) {
    
    if(length(x) != length(y)) {
      stop("Length of vectors are not equal.")
    }
    
    if(!is.null(x) & !is.null(y)) {
      m <<- m + length(x)
      SP <<- SP + sum(x * y)
      S1 <<- S1 + sum(x)
      S2 <<- S2 + sum(y)
    }
    
    CV1 <- SP / (m - df)
    CV2 <- S1 * S2 / (m * (m - df))
    
    CV1 - CV2
  }
}
