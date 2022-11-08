if(!is.loaded("src/dens.o")) {
  dyn.load("src/dens.o")
}

dens <- function(x,
                 kernel,
                 bandwidth,
                 points = 512L,
                 ...) {
  if(is.function(bandwidth)) {
    bw <- bandwidth(x, kernel, ...)
  } else {
    bw <- bandwidth
  }
  p <- seq(min(x) - bw, max(x) + bw, length.out = points)
  y <- .Call("C_dens", x, p, kernel, bw, environment())
  structure(
    list(
      x = p,
      y = y,
      bw = bw
    ),
    class = "dens"
  )
}

dens_cpp_partial <- function(x,
                             kernel,
                             bandwidth,
                             points = 512L,
                             ...) {
  if(is.function(bandwidth)) {
    bw <- bandwidth(x, kernel, ...)
  } else {
    bw <- bandwidth
  }
  p <- seq(min(x) - bw, max(x) + bw, length.out = points)
  y <- dens_rcpp_partial(x, kernel, bw, p)
  structure(
    list(
      x = p,
      y = y,
      bw = bw
    ),
    class = "dens"
  )
}

dens_cpp <- function(x,
                     kernel,
                     bandwidth,
                     points = 512L,
                     ...) {
  if(is.function(bandwidth)) {
    bw <- bandwidth(x, kernel, ...)
  } else {
    bw <- bandwidth
  }
  p <- seq(min(x) - bw, max(x) + bw, length.out = points)
  y <- dens_rcpp(x, bw, p)
  structure(
    list(
      x = p,
      y = y,
      bw = bw
    ),
    class = "dens"
  )
}

source("R_dens.R")
dens_R <- function(x,
                   kernel,
                   bandwidth,
                   points = 512L,
                   ...) {
  if(is.function(bandwidth)) {
    bw <- bandwidth(x, kernel, ...)
  } else {
    bw <- bandwidth
  }
  p <- seq(min(x) - bw, max(x) + bw, length.out = points)
  y <- R_dens(x, p, kernel, bw)
  structure(
    list(
      x = p,
      y = y,
      bw = bw
    ),
    class = "dens"
  )
}

dens_R1 <- function(x,
                   kernel,
                   bandwidth,
                   points = 512L,
                   ...) {
  if(is.function(bandwidth)) {
    bw <- bandwidth(x, kernel, ...)
  } else {
    bw <- bandwidth
  }
  p <- seq(min(x) - bw, max(x) + bw, length.out = points)
  y <- R_dens_for(x, p, kernel, bw)
  structure(
    list(
      x = p,
      y = y,
      bw = bw
    ),
    class = "dens"
  )
}
