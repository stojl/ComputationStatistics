batch_random <- function(N, ...) {
  function(x, y) {
    rand_indicies <- sample(seq_along(x), N, ...)
    x_rand <- x[rand_indicies]
    y_rand <- y[rand_indicies]
    
    list(x = x_rand, y = y_rand)
  }
}

batch_random_chunk <- function(N, shuffle = FALSE) {
  x_shuffle <- NULL
  y_shuffle <- NULL
  max_chunk <- 1
  force(shuffle)
  function(x, y) {
    if(shuffle & is.null(x_shuffle)) {
      rand_indicies <- sample(seq_along(x))
      x_shuffle <<- x[rand_indicies]
      y_shuffle <<- y[rand_indicies]
      max_chunk <<- max(1, floor(length(x_shuffle) / N) - 1)
    } else {
      x_shuffle <<- x
      y_shuffle <<- y
      max_chunk <<- max(1, floor(length(x_shuffle) / N) - 1)
    }
    chunk <- sample(1:max_chunk, 1)
    x_rand <- x_shuffle[((chunk - 1) * N + 1):(chunk * N)]
    y_rand <- y_shuffle[((chunk - 1) * N + 1):(chunk * N)]
    list(x = x_rand, y = y_rand)
  }
}

batch_linear_chunk <- function(N, shuffle = TRUE) {
  chunk <- 1
  x_shuffle <- NULL
  y_shuffle <- NULL
  force(shuffle)
  function(x, y) {
    if(shuffle & is.null(x_shuffle)) {
      rand_indicies <- sample(seq_along(x))
      x_shuffle <<- x[rand_indicies]
      y_shuffle <<- y[rand_indicies]
    } else {
      x_shuffle <<- x
      y_shuffle <<- y
    }
    max_index <- floor(length(x) / N)
    if(chunk != max_index) {
      x_rand <- x_shuffle[((chunk - 1) * N + 1):(chunk * N)]
      y_rand <- y_shuffle[((chunk - 1) * N + 1):(chunk * N)]
      chunk <<- chunk + 1
    } else {
      chunk <<- 1
      x_rand <- x_shuffle[((chunk - 1) * N + 1):(chunk * N)]
      y_rand <- y_shuffle[((chunk - 1) * N + 1):(chunk * N)]
    }
    list(x = x_rand, y = y_rand)
  }
}