batch_random <- function(batch_size = NULL, replace = TRUE) {
  if(!is.null(batch_size)) batch_size <- as.integer(batch_size)
  function(N) {
    if(is.null(batch_size)) batch_size <<- N
    if(batch_size > N & replace) {
      sample(N, replace = TRUE)
    } else {
      sample(N, batch_size, replace)
    }
  }
}

batch_random_chunk <- function(size, 
                               chunks = 1L, 
                               replace = TRUE, 
                               shuffle = FALSE) {
  shuffle_index <- NULL
  force(shuffle)
  function(N) {
    if(shuffle & is.null(shuffle_index)) {
      rand_indicies <- sample(N)
    }
    max_chunk <- max(1, floor(N / size) - 1)
    if(is.infinite(chunks)) chunks <- max_chunk
    if(chunks > max_chunk & replace) {
      chunk_id <- sample(1:max_chunk, replace = TRUE)
    } else {
      chunk_id <- sample(1:max_chunk, chunks, replace)
    }
    index <- numeric(chunks * size)
    if(is.null(shuffle_index)) {
        for(i in seq_along(chunk_id)) {
          index[((i - 1) * size + 1):(i * size)] <-
          ((chunk_id[i] - 1) * size + 1):(chunk_id[i] * size)
        }
    } else {
        for(i in seq_along(chunk_id)) {
          index[((i - 1) * size + 1):(i * size)] <-
            shuffle_index[((chunk_id[i] - 1) * size + 1):(chunk_id[i] * size)]
        }
    }
    index
  }
}

batch_linear_chunk <- function(size, chunks = 1L, shuffle = FALSE) {
  chunks <- as.integer(chunks)
  chunk <- 1
  shuffle_index <- NULL
  force(shuffle)
  function(N) {
    if(shuffle & is.null(shuffle_index)) {
      rand_indicies <- sample(N)
    }
    max_chunk <- max(1, floor(N / size) - 1)
    if(chunk + chunks >= max_chunk) chunk <- 1
    chunk_id <- chunk:(chunk + chunks)
    index <- numeric(chunks * size)
    if(is.null(shuffle_index)) {
      for(i in seq_along(chunk_id)) {
        index[((i - 1) * size + 1):(i * size)] <-
          ((chunk_id[i] - 1) * size + 1):(chunk_id[i] * size)
      }
      chunk <- chunk + chunks
    } else {
      for(i in seq_along(chunk_id)) {
        index[((i - 1) * size + 1):(i * size)] <-
          shuffle_index[((chunk_id[i] - 1) * size + 1):(chunk_id[i] * size)]
      }
      chunk <- chunk + chunks
    }
    index
  }
}