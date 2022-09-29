rng_stream <- function(batch_size = 10000L, sample_function, ...) {
  
  pool <- sample_function(batch_size, ...)
  pool_size <- batch_size
  pool_index <- 1
  
  function(n) {
    if(n <= pool_size) {
      pool_request <- pool[pool_index:(pool_index + n - 1)]
      pool_index <<- pool_index + n
      pool_size <<- pool_size - n
      
      return(pool_request)
    } else {

      new_batches <- ceiling(n / batch_size)
      
      pool <<- sample_function(batch_size * new_batches, ...)
      
      pool_request <- pool[1:n]
      pool_index <<- n - (new_batches - 1) + 1
      pool_size <<- batch_size * new_batches - n
      
      return(pool_request)
    }
  }
}
