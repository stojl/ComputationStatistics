mini_batch <- function(batch_size = 1) {
  function(par0, x, y, loss_gr, gamma) {
    M <- floor(length(x) / batch_size)
    par <- par0
    for(i in 1:M) {
      batch_indicies <- ((i - 1) * batch_size + 1):(i * batch_size)
      gr <- loss_gr(par, x[batch_indicies], y[batch_indicies])
      par <- par - gamma(i) * gr
    }
    par
  }
}
