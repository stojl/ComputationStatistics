stop_parm <- function(epsilon, k) {
  force(epsilon)
  force(k)
  n <- 0
  function(par0, par, loss_gr, loss) {
    if(norm(par0 - par, "2") < 
       epsilon * (norm(par0, "2") + epsilon)) {
      n <<- n + 1
    } else {
      n <<- 0
    }
    n == k
  }
}
