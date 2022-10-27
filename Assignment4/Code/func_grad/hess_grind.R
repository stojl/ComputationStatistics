get_hessian_grad <- function(x, y) {
  n <- length(x)
  force(y)
  
  function(par) {
    alpha <- par[1]
    beta <- par[2]
    gamma <- par[3]
    rho <- par[4]

    elogx <- exp(beta * log(x) - alpha)
    
    da <- (rho - gamma) * elogx / (1 + elogx)^2
    db <- -log(x) * da
    dr <- 1 / (1 + elogx)
    dg <- 1 - dr
    
    yf <- y - gamma - (rho - gamma) * dr
    
    gg <- 0
    rr <- 0
    gr <- 0
    aa <- 2 * da^2 / (1 + elogx) - da
    bb <- aa
    ab <- -log(x) * aa
    ar <- elogx * dr^2
    ag <- 1- ar
    br <- -log(x) * ar
    bg <- 1 - br
    
    Gaa <- mean(da^2)
    Gbb <- mean(db^2)
    Ggg <- mean(dg^2)
    Grr <- mean(dr^2)
    Gab <- mean(da * db)
    Gag <- mean(da * dg)
    Gar <- mean(da * dr)
    Gbg <- mean(db * dg)
    Gbr <- mean(db * dr)
    Ggr <- mean(dg * dr)
    
    Haa <- mean(yf * aa)
    Hbb <- mean(yf * bb)
    Hgg <- mean(yf * gg)
    Hrr <- mean(yf * rr)
    Hab <- mean(yf * ab)
    Hag <- mean(yf * ag)
    Har <- mean(yf * ar)
    Hbg <- mean(yf * bg)
    Hbr <- mean(yf * br)
    Hgr <- mean(yf * gr)
    
    AA <- Gaa - Haa
    AB <- Gab - Hab
    AG <- Gag - Hag
    AR <- Gar - Har
    BB <- Gbb - Hbb
    BG <- Gbg - Hbg
    BR <- Gbr - Hbr
    GG <- Ggg - Hgg
    GR <- Ggr - Hgr
    RR <- Grr - Hrr
    
    
    H <- c(AA, AB, AG, AR, 
           AB, BB, BG, BR,
           AG, BG, GG, GR,
           AR, BR, GR, RR) / 2
    
    dim(H) <- c(4, 4)
    grad <- c(-mean(yf * da, na.rm = TRUE) / 2,
              -mean(yf * db, na.rm = TRUE) / 2,
              -mean(yf * dg, na.rm = TRUE) / 2,
              -mean(yf * dr, na.rm = TRUE) / 2)
    
    list(H = H, gr = grad)
  }
}
