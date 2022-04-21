kern_est <- function(a.new, a, psi, bw, se.fit = FALSE, weights = NULL, int.mat = NULL, a.vals = NULL) {
    
  if(is.null(weights))
    weights <- rep(1, times = length(y))
  
  n <- length(a)
  
  ## LOESS Kernel
  
  # # subset index
  # a.std <- a - a.new
  # k <- floor(min(bw, 1)*length(a))
  # idx <- order(abs(a.std))[1:k]
  # 
  # # subset
  # a.std <- a.std[idx]
  # psi <- psi[idx]
  # max.a.std <- max(abs(a.std))
  # 
  # # construct kernel weight
  # k.std <- c((1 - abs(a.std/max.a.std)^3)^3)
  # g.std <- cbind(1, a.std)
  
  ## Gaussian Kernel
  a.std <- (a - a.new) / bw
  k.std <- dnorm(a.std) / bw
  g.std <- cbind(1, a.std)
  
  b <- lm(psi ~ -1 + g.std, weights = weights*k.std)$coefficients
  mu <- b[1]
  
  if (se.fit & !is.null(int.mat)) {
    
    eta <- c(g.std %*% b)
    
    ## LOESS
    # kern.mat <- matrix(rep(c((1 - abs((a.vals - a.new)/max.a.std)^3)^3), k), byrow = T, nrow = k)
    # kern.mat[matrix(rep(abs(a.vals - a.new)/max.a.std, k), byrow = T, nrow = k) > 1] <- 0
    # g.vals <- matrix(rep(c(a.vals - a.new), k), byrow = T, nrow = k)
    # intfn1.mat <- kern.mat * int.mat[idx,]
    # intfn2.mat <- g.vals * kern.mat * int.mat[idx,]
    # 
    # int1 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
    #                 (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2, 1, sum)
    # int2 <- apply(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), k), byrow = T, nrow = k)*
    #                 (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2, 1, sum)
    
    ## Gaussian
    
    int1 <- colMeans(k.std * t(int.mat))
    int2 <- colMeans(g.std * k.std * t(int.mat))
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- cbind(weights*(k.std * (psi - eta) + int1),
               weights*(a.std * k.std * (psi - eta) + int2))
    sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig = sig[1,1]))
    
  } else
    return(mu)
  
}