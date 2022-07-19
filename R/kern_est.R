kern_est <- function(a.new, a, psi, bw, se.fit = FALSE, weights = NULL, int.mat = NULL, a.vals = NULL) {
    
  if(is.null(weights))
    weights <- rep(1, times = length(a))
  
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
    
    ## Gaussian
    kern.mat <- matrix(rep(dnorm((a.vals - a.new) / bw) / bw, n), byrow = T, nrow = n)
    g.vals <- matrix(rep(c(a.vals - a.new) / bw, n), byrow = T, nrow = n)
    intfn1.mat <- kern.mat * int.mat
    intfn2.mat <- g.vals * kern.mat * int.mat
    
    int1 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn1.mat[,-1] + intfn1.mat[,-length(a.vals)])/2)
    int2 <- rowSums(matrix(rep((a.vals[-1] - a.vals[-length(a.vals)]), n), byrow = T, nrow = n)*
                      (intfn2.mat[,-1] + intfn2.mat[,-length(a.vals)])/2)
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- cbind(weights*(k.std * (psi - eta) + int1),
               weights*(a.std * k.std * (psi - eta) + int2))
    sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig = sig[1,1]))
    
  } else
    return(mu)
  
}