kern_est_dr <- function(a.new, psi, mhat.mat, a, s, a.vals, bw, 
                        se.fit = FALSE, weights = NULL, int.mat = NULL) {

  mhat <- mhat.mat[,which(a.vals == a.new)]
  
  n1 <- sum(s)
  n0 <- sum(1 - s)
  
  if(is.null(weights))
    weights <- rep(1, times = n0 + n1)
  
  a.std <- (a - a.new) / bw
  k.std.tmp <- dnorm(a.std) / bw
  k.std <- c(n1*k.std.tmp/sum(k.std.tmp), rep(1, n0))
  g.std <- as.matrix(Matrix::bdiag(cbind(1, a.std), rep(1, n0)))
  xi <- c(psi, mhat)
  
  b <- lm(xi ~ -1 + g.std, weights = weights*k.std)$coefficients
  mu <- unname(b[1] + b[3])
  
  if (se.fit) {
    
    eta <- c(g.std %*% b)
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- weights * k.std * g.std * (xi - eta)
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = Sig[1,1] + Sig[3,3]))
    
  } else
    return(mu)
  
}

kern_est_ipw <- function(a.new, y, a, a.vals, bw, se.fit = FALSE, weights = NULL, int.mat = NULL) {
  
  n1 <- length(y)
  
  if(is.null(weights))
    weights <- rep(1, times = n1)
  
  a.std <- (a - a.new) / bw
  k.std.tmp <- dnorm(a.std) / bw
  k.std <- n1*k.std.tmp/sum(k.std.tmp)
  g.std <- cbind(1, a.std)
  
  b <- lm(y ~ -1 + g.std, weights = weights*k.std)$coefficients
  mu <- unname(b[1])
  
  if (se.fit) {
    
    eta <- c(g.std %*% b)
    
    U <- solve(crossprod(g.std, weights*k.std*g.std))
    V <- weights * k.std * g.std * (xi - eta)
    Sig <- U %*% crossprod(V) %*% U
    
    return(c(mu = mu, sig2 = Sig[1,1]))
    
  } else
    return(mu)
  
}