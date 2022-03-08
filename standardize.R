standardize_quadprog <- function(X, target, lambda = NULL, lam_seq = seq(0, sqrt(nrow(X))/nrow(X), length.out = 100), neg_weights = FALSE) {
  
  # The system is effectively
  # minimize zeta * delta^2 + (1 - zeta) * ||gamma||^2
  # subject to
  #   sum gamma = 1
  #   -target_j - lambda_j <= - (X'gamma)_j
  #   target_j - lambda_j <= (X'gamma)_j
  #   gamma_i > 0
  
  Dmat = diag(1, nrow(X))
  dvec = rep(0, nrow(X))
  Amat = cbind(rep(1, nrow(X)), X, -X)
  
  if (!neg_weights)
    Amat = cbind(Amat, diag(rep(1, nrow(X))))
  
  if (is.null(lambda) & !is.null(lam_seq)) {
    
    gammas <- sapply(lam_seq, function(lambda, ...){
      
      bvec = c(1, target, -target) - c(0, rep(lambda, 2*ncol(X)))
      
      if (!neg_weights) {
        l = 1e-4/nrow(X)
        bvec = c(bvec, rep(l, nrow(X)))
      }
      
      balance.soln <- try(quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1), silent = TRUE)
      
      if (inherits(balance.soln, "try-error"))
        return(rep(1/nrow(X), nrow(X)))
      else
        return(balance.soln$solution)
      
    })
    
    idx <- sample(1:nrow(X), nrow(X), replace = TRUE)
    
    bal <- sapply(1:ncol(gammas), function(j, ...){
      
      gamma <- gammas[idx,j]
      X.tmp <- X[idx,]
      imbalance <- sqrt(sum(c(target - c(t(X.tmp) %*% gamma)/sum(gamma))^2))
      
    })
    
    gamma <- gammas[,which.min(bal)]
    
  } else if (!is.null(lambda)) {
    
    bvec = c(1, target, -target) - c(0, rep(lambda, 2*ncol(X)))
    
    if (!neg_weights) {
      l = 1e-4/nrow(X)
      bvec = c(bvec, rep(l, nrow(X)))
    }
    
    balance.soln <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
    gamma <- balance.soln$solution
    
  }
  
  imbalance <- abs(c(target - t(X) %*% gamma))
  names(imbalance) <- colnames(X)
  
  return(list(weights = gamma, imbalance = imbalance))
  
}
