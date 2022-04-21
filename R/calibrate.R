
### Entropy Balancing (Exact)

calibrate <- function(cmat, target, base_weights = NULL, coefs_init = NULL,
                      optim_ctrl = list(maxit = 500, reltol = 1e-10), ...) {
  
  if (!is.matrix(cmat))
    stop("cmat must be a matrix")
  
  if (!is.vector(target))
    stop("target must be a vector")
  
  fn <- match.fun(lagrange_ent)
  
  if (is.null(base_weights)) { # initialize base_weights
    base_weights <- rep(1, nrow(cmat))
  } else if (length(base_weights) != nrow(cmat)) { 
    stop("length(base_weights) != sample size")
  }
  
  # initialize coefs
  if (is.null(coefs_init)) {
    coefs_init <- rep(0, times = ncol(cmat)) 
  } else if (length(coefs_init) != ncol(cmat)) {
    stop("length(coefs_init) != ncol(cmat)")
  }
  
  extraArgs <- list(...)
  
  if (length(extraArgs)) {
    
    arg <- names(formals(stats::optim))
    indx <- match(names(extraArgs), arg, nomatch = 0)
    if (any(indx == 0)) 
      stop(paste("Argument", names(extraArgs)[indx == 0], "not matched"))
    
  }
  
  opt <- stats::optim(coefs_init, fn, method = "BFGS",
                      cmat = cmat,
                      base_weights = base_weights,
                      target = target,
                      control = optim_ctrl, ...)
  
  converged <- ifelse(opt$convergence == 0, TRUE, FALSE)
  coefs <- opt$par
  weights <- c( base_weights*exp(-cmat %*% coefs) )
  
  out <- list(weights = weights,
              coefs = coefs,
              converged = converged,
              cmat = cmat,
              target = target,
              base_weights = base_weights, 
              optim_ctrl = optim_ctrl)
  
  class(out) <- "cfit"
  return(out)
  
}

lagrange_ent <- function(coefs, cmat, target, base_weights) {
  
  temp <- sum(base_weights*exp(-cmat %*% coefs))
  out <- temp + sum(target * coefs)
  return(out)
  
}

### Quadratic Programming (Approximate):


# standardize <- function(X, target, lambda = NULL, lam_seq = seq(0, 0.25, length.out = 25)) {
#   
#   P <- diag(1, nrow(X))
#   q <- rep(0, nrow(X))
#   A <- cbind(rep(1, nrow(X)), X, -X, diag(1, nrow(X)))
#   
#   if (is.null(lambda) & !is.null(lam_seq)) {
#     
#     gammas <- sapply(lam_seq, function(lambda, ...){
#       
#       l <- 1e-4/nrow(X)
#       b <- c(1, target - lambda, -target - lambda, rep(l, nrow(X)))
#       
#       balance.soln <- try(quadprog::solve.QP(P, q, A, b, meq = 1), silent = TRUE)
#       
#       if (inherits(balance.soln, "try-error"))
#         return(rep(1/nrow(X), nrow(X)))
#       else
#         return(balance.soln$solution)
#       
#     })
#     
#     boot_idx <- replicate(100, sample(1:nrow(X), nrow(X), replace = TRUE))
#     
#     bal <- sapply(1:ncol(gammas), function(j, ...){
#       
#       gamma.tmp <- gammas[,j]
#       
#       imbalance <- apply(boot_idx, 2, function(k, ...)
#         sqrt(sum(c(target - c(t(X[k,]) %*% gamma.tmp[k])/sum(gamma.tmp[k]))^2)))
#       
#       mean(imbalance)
#       
#     })
#     
#     gamma <- gammas[,which.min(bal)]
#     
#   } else if (!is.null(lambda)) {
#     
#     l <- 1e-4/nrow(X)
#     b <- c(1, target - lambda, -target - lambda, rep(l, nrow(X)))
#     balance.soln <- quadprog::solve.QP(P, q, A, b, meq = 1)
#     gamma <- balance.soln$solution
#     
#   }
#   
#   imbalance <- abs(c(target - t(X) %*% gamma))
#   names(imbalance) <- colnames(X)
#   
#   return(list(weights = gamma, imbalance = imbalance))
#   
# }

### Based on https://github.com/ebenmichael/balancer

standardize <- function(X, target, lambda = 1, lowlim = 1e-6, uplim = 1,
                        exact = FALSE, eps_abs = 1e-6, eps_rel = 1e-6, ...) {
  
  # convert X to a matrix
  X <- as.matrix(X)
  
  # ensure that target is a vector
  target <- c(target)
  
  # dimension of auxiliary weights
  m <- ncol(X)
  n <- nrow(X)
  
  # creates P and q
  P <- Matrix::bdiag(Matrix::Diagonal(n, 0), Matrix::Diagonal(m, 1))
  q <- Matrix::sparseVector(-c(X %*% target), 1:n, n + m)
  I0 <- Matrix::bdiag(Matrix::Diagonal(n, 1), Matrix::Diagonal(m, 0))
  P <- P + lambda * I0
  
  # creates A, l, and u
  constraints <- create_constraints(X, target, lowlim, uplim, exact)
  
  settings <- do.call(osqp::osqpSettings,
                      c(list(eps_rel = eps_rel,
                             eps_abs = eps_abs,
                             verbose = FALSE),
                        list(...)))
  
  solution <- osqp::solve_osqp(P, q, constraints$A, constraints$l, 
                               constraints$u, pars = settings)
  
  weights <- solution$x[1:n]
  
  # compute imbalance matrix
  imbalance <- abs(c(target - t(X) %*% weights))
  names(imbalance) <- colnames(X)
  
  weights[weights < 0] <- 0
  weights <- n*weights
  
  return(list(weights = weights, imbalance = imbalance))
  
}

create_constraints <- function(X, target, lowlim, uplim, exact) {
  
  n <- nrow(X)
  m <- ncol(X)
  Xt <- t(X)
  
  # sum-to-one constraint for each group
  A1 <- Matrix::Matrix(1, nrow = 1, ncol = n)
  A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow = nrow(A1), ncol = m))
  l1 <- 1
  u1 <- 1
  
  # upper and lower bounds
  A2 <- Matrix::Diagonal(n)
  A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = m))
  l2 <- rep(lowlim, n)
  u2 <- rep(uplim, n)
  
  if(exact) {
    # Constrain the overall mean to be equal to the target
    A3 <- Matrix::Matrix(Xt)
    A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = m))
    l3 <- target
    u3 <- target
  } else {
    # skip this constraint and just make empty
    A3 <- NULL
    l3 <- NULL
    u3 <- NULL
  }
  
  # constrain the auxiliary weights to be sqrt(P)'gamma
  sqrtP <- Matrix::Matrix(Xt)
  A4 <- Matrix::cbind2(sqrtP, -Matrix::Diagonal(m))
  l4 <- rep(0, m)
  u4 <- rep(0, m)
  
  A <- rbind(A1, A2, A3, A4)
  l <- c(l1, l2, l3, l4)
  u <- c(u1, u2, u3, u4)
  
  return(list(A = A, l = l, u = u))
  
}

