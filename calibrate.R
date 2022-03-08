
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

### Quadratic Programming (Approximate) Based on:
### https://github.com/ebenmichael/balancer

standardize <- function(X, target, lambda = 0, lowlim = 0, uplim = 1,
                        return_data = T, exact_global = F, scale_sample_size = F,
                        eps_abs = 1e-5, eps_rel = 1e-5, ...) {
  
  # convert X to a matrix
  X <- as.matrix(X)
  
  # ensure that target is a vector
  target <- c(target)
  
  # dimension of auxiliary weights
  m <- ncol(X)
  n <- nrow(X)
  
  q <- create_q_vector(X, target, m)
  P <- create_P_matrix(n, m)
  
  I0 <- create_I0_matrix(X, n, m, scale_sample_size)
  P <- P + lambda * I0
  
  constraints <- create_constraints(X, target, lowlim, uplim, exact_global)
  
  settings <- do.call(osqp::osqpSettings,
                      c(list(eps_rel = eps_rel,
                             eps_abs = eps_abs),
                        list(...)))
  
  solution <- osqp::solve_osqp(P, q, constraints$A,
                               constraints$l, constraints$u,
                               pars = settings)
  
  weights <- solution$x[1:nrow(X)]
  
  # compute imbalance matrix
  imbalance <- as.matrix(target - t(X) %*% weights)
  
  if(return_data) {
    data_out <- list(P = P  - lambda * I0, q = q, constraints = constraints)
  } else {
    data_out <- NULL
  }
  
  return(list(weights = weights, imbalance = imbalance, data_out = data_out))
  
}

create_I0_matrix <- function(X, n, m, scale_sample_size = FALSE) {
  
  if(scale_sample_size) {
    I0 <- Matrix::Diagonal(n, n)
  } else {
    I0 <- Matrix::Diagonal(n)
  }
  I0 <- Matrix::bdiag(I0, Matrix::Diagonal(m, 0))
  return(I0)
  
}

create_q_vector <- function(X, target, m) {
  
  q <- -c(X %*% target)
  q <- Matrix::sparseVector(q, 1:length(q), length(q) + m)
  return(q)
  
}

create_P_matrix <- function(n, m) {
  
  return(Matrix::bdiag(Matrix::Diagonal(n, 0), Matrix::Diagonal(m, 1)))
  
}

create_constraints <- function(X, target, lowlim, uplim, exact_global) {
  
  n <- nrow(X)
  m <- ncol(X)
  Xt <- t(X)
  
  # sum-to-one constraint for each group
  A1 <- Matrix::Matrix(t(matrix(1, nrow = n, ncol = 1)))
  A1 <- Matrix::cbind2(A1, Matrix::Matrix(0, nrow=nrow(A1), ncol = m))
  l1 <- 1
  u1 <- 1
  
  # upper and lower bounds
  A2 <- Matrix::Diagonal(n)
  A2 <- Matrix::cbind2(A2, Matrix::Matrix(0, nrow = nrow(A2), ncol = m))
  l2 <- rep(lowlim, n)
  u2 <- rep(uplim, n)
  
  if(exact_global) {
    # Constrain the overall mean to be equal to the target
    A3 <- Xt * n
    A3 <- Matrix::cbind2(A3, Matrix::Matrix(0, nrow = nrow(A3), ncol = m))
    l3 <- n * target
    u3 <- n * target
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

