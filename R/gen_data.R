gen_data <- function(n0, n1, sig2 = 2, tau2 = 1, scenario, a.vals) {
  
  n <- n1 + n0
  
  # covariates
  x01 <- stats::rnorm(n0, 0, 1)
  x02 <- stats::rnorm(n0, -1, 2)
  x03 <- stats::rnorm(n0, 1, 2)
  x04 <- stats::rnorm(n0, -1, 2)
  
  x11 <- stats::rnorm(n1, 0, 1)
  x12 <- stats::rnorm(n1, 1, 2)
  x13 <- stats::rnorm(n1, -1, 2)
  x14 <- stats::rnorm(n1, 1, 2)
  
  x1 <- c(x01, x11)
  x2 <- c(x02, x12)
  x3 <- c(x03, x13)
  x4 <- c(x04, x14)
  
  # u11 <- as.numeric(scale(exp(x11/2)))
  # u12 <- as.numeric(scale(x12/(1 + exp(x11)) + 10))
  # u13 <- as.numeric(scale((x11*x13/25 + 0.6)^3))
  # u14 <- as.numeric(scale((x12 + x14 + 10)^2))
  # 
  # u01 <- as.numeric(scale(exp(x01/2)))
  # u02 <- as.numeric(scale(x02/(1 + exp(x01)) + 10))
  # u03 <- as.numeric(scale((x01*x03/25 + 0.6)^3))
  # u04 <- as.numeric(scale((x02 + x04 + 10)^2))
  # 
  # u1 <- c(u01, u11)
  # u2 <- c(u02, u12)
  # u3 <- c(u03, u13)
  # u4 <- c(u04, u14)
  
  u1 <- exp(x1/4)
  u2 <- x2/(1 + exp(x1)) 
  u3 <- abs(x1*x3) - 1
  u4 <- ((x2 + x4)/5)^2 - 1
  
  # create matrix
  X <- cbind(int = rep(1, n), x1, x2, x3, x4)
  U <- cbind(int = rep(1, n), u1, u2, u3, u4)
  
  s <- rep(c(0,1), c(n0, n1))
  
  # coefficients
  lambda <- c(10, 0.5, -0.5, -0.5, 0.5)
  beta <- c(1, -0.25, 0.75, -0.75, 0.25)
  alpha <- c(1, -0.75, -0.25, 0.25, 0.75)
  
  if (scenario == "base"){
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, sqrt(tau2)) # treatment
    mu <- c(X %*% beta) - 2*cos(pi*(a - 6)/4) +
      (a - 10)*(X %*% alpha)
  } else if (scenario == "ps-mis"){
    e_X <- c(U %*% lambda)
    a <- rnorm(n, e_X, sqrt(tau2))
    mu <- c(X %*% beta) - 2*cos(pi*(a - 6)/4) +
      (a - 10)*(X %*% alpha)
  } else if (scenario == "out-mis") {
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, sqrt(tau2))
    mu <- c(U %*% beta) - 2*cos(pi*(a - 6)/4) +
      (a - 10)*(U%*%alpha)
  } else if (scenario == "mis"){
    e_X <- c(U %*% lambda)
    a <- rnorm(n, e_X, sqrt(tau2))
    mu <- c(U %*% beta) - 2*cos(pi*(a - 6)/4) +
      (a - 10)*(U%*%alpha)
  }
    
  if (scenario == "out-mis" | scenario == "mis") {
    ERC <- sum(colMeans(U[s == 0,]) * beta) -
      2*cos(pi*(a.vals - 6)/4) +
      (a.vals - 10)*mean(U[s == 0,] %*% alpha)
  } else { # out_scen == "a"
    ERC <- sum(c(1, 0, -1, 1, -1) * beta) -
      2*cos(pi*(a.vals - 6)/4) +
      (a.vals - 10)*mean(X[s == 0,] %*% alpha)
  }
    
  # potential outcomes
  y <- rnorm(n, mu, sqrt(sig2))
  
  # create simulation dataset
  sim <- list(y = y, a = a, s = s, X = X, U = U, ERC = ERC)
  
  return(sim)
  
}
