gen_data <- function(n0, n1, sig2 = 5, tau2 = 1, scenario, a.vals) {
  
  n <- n1 + n0
  
  # covariates
  x01 <- stats::rnorm(n0, -1, 4)
  x02 <- stats::rnorm(n0, 1, 4)
  x03 <- stats::rnorm(n0, 0, 2)
  x04 <- stats::rbinom(n0, 1, 0.3)
  
  x11 <- stats::rnorm(n1, 1, 4)
  x12 <- stats::rnorm(n1, -1, 4)
  x13 <- stats::rnorm(n1, 0, 1)
  x14 <- stats::rbinom(n1, 1, 0.7)
  
  x1 <- c(x01, x11)
  x2 <- c(x02, x12)
  x3 <- c(x03, x13)
  x4 <- c(x04, x14)
  
  u1 <- as.numeric(scale(exp(x1/2)))
  u2 <- as.numeric(scale(x2/(1 + exp(x1)) + 10))
  u3 <- as.numeric(scale((x1*x3/25 + 0.6)^3))
  u4 <- as.numeric(scale((x2 + x4 + 10)^2))
  
  # create matrix
  X <- cbind(int = rep(1, n), x1, x2, x3, x4)
  U <- cbind(int = rep(1, n), u1, u2, u3, u4)
  
  s <- rep(c(0,1), c(n0, n1))
  
  # coefficients
  lambda <- c(10, 1, -1, -1, 1)
  beta <- c(2,-1,3,-3,1)
  
  # beta0 <- c(4, -3, -1, 1, 3)
  # beta1 <- c(2, -1, -3, 3, 1)
  # alpha <- c(2, 2, 2, -2, -2)
  # lambda <- c(0.5, -0.5, 0.5, -0.5, 0.5)
  # delta <- c(-0.25, 0, 0, 0, 0)
  # gamma <- c(-0.25, 0.25, -0.25, 0.25, -0.25)
  
  if (scenario == "base"){
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, sqrt(tau2)) # treatment
    mu <- c(X %*% beta) + (a - 10) - 
      2*cos(pi*(a - 6)/4) - (a - 10)*x1
  } else if (scenario == "ps-mis"){
    e_X <- c(U %*% lambda)
    a <- rnorm(n, e_X, sqrt(tau2))
    mu <- c(X %*% beta) + (a - 10) - 
      2*cos(pi*(a - 6)/4) - (a - 10)*x1
  } else if (scenario == "out-mis") {
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, sqrt(tau2))
    mu <- c(U %*% beta) + (a - 10) - 
      2*cos(pi*(a - 6)/4) - (a - 10)*u1
  } else if (scenario == "ps-overlap"){
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, (1/2)*sqrt(tau2))
    mu <- c(U %*% beta) + (a - 10) - 
      2*cos(pi*(a - 6)/4) - (a - 10)*u1
  } else if (scenario == "samp-overlap"){
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, sqrt(tau2))
    mu <- c(U %*% beta) + (a - 10) - 
      2*cos(pi*(a - 6)/4) - (a - 10)*u1
  }
  
  ERC <- rep(NA, length(a.vals))
  
  for(i in 1:length(a.vals)) {
    
    a.vec <- rep(a.vals[i], sum(s == 0))
    
    if (scenario %in% c("out-mis", "ps-overlap", "samp-overlap")) {
      mu_out <- U[s == 0,]%*% beta + 0.25*(a.vec - 10) - 0.75*cos(pi*(a.vec - 6)/4) - 0.25*(a.vec - 10)*X[s == 0,2]
    } else { # out_scen == "a"
      mu_out <- X[s == 0,] %*% beta + 0.25*(a.vec - 10) - 0.75*cos(pi*(a.vec - 6)/4) - 0.25*(a.vec - 10)*U[s == 0,2]
    }
    
    ERC[i] <- mean(mu_out)
    
  }
  
  # potential outcomes
  y <- rnorm(n, mu, sqrt(sig2))
  
  # create simulation dataset
  sim <- list(y = y, a = a, s = s, X = X, U = U, ERC = ERC)
  
  return(sim)
  
}
