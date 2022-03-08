gen_data <- function(n, sig2 = 5, scenario = c("base", "exchange", "ps-mis", "out-mis",
                                               "sample-overlap", "treat-overlap")){
  
  # covariates
  x1 <- stats::rnorm(n, 0, 1)
  x2 <- stats::rnorm(n, 0, 1)
  x3 <- stats::rnorm(n, 0, 1)
  x4 <- stats::rnorm(n, 0, 1)
  
  # transformed predictors
  u1 <- as.numeric(scale(exp((x1 + x4)/2)))
  u2 <- as.numeric(scale(x2/(1 + exp(x1))))
  u3 <- as.numeric(scale(log(abs(x2*x3))))
  u4 <- as.numeric(scale((x3 + x4)^2))
  
  # create matrix
  X <- cbind(int = rep(1, n), x1, x2, x3, x4)
  U <- cbind(int = rep(1, n), u1, u2, u3, u4)
  # V <- cbind(int = rep(1, n), v1, v2, v3, v4)
  
  # coefficients
  beta0 <- c(2, -3, -1, 1, 3)
  beta1 <- c(0, 2, -2, -2, 2)
  alpha <- c(-2, -1, 3, -3, 1)
  lambda <- c(0, 1, -1, 1, -1)
  gamma <- c(0.5, -0.5, 0.5, -0.5, 0.5)
  
  # beta0 <- c(4, -3, -1, 1, 3)
  # beta1 <- c(2, -1, -3, 3, 1)
  # alpha <- c(2, 2, 2, -2, -2)
  # lambda <- c(0.5, -0.5, 0.5, -0.5, 0.5)
  # delta <- c(-0.25, 0, 0, 0, 0)
  # gamma <- c(-0.25, 0.25, -0.25, 0.25, -0.25)
  
  if (scenario == "base"){
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, 2) 
    mu_0 <- c(X %*% beta0)
    mu_1 <- mu_0 + a*c(X %*% alpha)
    PATE <- mean(X[s == 0,] %*% alpha)
  } else if (scenario == "exchange"){
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, 2) 
    mu_0 <- c(s*(X %*% beta0) + (1 - s)*(X %*% beta1))
    mu_1 <- mu_0 + a*c(X %*% alpha)
    PATE <- mean(X[s == 0,] %*% alpha)
  } else if (scenario == "ps-mis") {
    f_X <- c(1/(1 + exp(-U %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(U %*% lambda)
    a <- rnorm(n, e_X, 2) 
    mu_0 <- c(s*(X %*% beta0) + (1 - s)*(X %*% beta1))
    mu_1 <- mu_0 + a*c(X %*% alpha)
    PATE <- mean(X[s == 0,] %*% alpha)
  } else if (scenario == "out-mis"){
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, 2) 
    mu_0 <- c(s*(U %*% beta0) + (1 - s)*(U %*% beta1))
    mu_1 <- mu_0 + a*c(U %*% alpha)
    PATE <- mean(U[s == 0,] %*% alpha)
  } else if (scenario == "sample-overlap") {
    f_X <- c(1/(1 + exp(-X %*% (4*gamma))))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, 2) 
    mu_0 <- c(U %*% beta0)
    mu_1 <- mu_0 + a*c(U %*% alpha)
    PATE <- mean(U[s == 0,] %*% alpha)
  } else if (scenario == "treat-overlap"){
    f_X <- c(1/(1 + exp(-X %*% gamma)))
    s <- rbinom(n, 1, f_X) # sample
    e_X <- c(X %*% lambda)
    a <- rnorm(n, e_X, 1) # treatment
    mu_0 <- c(U %*% beta0)
    mu_1 <- mu_0 + a*c(U %*% alpha)
    PATE <- mean(U[s == 0,] %*% alpha)
  }
  
  # observed outcome
  y <- rnorm(n, mu_0, sqrt(sig2))
  
  # create simulation dataset
  sim <- list(y = y, z = z, s = s, X = X, U = U, PATE = PATE)
  
  return(sim)
  
}
