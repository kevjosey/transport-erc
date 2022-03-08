glm_transport <- function(a1, y1, x1, x0, offset1, family, degree = 4,
                          a.vals = seq(min(a), max(a), length.out = 100)) {
  
  colnames(x0) <- colnames(x1)
  target <- colMeans(x0)
  n1 <- nrow(x1)

  # estimate nuisance outcome model with splines
  mumod <- glm(y ~ ns(a, df = degree) + ., data = x1, offset = offset1, family = family)

  estimate <- sapply(a.vals, function(a.tmp, ...) {
    xa.tmp <- data.frame(x0, a = a.tmp)
    colnames(xa.tmp) <- colnames(xa) 
    return(mean(predict(mumod, newdata = xa.tmp, type = "response")))
  })
  
  return(estimate)
  
}

ipw_transport <- function(a1, y1, x1, x0, offset1, family, degree = 4, 
                          a.vals = seq(min(a), max(a), length.out = 100),
                          sl.lib = c("SL.mean", "SL.glm", "SL.earth", "SL.glmnet", "SL.ranger")) {
  
  # estimate nuisance GPS parameters with lm
  pimod <- SuperLearner(Y = a1, X = x1, SL.library = sl.lib, family = gaussian())
  pimod.vals1 <- c(pimod$SL.predict)
  pimod.vals0 <- predict(pimod, newdata = data.frame(x0))$pred
  pi2mod.vals <- var(a1 - pimod.vals1)
  
  # estimate sampling weights
  rhomod <- standardize(X = X1, target = target, lambda = 1)
  rhohat1 <- rhomod$weights
  rhohat1[rhohat1 < 0] <- 0
  
  # parametric density
  # pihat1 <- dnorm(a1, pimod.vals1, sqrt(pi2mod.vals))
  # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  #   dnorm(a.tmp, c(pimod.vals1, pimod.vals0), sqrt(pi2mod.vals))
  # })
  # phat1 <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n1,], na.rm = T)), x = a1)$y
  # phat1[phat1<0] <- 1e-4
  
  # nonparametric density
  a.std1 <- c(a1 - pimod.vals1) / sqrt(pi2mod.vals)
  dens <- density(a.std1)
  pihat1 <- approx(x = dens$x, y = dens$y, xout = a.std1)$y / sqrt(pi2mod.vals)
  
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - c(pimod.vals1, pimod.vals0)) / sqrt(pi2mod.vals)
    approx(x = dens$x, y = dens$y, xout = std)$y / sqrt(pi2mod.vals)
  })
  
  phat1 <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n1,], na.rm = T)), x = a1)$y
  phat1[phat1<0] <- 1e-4  
  
  ipw1 <- phat1/pihat1
  wts <- ipw1*rhohat1
  new_mod <- glm(y1 ~ ns(a, df = degree), offset = offset1, family = family, weights = weights)
  estimate <- predict(new_mod, newdata = data.frame(a = a.vals), type = "response")
  
  return(estimate)
  
}

tmle_transport <- function(a1, y1, x1, x0, offset1, family, degree = 4,
                           a.vals = seq(min(a), max(a), length.out = 100),
                           sl.lib = c("SL.mean", "SL.glm", "SL.earth", "SL.glmnet", "SL.ranger")){
  
  colnames(x0) <- colnames(x1)
  target <- colMeans(x0)
  n1 <- nrow(x1)
  
  # estimate nuisance outcome model with splines
  mumod <- glm(y ~ ns(a, df = degree) + ., data = x1, offset = offset1, family = family)
  muhat1 <- predict(mumod, newdata = data.frame(x1, a = a1), type = "response")
  
  # outcomes models given a.vals
  muhat.mat0 <- sapply(a.vals, function(a.tmp, ...) {
    xa.tmp <- data.frame(x0, a = a.tmp)
    return(predict(mumod, newdata = xa.tmp, type = "response"))
  })
  
  # estimate nuisance GPS parameters with lm
  pimod <- SuperLearner(Y = a1, X = x1, SL.library = sl.lib, family = gaussian())
  pimod.vals1 <- c(pimod$SL.predict)
  pimod.vals0 <- predict(pimod, newdata = data.frame(x0))$pred
  pi2mod.vals <- var(a1 - pimod.vals1)
  
  # estimate sampling weights
  rhomod <- standardize(X = X1, target = target, lambda = 1)
  rhohat1 <- rhomod$weights*n1
  rhohat1[rhohat1 < 0] <- 1e-4
  
  # parametric density
  # pihat1 <- dnorm(a1, pimod.vals1, sqrt(pi2mod.vals))
  # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  #   dnorm(a.tmp, c(pimod.vals1, pimod.vals0), sqrt(pi2mod.vals))
  # })
  # phat1 <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n1,], na.rm = T)), x = a1)$y
  # phat1[phat1<0] <- 1e-4
  
  # nonparametric density
  a.std1 <- c(a1 - pimod.vals1) / sqrt(pi2mod.vals)
  dens <- density(a.std1)
  pihat1 <- approx(x = dens$x, y = dens$y, xout = a.std1)$y / sqrt(pi2mod.vals)
  
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    std <- c(a.tmp - c(pimod.vals1, pimod.vals0)) / sqrt(pi2mod.vals)
    approx(x = dens$x, y = dens$y, xout = std)$y / sqrt(pi2mod.vals)
  })
  
  phat1 <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n1,], na.rm = T)), x = a1)$y
  phat1[phat1<0] <- 1e-4  
  
  # TMLE update
  nsa <- ns(a, df = degree + 1, intercept = TRUE)
  ipw1 <- phat1/pihat1
  base <- nsa*ipw1*rhohat1
  new_mod <- glm(y1 ~ 0 + base, offset = log(muhat1) + offset1, family = family)
  param <- coef(new_mod)
  
  # predict spline basis and impute
  estimate <- sapply(1:length(a.vals), function(k, ...) {
    
    muhat.tmp <- muhat.mat0[,k]
    pihat.tmp <- pihat.mat[,k]
    a.tmp <- a.vals[k]
    wts <- c(mean(pihat.tmp[1:n1], na.rm = TRUE)/pihat.tmp[-(1:n1)])
    mat <- predict(nsa, newx = rep(a.tmp, n))*wts
    return(mean(exp(log(muhat.tmp) + c(mat%*%param)), na.rm = TRUE))
    
  })
  
  return(list(estimate = estimate, ipw = ipw1, ios = rhohat1))
  
}