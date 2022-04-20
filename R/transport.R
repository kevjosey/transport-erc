glm_transport <- function(a, y, x1, x0, offset, family, df = 4,
                          a.vals = seq(min(a), max(a), length.out = 100)) {
  
  colnames(x0) <- colnames(x1)
  target <- colMeans(x0)
  n1 <- nrow(x1)

  # estimate nuisance outcome model with splines
  mumod <- glm(y ~ ns(a, df) + . + a:., data = data.frame(x1[,-1]), offset = offset, family = family)

  estimate <- sapply(a.vals, function(a.tmp, ...) {
    xa.tmp <- data.frame(x0, a = a.tmp)
    return(mean(predict(mumod, newdata = xa.tmp, type = "response")))
  })
  
  return(estimate)
  
}

ipw_transport <- function(a, y, x1, x0, offset, family, df = 4, 
                          a.vals = seq(min(a), max(a), length.out = 100),
                          sl.lib = c("SL.mean", "SL.glm", "SL.earth", "SL.glmnet", "SL.ranger")) {
  
  colnames(x0) <- colnames(x1)
  n1 <- nrow(x1)
  m <- ncol(x1)
  
  # estimate sampling weights
  xstar <- apply(x1, 2, function(x) x - mean(x))[,-1]
  astar <- c(a - mean(a))/var(a)
  astar2 <- -1 + c(a - mean(a))^2/var(a)
  mod <- calibrate(cmat = cbind(x1, x1*astar, x1*astar2), 
                      target = c(colMeans(x0), rep(0, m), rep(0, m)))
  wts <- mod$weights
  new_mod <- glm(y ~ ns(a, df = df), offset = offset, family = family, weights = wts)
  estimate <- predict(new_mod, newdata = data.frame(a = a.vals), type = "response")
  
  return(estimate)
  
}

dr_transport <- function(a, y, x1, x0, offset, family, df = 4,
                           a.vals = seq(min(a), max(a), length.out = 100)){
  
  colnames(x0) <- colnames(x1)
  n1 <- nrow(x1)
  m <- ncol(x1)
  
  # estimate nuisance outcome model with splines
  mumod <- glm(y ~ ns(a, df) + . + a:., data = data.frame(x1[,-1]), offset = offset, family = family)
  muhat <- predict(mumod, newdata = data.frame(x1[,-1], a = a), type = "response")
  
  # outcomes models given a.vals
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    xa.tmp <- data.frame(x0[,-1], a = a.tmp)
    return(predict(mumod, newdata = xa.tmp, type = "response"))
  })
  
  mhat <- predict(smooth.spline(a.vals, colMeans(muhat.mat)), x = a)$y
  
  # estimate sampling weights
  xstar <- apply(x1, 2, function(x) x - mean(x))[,-1]
  astar <- c(a - mean(a))/var(a)
  astar2 <- -1 + c(a - mean(a))^2/var(a)
  mod <- calibrate(cmat = cbind(x1, x1*astar, x1*astar2), 
                   target = c(colMeans(x0), rep(0, m), rep(0, m)))
  wts <- mod$weights
  psi <- c((y - muhat)*wts + mhat)
  
  new_mod <- glm(psi ~ ns(a, df = df), offset = offset, family = family)
  estimate <- predict(new_mod, newdata = data.frame(a = a.vals), type = "response")
  
  return(estimate)
  
}

tmle_transport <- function(a, y, x1, x0, offset, family, df = 4,
                           a.vals = seq(min(a), max(a), length.out = 100),
                           sl.lib = c("SL.mean", "SL.glm", "SL.earth", "SL.glmnet", "SL.ranger")){
  
  colnames(x0) <- colnames(x1)
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  
  # estimate nuisance outcome model with splines
  mumod <- glm(y ~ ns(a, df) + . + a:., data = data.frame(x1[,-1]), offset = offset, family = family)
  muhat <- predict(mumod, newdata = data.frame(x1[,-1], a = a), type = "response")
  
  # outcomes models given a.vals
  muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    xa.tmp <- data.frame(x0[,-1], a = a.tmp)
    return(predict(mumod, newdata = xa.tmp, type = "response"))
  })
  
  # estimate nuisance GPS parameters with lm
  pimod <- SuperLearner(Y = a, X = data.frame(x1[,-1]), SL.library = sl.lib, family = gaussian())
  pimod.vals1 <- c(pimod$SL.predict)
  pimod.vals0 <- predict(pimod, newdata = data.frame(x0[,-1]))$pred
  pi2mod.vals <- mean((a - pimod.vals1)^2)
  
  # estimate sampling weights
  rhomod <- calibrate(cmat = x1, target = colSums(x0))
  rhohat <- rhomod$weights
  
  # parametric density
  pihat <- dnorm(a, pimod.vals1, sqrt(pi2mod.vals))
  pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
    dnorm(a.tmp, c(pimod.vals1, pimod.vals0), sqrt(pi2mod.vals))
  })
  phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n1,], na.rm = T)), x = a)$y
  phat[phat<0] <- 1e-4
  
  # nonparametric density
  # a.std <- c(a - pimod.vals1) / sqrt(pi2mod.vals)
  # dens <- density(a.std)
  # pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / sqrt(pi2mod.vals)
  # 
  # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
  #   std <- c(a.tmp - c(pimod.vals1, pimod.vals0)) / sqrt(pi2mod.vals)
  #   approx(x = dens$x, y = dens$y, xout = std)$y / sqrt(pi2mod.vals)
  # })
  # 
  # phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n1,], na.rm = T)), x = a)$y
  # phat[phat<0] <- 1e-4  
  
  # TMLE update
  nsa <- ns(a, df = df + 1, intercept = TRUE)
  ipw <- phat/pihat
  trunc0 <- quantile(ipw, 0.01)
  trunc1 <- quantile(ipw, 0.99)
  ipw[ipw < trunc0] <- trunc0
  ipw[ipw > trunc1] <- trunc1
  base <- nsa*ipw*rhohat
  new_mod <- glm(y ~ 0 + base, offset = family$linkfun(muhat) + offset, family = family)
  param <- coef(new_mod)
  
  # predict spline basis and impute
  estimate <- sapply(1:length(a.vals), function(k, ...) {
    
    muhat.tmp <- muhat.mat[,k]
    pihat.tmp <- pihat.mat[,k]
    a.tmp <- a.vals[k]
    wts <- c(mean(pihat.tmp[-(1:n1)], na.rm = TRUE)/pihat.tmp[-(1:n1)])
    wts[wts < trunc0] <- trunc0
    wts[wts > trunc1] <- trunc1
    wts2 <- wts*exp(-c(x0%*%rhomod$coefs))
    mat <- predict(nsa, newx = rep(a.tmp, n0))*wts2
    return(mean(family$linkinv(family$linkfun(muhat.tmp) + c(mat%*%param)), na.rm = TRUE))
    
  })
  
  return(list(estimate = estimate, ipw = ipw, ios = rhohat))
  
}
