glm_transport <- function(a, y, x1, x0, offset = NULL, weights = NULL, family = gaussian(),
                          df = 4, a.vals = seq(min(a), max(a), length.out = 100)) {
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
  if(is.null(weights))
    weights <- rep(1, times = length(y))
  
  colnames(x0) <- colnames(x1)
  target <- colMeans(x0)
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  
  mumod <- glm(y ~ a + cos(pi*(a - 6)/4) + x1 + x2 + x3 + x4 + a:(x1 + x2 + x3 + x4),
               data = data.frame(a = a, x1), weights = weights,
               offset = offset, family = family)
  muhat <- mumod$fitted.values
  
  out <- sapply(a.vals, function(a.tmp, ...) {
    xa.tmp <- data.frame(x0, a = a.tmp)
    return(mean(predict(mumod, newdata = xa.tmp, type = "response")))
  })
  
  return(out)
  
}

ipw_transport <- function(a, y, x1, x0, offset = NULL, weights = NULL, family = gaussian(),
                          df = 4, bw = 1, a.vals = seq(min(a), max(a), length.out = 100),
                          sl.lib = c("SL.mean", "SL.glm", "SL.earth", "SL.glmnet", "SL.ranger")) {
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
  if(is.null(weights))
    weights <- rep(1, times = length(y))
  
  x1.mat <- model.matrix(~ ., data = data.frame(x1))
  x0.mat <- model.matrix(~ ., data = data.frame(x0))
  
  colnames(x0) <- colnames(x1)
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  m <- ncol(x1)
  
  # estimate sampling weights
  mod <- calibrate(cmat = cbind(x1.mat), target = n1*colMeans(x0.mat))
  rhohat <- mod$weights
  pimod <- SuperLearner(Y = a, X = data.frame(x1), SL.lib = sl.lib, family = gaussian())
  pimod.vals <- c(pimod$SL.predict)
  
  sd.a <- sd(a - pimod.vals)
  a.std <- c(a - pimod.vals)/sd.a
  dens <- density(a.std)
  pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / sd.a
  
  phat <- sapply(a, function(a.tmp, ...) {
    a.std.tmp <- c(a.tmp - pimod.vals)/sd.a
    return(mean(approx(x = dens$x, y = dens$y, xout = a.std.tmp)$y / sd.a, na.rm = T))
  })
  
  wts <- rhohat*c(phat/pihat)
  
  # estimate sampling weights
  # astar <- c(a - mean(a))/var(a)
  # astar2 <- -1 + c(a - mean(a))^2/var(a)
  # mod <- calibrate(cmat = cbind(x1, x1*astar, astar2), 
  #                  target = c(n1*colMeans(x0), rep(0, m+1)))
  # wts <- mod$weights

  out <- sapply(a.vals, kern_est_ipw, y = y, a = a, bw = bw, weights = wts, se.fit = TRUE)

  return(out)
  
}

dr_transport <- function(a, y, x1, x0, s, offset = NULL, weights = NULL, family = gaussian(), 
                         df = 4, bw = 1,  a.vals = seq(min(a), max(a), length.out = 100),
                         sl.lib = c("SL.mean", "SL.glm", "SL.earth", "SL.glmnet", "SL.ranger")){
  
  if(is.null(offset))
    offset <- rep(0, times = length(y))
  
  if(is.null(weights))
    weights <- rep(1, times = length(y))
  
  colnames(x0) <- colnames(x1)
  n1 <- nrow(x1)
  n0 <- nrow(x0)
  m <- ncol(x1)
  
  x1.mat <- model.matrix(~ ., data = data.frame(x1))
  x0.mat <- model.matrix(~ ., data = data.frame(x0))
  
  # estimate nuisance outcome model with splines
  mumod <- gam(y ~ s(a, 5) + . + a:. - a, data = data.frame(a = a, x1), offset = offset, family = family)
  muhat <- mumod$fitted.values
  
  mhat.mat <- sapply(a.vals, function(a.tmp, ...) {
    xa.tmp <- data.frame(x0, a = a.tmp)
    return(predict(mumod, newdata = xa.tmp, type = "response"))
  })
  
  # weights
  mod <- calibrate(cmat = cbind(x1.mat), target = n1*colMeans(x0.mat))
  rhohat <- mod$weights
  pimod <- SuperLearner(Y = a, X = data.frame(x1), SL.lib = sl.lib, family = gaussian())
  pimod.vals <- c(pimod$SL.predict)

  sd.a <- sd(a - pimod.vals)
  a.std <- c(a - pimod.vals)/sd.a
  dens <- density(a.std)
  pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / sd.a

  phat <- sapply(a, function(a.tmp, ...) {
    a.std.tmp <- c(a.tmp - pimod.vals)/sd.a
    return(mean(approx(x = dens$x, y = dens$y, xout = a.std.tmp)$y / sd.a, na.rm = T))
  })

  wts <- c(rhohat*phat/pihat, rep(1, n0))
  
  # full calibration weights
  # astar <- c(a - mean(a))/var(a)
  # astar2 <- -1 + c(a - mean(a))^2/var(a)
  # mod <- calibrate(cmat = cbind(x1, x1*astar, astar2),
  #                  target = c(n0*colMeans(x0), rep(0, m+1)))
  # wts <- mod$weights
  
  # regression
  psi <- (y - muhat)
  
  out <- sapply(a.vals, kern_est_dr, psi = psi, mhat.mat = mhat.mat, a = a,
                weights = wts, s = s, a.vals = a.vals, bw = bw, se.fit = TRUE)
  
  return(out)
  
}

# tmle_transport <- function(a, y, x1, x0, offset, family, df = 4,
#                            a.vals = seq(min(a), max(a), length.out = 100),
#                            sl.lib = c("SL.mean", "SL.glm", "SL.earth", "SL.glmnet", "SL.ranger")){
#   
#   colnames(x0) <- colnames(x1)
#   n1 <- nrow(x1)
#   n0 <- nrow(x0)
#   
#   # estimate nuisance outcome model with splines
#   mumod <- glm(y ~ ns(a, df) + . + a:., data = data.frame(x1[,-1]), offset = offset, family = family)
#   muhat <- predict(mumod, newdata = data.frame(x1[,-1], a = a), type = "response")
#   
#   # outcomes models given a.vals
#   muhat.mat <- sapply(a.vals, function(a.tmp, ...) {
#     xa.tmp <- data.frame(x0[,-1], a = a.tmp)
#     return(predict(mumod, newdata = xa.tmp, type = "response"))
#   })
#   
#   # estimate nuisance GPS parameters with lm
#   pimod <- SuperLearner(Y = a, X = data.frame(x1[,-1]), SL.library = sl.lib, family = gaussian())
#   pimod.vals1 <- c(pimod$SL.predict)
#   pimod.vals0 <- predict(pimod, newdata = data.frame(x0[,-1]))$pred
#   pi2mod.vals <- mean((a - pimod.vals1)^2)
#   
#   # estimate sampling weights
#   rhomod <- calibrate(cmat = x1, target = colSums(x0))
#   rhohat <- rhomod$weights
#   
#   # parametric density
#   pihat <- dnorm(a, pimod.vals1, sqrt(pi2mod.vals))
#   pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
#     dnorm(a.tmp, c(pimod.vals1, pimod.vals0), sqrt(pi2mod.vals))
#   })
#   phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n1,], na.rm = T)), x = a)$y
#   phat[phat<0] <- 1e-4
#   
#   # nonparametric density
#   # a.std <- c(a - pimod.vals1) / sqrt(pi2mod.vals)
#   # dens <- density(a.std)
#   # pihat <- approx(x = dens$x, y = dens$y, xout = a.std)$y / sqrt(pi2mod.vals)
#   # 
#   # pihat.mat <- sapply(a.vals, function(a.tmp, ...) {
#   #   std <- c(a.tmp - c(pimod.vals1, pimod.vals0)) / sqrt(pi2mod.vals)
#   #   approx(x = dens$x, y = dens$y, xout = std)$y / sqrt(pi2mod.vals)
#   # })
#   # 
#   # phat <- predict(smooth.spline(a.vals, colMeans(pihat.mat[1:n1,], na.rm = T)), x = a)$y
#   # phat[phat<0] <- 1e-4  
#   
#   # TMLE update
#   nsa <- ns(a, df = df + 1, intercept = TRUE)
#   ipw <- phat/pihat
#   trunc0 <- quantile(ipw, 0.01)
#   trunc1 <- quantile(ipw, 0.99)
#   ipw[ipw < trunc0] <- trunc0
#   ipw[ipw > trunc1] <- trunc1
#   base <- nsa*ipw*rhohat
#   new_mod <- glm(y ~ 0 + base, offset = family$linkfun(muhat) + offset, family = family)
#   param <- coef(new_mod)
#   
#   # predict spline basis and impute
#   estimate <- sapply(1:length(a.vals), function(k, ...) {
#     
#     muhat.tmp <- muhat.mat[,k]
#     pihat.tmp <- pihat.mat[,k]
#     a.tmp <- a.vals[k]
#     wts <- c(mean(pihat.tmp[-(1:n1)], na.rm = TRUE)/pihat.tmp[-(1:n1)])
#     wts[wts < trunc0] <- trunc0
#     wts[wts > trunc1] <- trunc1
#     wts2 <- wts*exp(-c(x0%*%rhomod$coefs))
#     mat <- predict(nsa, newx = rep(a.tmp, n0))*wts2
#     return(mean(family$linkinv(family$linkfun(muhat.tmp) + c(mat%*%param)), na.rm = TRUE))
#     
#   })
#   
#   return(list(estimate = estimate, ipw = ipw, ios = rhohat))
#   
# }
