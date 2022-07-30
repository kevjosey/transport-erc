library(SuperLearner)

source("~/Github/transport-erc/R/gen_data.R")
source("~/Github/transport-erc/R/transport.R")
source("~/Github/transport-erc/R/calibrate.R")
source("~/Github/transport-erc/R/kern_est.R")

iter <- 100
n0 <- 2000
n1 <- 1000
sig2 <- 1
tau2 <- 2
df <- 6
bw <- 0.5
scenario <- "base"
a.vals <- seq(6, 14, length.out = 101)
sl.lib = c("SL.mean", "SL.glm")

rslt_mat <- sapply(1:iter, function(i, ...){
  
  print(i)
  
  simDat <- gen_data(n0 = n0, n1 = n1, sig2 = sig2, tau2 = tau2, a.vals = a.vals, scenario = scenario)
  a <- simDat$a
  s <- simDat$s
  x <- simDat$X
  u <- simDat$U
  y <- simDat$y
  ERC <- simDat$ERC
  
  a <- a[s == 1]
  x1 <- x[s == 1,-1]
  x0 <- x[s == 0,-1]
  y <- y[s == 1]
  family <- gaussian()
  offset <- rep(0, length(y))
  weights <- rep(1, length(y))
  
  out <- glm_transport(a = a, y = y, x1 = x1, x0 = x0, df = df,
                       family = family, a.vals = a.vals)
  ipw <- ipw_transport(a = a, y = y, x1 = x1, x0 = x0, df = df, 
                       bw = bw, family = family, a.vals = a.vals, sl.lib = sl.lib)
  dr <- dr_transport(a = a, y = y, x1 = x1, x0 = x0, s = s, df = df, 
                      bw = bw, family = family, a.vals = a.vals, sl.lib = sl.lib)
  
  data.frame(ERC = ERC, out = out, ipw = ipw[1,], dr = dr[1,], dr_var = dr[2,])
  
})

plot(a.vals, rowMeans(do.call(cbind, rslt_mat[1,])), type = "l", lwd = 2,
     col = "green",  ylim = c(-8, 8), main = "Propensity Score Misspecification")
lines(a.vals, rowMeans(do.call(cbind, rslt_mat[2,])), col = "red")
lines(a.vals, rowMeans(do.call(cbind, rslt_mat[3,])), col = "purple")
lines(a.vals, rowMeans(do.call(cbind, rslt_mat[4,])), col = "blue")
lines(a.vals, rowMeans(do.call(cbind, rslt_mat[4,]) - 1.96*sqrt(do.call(cbind, rslt_mat[5,]))), col = "blue", lty = 2)
lines(a.vals, rowMeans(do.call(cbind, rslt_mat[4,]) + 1.96*sqrt(do.call(cbind, rslt_mat[5,]))), col = "blue", lty = 2)
grid(lwd = 2)
legend(x = 6, y = 8, legend = c("TRUE ERF", "G-COMP", "IPW", "DR", "DR 95% CI"),
       col = c("green", "red", "purple", "blue", "blue"), lty = c(1,1,1,1,2))

erc.mat <- matrix(rep(rowMeans(do.call(cbind, rslt_mat[1,])), iter), nrow = length(a.vals))
upper <- do.call(cbind, rslt_mat[4,]) + 1.96*sqrt(do.call(cbind, rslt_mat[5,]))
lower <- do.call(cbind, rslt_mat[4,]) - 1.96*sqrt(do.call(cbind, rslt_mat[5,]))
cp <- rowMeans(erc.mat < upper, erc.mat > lower)
cp
