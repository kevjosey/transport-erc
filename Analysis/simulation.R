library(SuperLearner)

source("~/Github/transport-erc/R/gen_data.R")
source("~/Github/transport-erc/R/transport.R")
source("~/Github/transport-erc/R/calibrate.R")
source("~/Github/transport-erc/R/kern_est.R")

iter <- 100
n0 <- 500
n1 <- 500
sig2 <- 1
tau2 <- 2
df <- 6
bw <- 1
scenario <- "out-mis"
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
  out1 <- glm_transport(a = a, y = y, x1 = x1, x0 = x1, df = df,
                       family = family, a.vals = a.vals)
  ipw <- ipw_transport(a = a, y = y, x1 = x1, x0 = x0, df = df,
                       family = family, a.vals = a.vals, sl.lib = sl.lib)
  dr <- dr_transport(a = a, y = y, x1 = x1, x0 = x0, df = df,
                      bw = bw, family = family, a.vals = a.vals, sl.lib = sl.lib)
  dr1 <- dr_transport(a = a, y = y, x1 = x1, x0 = x1, df = df,
                     bw = bw, family = family, a.vals = a.vals, sl.lib = sl.lib)
  
  data.frame(ERC = ERC, out = out, ipw = ipw, dr = dr, out1 = out1, dr1 = dr1)
  
  
})

plot(a.vals, rowMeans(do.call(cbind, rslt_mat[1,])), type = "l", lwd = 2,
     col = "green",  ylim = c(-8, 8), main = "Propensity Score Misspecification")
lines(a.vals, rowMeans(do.call(cbind, rslt_mat[2,])), col = "red")
lines(a.vals, rowMeans(do.call(cbind, rslt_mat[3,])), col = "blue")
lines(a.vals, rowMeans(do.call(cbind, rslt_mat[4,])), col = "purple")
lines(a.vals, rowMeans(do.call(cbind, rslt_mat[5,])), col = "red", lty = 2)
lines(a.vals, rowMeans(do.call(cbind, rslt_mat[6,])), col = "purple", lty = 2)
grid(lwd = 2)
legend(x = 6, y = 8, legend = c("TRUE ERF", "G-COMP", "IPW", "DR", "Ignore Transport"),
       col = c("green", "red", "blue", "purple", "black"), lty = c(1,1,1,1,2))

