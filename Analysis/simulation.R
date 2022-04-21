library(SuperLearner)

source("~/Github/transport-erc/R/gen_data.R")
source("~/Github/transport-erc/R/transport.R")
source("~/Github/transport-erc/R/calibrate.R")
source("~/Github/transport-erc/R/kern_est.R")

iter <- 1000
n0 <- 2000
n1 <- 2000
sig2 <- 1
tau2 <- 2
df <- 5
bw <- 0.5
scenario <- "ps-mis"
a.vals <- seq(6, 14, length.out = 201)
sl.lib = c("SL.mean", "SL.glm", "SL.earth", "SL.glmnet", "SL.ranger")

simDat <- gen_data(n0 = n0, n1 = n1, sig2 = sig2, tau2 = tau2, a.vals = a.vals, scenario = scenario)
a <- simDat$a
s <- simDat$s
x <- simDat$X
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
                     family = family, a.vals = a.vals)
dr <- dr_transport(a = a, y = y, x1 = x1, x0 = x0, df = df,
                   bw = bw, family = family, a.vals = a.vals)

plot(a.vals, ERC, type = "l", col = "red", ylim = c(-10,10))
lines(a.vals, out, col = "green")
lines(a.vals, ipw, col = "blue")
lines(a.vals, dr, col = "purple")
