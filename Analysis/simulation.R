library(SuperLearner)

source("~/Github/transport-erc/R/gen_data.R")
source("~/Github/transport-erc/R/transport.R")
source("~/Github/transport-erc/R/calibrate.R")

iter <- 1000
n0 <- 500
n1 <- 500
sig2 <- 4
tau2 <- 2
scenario <- "out-mis"
a.vals <- seq(5, 15, length.out = 61)

set.seed(42)

simDat <- gen_data(n0 = n0, n1 = n1, sig2 = sig2, tau2 = tau2, a.vals = a.vals, scenario = scenario)
a <- simDat$a
s <- simDat$s
x <- simDat$X
y <- simDat$y
ERC <- simDat$ERC

a <- a[s == 1]
x1 <- x[s == 1,]
x0 <- x[s == 0,]
y <- y[s == 1]
family <- gaussian()
offset <- rep(0, length(y))

out <- glm_transport(a = a, y = y, x1 = x1, x0 = x0, offset = offset, 
                     family = family, df = 5, a.vals = a.vals)
ipw <- ipw_transport(a = a, y = y, x1 = x1, x0 = x0, offset = offset, 
                     family = family, df = 5, a.vals = a.vals)
tmle <- tmle_transport(a = a, y = y, x1 = x1, x0 = x0, offset = offset, 
                       family = family, df = 5, a.vals = a.vals)

plot(a.vals, ERC, type = "l", col = "red", ylim = c(-2,6))
lines(a.vals, out, col = "green")
lines(a.vals, ipw, col = "blue")
lines(a.vals, tmle$estimate, col = "purple")
