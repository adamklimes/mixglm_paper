## simulations
library(mixglm)

# functions
overlap <- function(dist){
  fun <- function(x, dist) pmin(dnorm(x, 0, 1), dnorm(x, dist, 1))
  integrate(fun, -Inf, Inf, dist)$value
}
# overlap(4)
simData <- function(n = 200, dist = 4, shift = 0){
  group <- rbinom(n, 1, 0.5)
  xx <- rep(NA, n)
  xx[group == 0] <- runif(sum(group == 0), 0, 1)
  xx[group == 1] <- runif(sum(group == 1), 0 + shift, 1 + shift)
  yy <- rnorm(n, group * dist)
  data.frame(yy, xx, group)
}
runMod <- function(dat){
  numStates <- 2
  mod <- mixglm(
    stateValModels = yy ~ xx,
    stateProbModels = ~ xx,
    statePrecModels = ~ xx,
    inputData = dat,
    numStates = numStates,
    mcmcChains = 2,
    setInit = list(intercept_stateVal = runif(numStates, 0, 2),
                   intercept_statePrec = rnorm(numStates, 2, 0.5),
                   intercept_stateProb = c(0, rnorm(numStates - 1, 0, 0.1)),
                   xx_stateVal = rnorm(numStates, 0, 0.1),
                   xx_statePrec = rnorm(numStates, 0, 0.1),
                   xx_stateProb = rnorm(numStates, 0, 0.1))
  )
  mod
}
oneRun <- function(n = 200, dist = 4, shift = 0, doPlot = FALSE){
  dat <- simData(n, dist, shift)
  mod <- runMod(dat)
  sel <- (dat$xx >= shift) & (dat$xx <= 1)
  fails <- sum(is.na(predict(mod, threshold = 0)$obsDat$distToTip)[sel]) / sum(sel)
  if (doPlot) landscapeMixglm(yy ~ xx, mod)
  c(fails = fails, n = n, dist = dist, shift = shift)
}

# simulations
R <- 100
ns <- c(50, 100, 200, 300)
dists <- c(2.5,3,4,5)
shifts <- c(0,0.2,0.5,0.8)

for (i in 1:20){
  write(oneRun(shift = shifts[4]), "data/simulations/results.txt", append = TRUE)
  print(i)
}

sim1 <- do.call(rbind, Map(oneRun, n = rep(ns[1], each = R)))
simTry <- do.call(rbind, Map(oneRun, n = rep(c(100,150,200), each = 3)))

# figures
simRes <- read.table(file = "data/simulations/results.txt", header = TRUE)
auxPowerN <- tapply(simRes$fails == 0, paste(simRes$n, simRes$dist, simRes$shift), mean)
powerN <- data.frame(auxPowerN, do.call(rbind, lapply(strsplit(names(auxPowerN), " "), as.numeric)))
colnames(powerN) <- c("power", "n", "dist", "shift")
rownames(powerN) <- NULL

powerPlot <- function(powerN, selVar, xlab, ylab = ""){
  nSel <- if (selVar != "n") powerN$n == 200 else TRUE
  distSel <- if (selVar != "dist") powerN$dist == 4 else TRUE
  shiftSel <- if (selVar != "shift") powerN$shift == 0 else TRUE
  sel <- nSel & distSel & shiftSel
  ord <- order(powerN[sel, ][[selVar]])
  plot(range(powerN[sel, ][[selVar]]), 0:1, type = "n", xlab = xlab, ylab = ylab, axes = FALSE, yaxs = "i")
  abline(h = 1:4/5, lty = 2, col = "grey")
  box()
  axis(1)
  lines(powerN[sel, ][[selVar]][ord], powerN[sel, ]$power[ord])
  points(powerN[sel, ][[selVar]][ord], powerN[sel, ]$power[ord], pch = 16, cex = 1.3)
}
datSim <- list(A1 = simData(50, 4, 0)[,c(2:1, 3)], A2 = simData(100, 4, 0)[,c(2:1, 3)], A3 = simData(200, 4, 0)[,c(2:1, 3)], A4 = simData(300, 4, 0)[,c(2:1, 3)],
  B1 = simData(200, 2.5, 0)[,c(2:1, 3)], B2 = simData(200, 3, 0)[,c(2:1, 3)], B3 = simData(200, 4, 0)[,c(2:1, 3)], B4 = simData(200, 5, 0)[,c(2:1, 3)],
  C1 = simData(200, 4, 0)[,c(2:1, 3)], C2 = simData(200, 4, 0.2)[,c(2:1, 3)], C3 = simData(200, 4, 0.5)[,c(2:1, 3)], C4 = simData(200, 4, 0.8)[,c(2:1, 3)])
lineL <- function(){
  abline(v = c(0, 1.6,3.2,4.8))
  lines(c(0, 1.2), rep(par("usr")[3], 2))
  lines(c(1.6, 2.8), rep(par("usr")[3], 2))
  lines(c(3.2, 4.4), rep(par("usr")[3], 2))
  lines(c(4.8, 6), rep(par("usr")[3], 2))
}
plotIlust <- function(dat1, dat2, dat3, dat4, cf = c(1,1,1,1)){
  cols = c("grey50", "grey20")
  plot(c(0,6), range(unlist(datSim)), type = "n", axes = FALSE, ann = FALSE)
  points(dat1$xx/cf[1]+0.1, dat1$yy, cex = 0.5, pch = 16, col = cols[dat1$group + 1])
  points(dat2$xx/cf[2]+1.7, dat2$yy, cex = 0.5, pch = 16, col = cols[dat2$group + 1])
  points(dat3$xx/cf[3]+3.3, dat3$yy, cex = 0.5, pch = 16, col = cols[dat3$group + 1])
  points(dat4$xx/cf[4]+4.9, dat4$yy, cex = 0.5, pch = 16, col = cols[dat4$group + 1])
  lineL()
}

# png("figures/FigS4_power.png", height = 480*7, width = 4800, res = 720)
set.seed(10)
layout(matrix(c(1,4,4,2,5,5,3,6,6), ncol = 3))
par(mai = c(0.1,0.6,0.1,0))
plotIlust(datSim$A1, datSim$A2, datSim$A3, datSim$A4)
axis(2, las = 2, at = mean(par("usr")[3:4]), labels = "Example\ndatasets", tick = FALSE, line = -0.5)
par(mai = c(0.1,0.3,0.1,0.3))
plotIlust(datSim$B1, datSim$B2, datSim$B3, datSim$B4)
par(mai = c(0.1,0,0.1,0.6))
plotIlust(datSim$C1, datSim$C2, datSim$C3, datSim$C4, cf = c(1,1.2,1.5,1.8))
par(mai = c(0.6,0.6,0.1,0))
powerPlot(powerN, "n", "Sample size", "Power")
axis(2, las = 2)
par(mai = c(0.6,0.3,0.1,0.3))
powerPlot(powerN, "dist", "Distance (SDs)")
par(mai = c(0.6,0,0.1,0.6))
powerPlot(powerN, "shift", "Shift")

#_
