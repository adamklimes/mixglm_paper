## comparison_packages
library(earlywarnings)

simData <- function(n = 200, dist = 4, shift = 0){
  group <- rbinom(n, 1, 0.5)
  xx <- rep(NA, n)
  xx[group == 0] <- runif(sum(group == 0), 0, 1)
  xx[group == 1] <- runif(sum(group == 1), 0 + shift, 1 + shift)
  yy <- rnorm(n, group * dist)
  data.frame(yy, xx, group)
}
R <- 100
ns <- c(50, 100, 200, 300)
dists <- c(2.5,3,5) # 4 is default
shifts <- c(0.2,0.5,0.8) # 0 is default
oneRunMP <- function(n = 200, dist = 4, shift = 0){
  gridSize <- 50 # 50 is default
  dat <- simData(n, dist, shift)
  res <- movpotential_ews(dat$yy, dat$xx, grid.size = gridSize, plot.contours = FALSE)
  sel <- (res$res$pars[, 1] >= shift) & (res$res$pars[, 1] <= 1)
  fails <- sum(rowSums(res$res$mins[sel, ]) == 0) / sum(sel)
  c(fails = fails, n = n, dist = dist, shift = shift)
}

sim1 <- do.call(rbind, Map(oneRunMP, n = rep(ns, each = R)))
sim2 <- do.call(rbind, Map(oneRunMP, dist = rep(dists, each = R)))
sim3 <- do.call(rbind, Map(oneRunMP, shift = rep(shifts, each = R)))
simResMP <- data.frame(rbind(sim1, sim2, sim3))
# write.table(simResMP, file = "data/simulations/resultsMP.txt")
simResMP <- read.table(file = "data/simulations/resultsMP.txt", header = TRUE)
simRes <- read.table(file = "data/simulations/results.txt", header = TRUE)

auxPowerNMP <- tapply(simResMP$fails <= 0.1, paste(simResMP$n, simResMP$dist, simResMP$shift), mean)
powerNMP <- data.frame(auxPowerNMP, do.call(rbind, lapply(strsplit(names(auxPowerNMP), " "), as.numeric)))
colnames(powerNMP) <- c("power", "n", "dist", "shift")
rownames(powerNMP) <- NULL
auxPowerNMP2 <- tapply(simResMP$fails == 0, paste(simResMP$n, simResMP$dist, simResMP$shift), mean)
powerNMP2 <- data.frame(auxPowerNMP2, do.call(rbind, lapply(strsplit(names(auxPowerNMP2), " "), as.numeric)))
colnames(powerNMP2) <- c("power", "n", "dist", "shift")
rownames(powerNMP2) <- NULL
auxPowerN <- tapply(simRes$fails == 0, paste(simRes$n, simRes$dist, simRes$shift), mean)
powerN <- data.frame(auxPowerN, do.call(rbind, lapply(strsplit(names(auxPowerN), " "), as.numeric)))
colnames(powerN) <- c("power", "n", "dist", "shift")
rownames(powerN) <- NULL

powerPlot <- function(powerN, selVar, xlab, ylab = "", mainPlot = TRUE, col = "black"){
  nSel <- if (selVar != "n") powerN$n == 200 else TRUE
  distSel <- if (selVar != "dist") powerN$dist == 4 else TRUE
  shiftSel <- if (selVar != "shift") powerN$shift == 0 else TRUE
  sel <- nSel & distSel & shiftSel
  ord <- order(powerN[sel, ][[selVar]])
  if (mainPlot){
    plot(range(powerN[sel, ][[selVar]]), 0:1, type = "n", xlab = xlab, ylab = ylab, axes = FALSE, yaxs = "i")
    abline(h = 1:4/5, lty = 2, col = "grey")
    box()
    axis(1)
  }
  lines(powerN[sel, ][[selVar]][ord], powerN[sel, ]$power[ord], col = col)
  points(powerN[sel, ][[selVar]][ord], powerN[sel, ]$power[ord], pch = 16, cex = 1.3, col = col)
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

# png("figures/FigS5_power.png", height = 480*7, width = 4800, res = 720)
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
powerPlot(powerNMP, "n", "Sample size", "Power", col = "orange")
powerPlot(powerN, "n", "", mainPlot = FALSE)
powerPlot(powerNMP2, "n", "", mainPlot = FALSE, col = "red")
axis(2, las = 2)
par(mai = c(0.6,0.3,0.1,0.3))
powerPlot(powerNMP, "dist", "Distance (SDs)", col = "orange")
powerPlot(powerN, "dist", "", mainPlot = FALSE)
powerPlot(powerNMP2, "dist", "", mainPlot = FALSE, col = "red")
par(mai = c(0.6,0,0.1,0.6))
powerPlot(powerNMP, "shift", "Shift", col = "orange")
powerPlot(powerN, "shift", "", mainPlot = FALSE)
powerPlot(powerNMP2, "shift", "", mainPlot = FALSE, col = "red")

#_
