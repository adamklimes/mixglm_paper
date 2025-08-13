# figure 1
fn <- function(x, cf = 1) {
  -(dnorm(1:(100*cf), 10*cf, 5*cf)/1.3 + dnorm(1:(100*cf), 20*cf, 5*cf)/1.4 + dnorm(1:(100*cf), (90 - 7/8 * x)*cf, 5*cf))
}
sc <- function(x) (x - min(x)) / max(x - min(x))
mat <- t(sapply(100:1, fn))
cols <- "purple"
matC <- matrix(0, 99, 99)
matC[60:61, ] <- cols

frS <- function(x, a = 1){
  p <- par("usr")
  if (a == 1) diff(p[1:2])*x else diff(p[3:4])*x
}
lineSq <- function(x, y, cols, cfx = 1, i = 6:1){
  id <- 1
  for (ii in i){
    lines(x+frS(ii/60)*cfx, y+frS(ii/60, 2), col = cols[id], lwd = 3)
    id <- id + 1
  }
}
fnAngle <- function(x, y){
  rad2deg <- function(rad) (rad * 180) / (pi)
  rad2deg(acos((x*x)/(2*x*sqrt((x*x+y*y)/4))))
}
plotSC <- function(ylab = "Probability density"){
  par(mai = c(0.4,0.4,0.3,0.1))
  plot(1:1000, fn(61, 10), type = "l", col = cols, lwd = 4, axes = FALSE,
    xlim = c(0, 550), ylim = c(-0.0087, 0), xlab = "",
    ylab = "", main = "Stability curve")
  axis(1, labels = "System state variable", at = mean(par("usr")[1:2]), tick = FALSE)
  axis(2, labels = ylab, at = mean(par("usr")[3:4]), tick = FALSE)
  points(230, fn(61, 10)[230], pch = 16, cex = 3)
  lines(rep(230, 2), c(-0.0015, fn(61, 10)[230]), lwd = 2, lty = 2)
  lines(rep(127, 2), c(-0.003, fn(61, 10)[127]), lwd = 2, lty = 2)
  lines(rep(276, 2), c(-0.0015, fn(61, 10)[276]), lwd = 2, lty = 2)
  arrows(230, -0.0027, x1 = 127, code = 3, length = 0.06, lwd = 1)
  arrows(230, -0.0012, x1 = 276, code = 3, length = 0.06, lwd = 1)
  text(253, -0.0006, "Distance to tipping point")
  rect(50, -0.0019, 320, -0.0023, border = "white", col = "white")
  text(178.5, -0.0021, "Distance to stable state", bg = "white")
  text(127, -0.0078, "Stable state", col = cols)
  text(276, -0.0062, "Tipping\npoint", col = cols) # -0.0052
  text(366, -0.0085, "Stable state", col = cols)
  if (ylab == "Probability density") arrows(par("usr")[1], par("usr")[4], y1 = par("usr")[3], xpd = TRUE, length = 0.1)
  if (ylab == "Potential energy") arrows(par("usr")[1], par("usr")[3], y1 = par("usr")[4], xpd = TRUE, length = 0.1)
  arrows(par("usr")[1] + 5, par("usr")[3], x1 = par("usr")[2], xpd = TRUE, length = 0.1)
}
plotSL <- function(mai = c(0,0.4,0.1,0.1), ylab = "Probability density"){
  par(mai = mai)
  projMat <- persp(1:nrow(mat), 1:ncol(mat), mat, theta = 35, phi = 30,
    xlab = "Climate", ylab = "System state variable",
    zlab = ylab, zlim = c(-0.15,0.07),
    main = "Stability landscape", col = 0)
  xy <- trans3d(61, 1:100, fn(39), projMat)
  lines(xy$x[1:22], xy$y[1:22], col = cols, lwd = 3, lty = 3)
  lines(xy$x[22:44], xy$y[22:44], col = cols, lwd = 3, lty = 1)
  lines(xy$x[44:64], xy$y[44:64], col = cols, lwd = 3, lty = 3)
  lines(xy$x[64:100], xy$y[64:100], col = cols, lwd = 3, lty = 1)
  if (ylab == "Probability density") {
    rect(-0.48,0.19,-0.46,0.15, border = "white", col = "white")
    arrows(-0.3706, -0.24, -0.368, -0.25, length = 0.1, angle = 10)
  }
}
addL <- function(posx, posy, C = TRUE, D = FALSE){
  par(new = TRUE, mfrow = c(1,1), mai = c(0,0,0,0))
  plot(0:1, 0:1, ann = FALSE, axes = FALSE, type = "n")
  text(posx[1], posy[1], "A", cex = 1.5)
  text(posx[2], posy[2], "B", cex = 1.5)
  if (C) text(posx[3], posy[3], "C", cex = 1.5)
  if (D) text(posx[4], posy[4], "D", cex = 1.5)
}

# png("figures/fig1_2panels.png", width = 4800*2, height = 4800, res = 72*17)
layout(matrix(c(1,1,2,2), ncol = 2))
plotSC(ylab = "Potential energy")
plotSL(c(0.0,0.4,0.4,0.1), ylab = "Potential energy")
addL(c(-0.015, 0.525), c(1.01, 1.01))
#_


