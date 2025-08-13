### example_treecover
# Recreation of Hirota et al. 2011/ Flores et al. 2024 analyses of tree cover stable states
library(terra)
library(mixglm)

treeCover

# plot(treeCover$precip, treeCover$treeCover, cex = 0.2)
# plot(treeCover$temp, treeCover$treeCover, cex = 0.2)
 # exclusion of extreme values - there are not enough datapoints to estimates stable states for:
 #  - precipitation > 4000 mm/year
 #  - mean annual temperature < 0°C
datTC <- treeCover[treeCover$precip < 4000 & treeCover$temp > 0, ] # 4860 observations
 # tree cover values are squeezed because beta distribution does not allow 0 and 1
squeeze <- function(x) x * 0.98 + 0.01
datTC$treeCoverProp <- squeeze(datTC$treeCover / 100)
 # standardization of predictors
st <- function(x, y = x) (x - mean(y)) / sd(y)
datTC$precipSt <- st(datTC$precip)
datTC$tempSt <- st(datTC$temp)

## analyses
# model - precipitation only
numStates <- 7
set.seed(21)
modP <- mixglm(
  stateValModels = treeCoverProp ~ precipSt,
  stateProbModels = ~ precipSt,
  statePrecModels = ~ precipSt,
  stateValError = "beta",
  inputData = datTC,
  numStates = numStates,
  mcmcChains = 2
)
# save(modP, file = paste0("data/analyses/modP", numStates, ".RData"))
load(file = "data/analyses/modP7.RData")

# model - precipitation + temperature
set.seed(21)
modPT <- mixglm(
  stateValModels = treeCoverProp ~ precipSt + tempSt,
  stateProbModels = ~ precipSt + tempSt,
  statePrecModels = ~ precipSt + tempSt,
  stateValError = "beta",
  inputData = datTC,
  numStates = nStates,
  mcmcChains = 2
)
# save(modPT, file = "data/analyses/modPT7.RData")

## figures
pred <- predict(modP, threshold = 0.1)
assignColor <- function(dat, col){
  sq <- seq(min(dat, na.rm = TRUE), max(dat, na.rm = TRUE), length.out = length(col)+1)
  col[findInterval(dat, sq, all.inside = TRUE)]
}
colDom <- c("darkgreen", "lightgreen")
colD <- c("orange", "darkgreen")
colEqi <- c("blue", "#fb8072")

colsDens <- assignColor(- pred$obsDat$potentEn + 1, heat.colors(130)[1:100])
colsState <- assignColor(pred$obsDat$distToState, rev(heat.colors(130)[1:100]))
colsTip <- assignColor(pred$obsDat$distToTip, heat.colors(130)[1:100])
colsDepth <- assignColor(pred$obsDat$potentialDepth, heat.colors(130)[1:100])
axesLas <- function(xat = -5:(-1)*20, yat = -7:1*10, drawX = TRUE, drawY = TRUE) {
  if (drawX) axis(1, labels = paste0(xat, "°"), at = xat)
  if (drawY) axis(2, labels = paste0(yat, "°"), at = yat, las = 2)
  box()
}
assignColDir <- function(tab, dat, cols = colD){
  tipResp <- tab$resp[tab$state < 0.5 & tab$catSt == 1]
  tipResp <- tipResp[!is.na(tipResp)]
  if (length(tipResp) < 1) NA else
    if (all(tipResp < dat)) cols[1] else
      if (all(tipResp > dat)) cols[2] else stop("An observation does not fall into any of predefined categories.")
}
assignColState <- function(tab, dat, cols = colDom){
  stateResp <- tab$resp[tab$state > 0.5 & tab$catSt == 1]
  tipResp <- tab$resp[tab$state < 0.5 & tab$catSt == 1]
  stateResp <- stateResp[!is.na(stateResp)]
  tipResp <- tipResp[!is.na(tipResp)]
  minE <- function(x) if (length(x) == 0) Inf else min(x)
  stateAttr <- c(stateResp[stateResp <= dat & dat - stateResp < minE((dat - tipResp)[dat - tipResp > 0])],
                 stateResp[stateResp > dat & stateResp - dat < minE((tipResp - dat)[tipResp - dat > 0])])
  if (length(stateAttr) < 1) NA else
    if (stateAttr > 0.6) cols[1] else cols[2]
}
plotSA <- function(col, main = "", drawX = TRUE, drawY = TRUE){
  plot(datTC$x, datTC$y, type = "n",
       xlab = "", ylab = "", main = main, asp = 1, axes = FALSE)
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = "#81B5FF")
  terra::polys(crop(countriesSA, usr), col = "white", border = NA)
  points(datTC$x, datTC$y, col = col, pch = 15, cex = 0.5)
  terra::polys(crop(countriesSA, usr))
  axesLas(drawX = drawX, drawY = drawY)
}
colDir <- unlist(Map(assignColDir, pred$tipStable, datTC$treeCoverProp))
colState <- unlist(Map(assignColState, pred$tipStable, datTC$treeCoverProp))
countriesSA <- geodata::gadm(c("BRA", "ARG", "PRY", "BOL", "SUR", "GUY", "VEN", "URY", "GUF", "ECU", "CUW", "COL", "CHL", "ABW", "TTO", "PER", "FLK", "PAN", "CRI", "NIC", "GRD", "HND", "SLV", "BRB", "VCT", "LCA", "MTQ", "DMA", "SGS"), path = "data/gadm/", level = 0, resolution = 2)
addL <- function(posx, posy, C = TRUE, D = FALSE, E = FALSE, cex = 1.5){
  par(new = TRUE, mfrow = c(1,1), mai = c(0,0,0,0))
  plot(0:1, 0:1, ann = FALSE, axes = FALSE, type = "n")
  text(posx[1], posy[1], "A", cex = cex)
  text(posx[2], posy[2], "B", cex = cex)
  if (C) text(posx[3], posy[3], "C", cex = cex)
  if (D) text(posx[4], posy[4], "D", cex = cex)
  if (E) text(posx[5], posy[5], "E", cex = cex)
}
colTr <- function(x, x2, col){
  rect(x, -26, x2, -23, col = col, border = NA, xpd = TRUE)
}
colState <- unlist(Map(assignColState, pred$tipStable, datTC$treeCoverProp))
 # future precip is used for domain prediction which is forest for all values over 4000
datTC$precip2100[datTC$precip2100 > 4000] <- 4000
predF <- predict(modP, newdata = data.frame(precipSt = st(datTC$precip2100, datTC$precip)), threshold = 0.1)
colStateF <- unlist(Map(assignColState, predF$tipStable, datTC$treeCoverProp))

# info
 # Fig. 2
sum(colDir == colD[1], na.rm = TRUE) / sum(colState == colDom[1])
sum(colDir == colD[2], na.rm = TRUE) / sum(colState == colDom[2])
sum(colState == colDom[2] & colStateF == colDom[1])/ sum(colState == colDom[2])
sum(colState == colDom[1] & colStateF == colDom[2])/ sum(colState == colDom[1])
 # Fig. S7 - alt. stable state
sum(!is.na(colDir)) / length(colDir)
 # Supplements
sum(colState == colDom[1]) / length(colState)

# figure 2
# png("figures/Fig2_results.png", width = 4800, height = 4800, res = 720)
layout(matrix(c(1,2,3,1,4,5), 3))
par(mai = c(0.6, 0.55, 0.05, 0.01))
landscapeMixglm(modP, treeCoverProp ~ precipSt, 0.1, axes = FALSE,
                xlab = "Precipitation (mm/yr)", ylab = "Tree cover (%)",
                ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE), eqiCol = colEqi)
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(1, labels = c(0:3*1000), at = st(c(0:3*1000), datTC$precip))
legend("topleft", lwd = 4, col = colEqi, legend = c("Stable states", "Tipping points"), bty = "n", cex = 0.9)
par(mai = c(0.3, 0.55, 0.05, 0.85))
plotSA(colState, drawX = FALSE)
legend(-20, -16, col = colDom, pch = 15, legend = c("Forest", "Savanna"), bty = "n", cex = 1.0, title = "  Domain", title.adj = 0, xpd = TRUE)
plotSA(colsDepth)
text(-17.5, -21, "Potential depth", xpd = TRUE, adj = c(0,0))
xsq <- seq(-16.5, 11, length.out = 101)
invisible(Map(colTr, head(xsq, -1), tail(xsq, -1), heat.colors(130)[1:100]))
text(-14.5, -30, "Low", xpd = TRUE, adj = c(0,0))
text(9, -30, "High", xpd = TRUE, adj = c(1,0))
par(mai = c(0.3, 0.1, 0.05, 1.3))
plotSA(colDir, drawX = FALSE, drawY = FALSE)
legend(-20, -16, col = colD, pch = 15, legend = c("Forest -> Savanna", "Savanna -> Forest"), bty = "n", cex = 1.0, title = "  Alternative states", title.adj = 0, xpd = TRUE)
plotSA(colStateF, drawY = FALSE)
legend(-20, -16, col = colDom, pch = 15, legend = c("Forest              ", "Savanna"), bty = "n", cex = 1.0, title = "  Predicted domain", title.adj = 0, xpd = TRUE)
addL(c(-0.02,-0.02, 0.48, -0.02, 0.48), c(1.02,0.67, 0.67,0.31,0.31), D = TRUE, E = TRUE)

# figure S1
# png("figures/FigS1_uncertainty.png", width = 4800, height = 4800, res = 720)
set.seed(9)
layout(matrix(c(1,1,1,1,1,1,1,1,2), nrow = 1))
par(mai = c(0.8,0.8,0.1,0.1))
landscapeMixglm(modP, treeCoverProp ~ precipSt, 0.1, axes = FALSE,
                xlab = "Precipitation (mm/yr)", ylab = "Tree cover (%)", randomSample = 20,
                col = terrain.colors(12), addPoints = FALSE)
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(1, labels = c(0:3*1000), at = st(c(0:3*1000), datTC$precip))
par(mai = c(2.5,0,1.8,0.5))
image(1, 1:12, matrix(1:12, nrow = 1), col = terrain.colors(12), axes = FALSE, ann = FALSE)
box()
axis(4, labels = c("Low", "High"), at = c(1.2,11.8), tick = FALSE, line = -0.5, las = 2)
axis(4, labels = "Uncertainty", at = 6, tick = FALSE, line = -0.5, font = 3)

# figure S2
# png("figures/FigS2_modelplot.png", width = 4800, height = 4800, res = 720)
par(mai = c(0.8,0.8,0.1,0.1))
plot(modP, treeCoverProp ~ precipSt, cex = 0.2, byChains = FALSE, xlab = "Precipitation (mm/yr)", ylab = "Tree cover (%)", drawAxes = FALSE)
axis(1, labels = 0:3*1000, at = st(c(0:3*1000), datTC$precip))
axis(2, labels = 0:5*20, at = 0:5/5, las = 2)

# figure S3
# png("figures/FigS3_landscape.png", width = 4800, height = 4800, res = 720)
layout(matrix(c(1,2,1,3), 2))
landscapeMixglm(modP, treeCoverProp ~ precipSt, 0.1, axes = FALSE, eqiCol = colEqi,
                xlab = "Precipitation (mm/yr)", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(1, labels = c(0:3*1000), at = st(c(0:3*1000), datTC$precip))
legend("topleft", lwd = 4, col = colEqi, legend = c("Stable states", "Tipping points"), bty = "n", cex = 0.9)
lines(rep(st(1000, datTC$precip), 2), c(1.05, -0.01), lwd = 2, lty = 2, xpd = TRUE)
lines(rep(st(2000, datTC$precip), 2), c(1.05, -0.01), lwd = 2, lty = 2, xpd = TRUE)
sc1 <- sliceMixglm(modP, treeCoverProp ~ precipSt, byChains = FALSE, plotEst = FALSE,
                   value = st(1000, datTC$precip), doPlot = FALSE)
plot(sc1$resp, sc1[[1]]$mean, type = "l", axes = FALSE, xlab = "Tree cover (%)",
     ylab = "Probability density", ylim = rev(range(sc1[[1]]$mean)))
box(bty = "l")
axis(2, las = 2)
axis(1, label = 0:5*20, at = 0:5/5)
axis(3, label = "Precipitation 1000 mm/yr", at = mean(par("usr")[1:2]), line = 1,
     cex = 2, tick = FALSE, font = 2)
sc2 <- sliceMixglm(modP, treeCoverProp ~ precipSt, byChains = FALSE, plotEst = FALSE,
                   value = st(2000, datTC$precip), doPlot = FALSE)
plot(sc2$resp, sc2[[1]]$mean, type = "l", axes = FALSE, xlab = "Tree cover (%)",
     ylab = "Probability density", ylim = rev(range(sc2[[1]]$mean)))
box(bty = "l")
axis(2, las = 2)
axis(1, label = 0:5*20, at = 0:5/5)
axis(3, label = "Precipitation 2000 mm/yr", at = mean(par("usr")[1:2]), line = 1,
     cex = 2, tick = FALSE, font = 2)
par(new = TRUE, mfrow = c(1,1), mai = c(0,0,0,0))
plot(0:1, 0:1, ann = FALSE, axes = FALSE, type = "n")
addL(c(-0.01,-0.01,0.525), c(0.99,0.46,0.46))

# figure S4
# power analysis; script "simulations.R"

# figure S5
# power analysis - comparison; script "comparison_packages.R"

# figure S6
# png("figures/FigS6_thresholds.png", width = 480*7, height = 480*14, res = 720)
par(mfrow = c(4,1), mai = c(0,0.6,0.6,0.1))
landscapeMixglm(modP, treeCoverProp ~ precipSt, 0, axes = FALSE,
                xlab = "", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(3, labels = "Threshold = 0.0", at = mean(par("usr")[1:2]), tick = FALSE, line = -0.5, cex.axis = 1.3, font = 2)
par(mai = c(0.2,0.6,0.4,0.1))
landscapeMixglm(modP, treeCoverProp ~ precipSt, 0.1, axes = FALSE,
                xlab = "", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(3, labels = "Threshold = 0.1", at = mean(par("usr")[1:2]), tick = FALSE, line = -0.5, cex.axis = 1.3, font = 2)
par(mai = c(0.4,0.6,0.2,0.1))
landscapeMixglm(modP, treeCoverProp ~ precipSt, 0.2, axes = FALSE,
                xlab = "", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(3, labels = "Threshold = 0.2", at = mean(par("usr")[1:2]), tick = FALSE, line = -0.5, cex.axis = 1.3, font = 2)
par(mai = c(0.6,0.6,0,0.1))
landscapeMixglm(modP, treeCoverProp ~ precipSt, 0.3, axes = FALSE,
                xlab = "Precipitation (mm/yr)", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(3, labels = "Threshold = 0.3", at = mean(par("usr")[1:2]), tick = FALSE, line = -0.5, cex.axis = 1.3, font = 2)
axis(1, labels = c(0:3*1000), at = st(c(0:3*1000), datTC$precip))

# figure S7
# png("figures/FigS7_predictions.png", width = 4800, height = 4800*1.2, res = 720)
layout(matrix(c(1,1,2,2,5,3,3,4,4,5), ncol = 2))
par(mai = c(0.45,0.45,0.4,0.1))
plotSA(colsDens, "Probability density (scaled)")
plotSA(colsState, "Distance to stable state")
plotSA(colsTip, "Distance to tipping point")
plotSA(colsDepth, "Potential depth")
par(mai = c(1.1,2.25,0.3,1.9))
image(1:100, 1, matrix(1:100, ncol = 1), col = rev(heat.colors(130)[1:100]), axes = FALSE, ann = FALSE)
box()
axis(3, labels = c("Legend"), at = 50, tick = FALSE, font = 2)
axis(1, labels = c("High", "Low"), at = c(5,95), tick = FALSE, line = -0.5)
axis(1, labels = c("Low", "High"), at = c(5,95), tick = FALSE, line = 1)
axis(1, labels = c("High", "Low"), at = c(5,95), tick = FALSE, line = 2.5)
axis(1, labels = c("High", "Low"), at = c(5,95), tick = FALSE, line = 4)
axis(1, labels = "Probability density", at = 50, tick = FALSE, line = -0.5, font = 3)
axis(1, labels = "Distance to stable state", at = 50, tick = FALSE, line = 1, font = 3)
axis(1, labels = "Distance to tipping point", at = 50, tick = FALSE, line = 2.5, font = 3)
axis(1, labels = "Potential depth", at = 50, tick = FALSE, line = 4, font = 3)
addL(c(-0.01,0.53,-0.01,0.53), c(1.005,1.005,0.575,0.575), D = TRUE)

# figure S8
load("data/analyses/modPT7.RData")

predPT <- predict(modPT, threshold = 0.1)
colsDensPT <- assignColor(- predPT$obsDat$potentEn + 1, heat.colors(130)[1:100])
colsStatePT <- assignColor(predPT$obsDat$distToState, rev(heat.colors(130)[1:100]))
colsTipPT <- assignColor(predPT$obsDat$distToTip, heat.colors(130)[1:100])
colsDepthPT <- assignColor(predPT$obsDat$potentialDepth, heat.colors(130)[1:100])

# png("figures/FigS8_PTpredictions.png", width = 4800, height = 4800*1.2, res = 720)
layout(matrix(c(1,1,2,2,5,3,3,4,4,5), ncol = 2))
par(mai = c(0.45,0.45,0.4,0.1))
plotSA(colsDensPT, "Probability density (scaled)")
plotSA(colsStatePT, "Distance to stable state")
plotSA(colsTipPT, "Distance to tipping point")
plotSA(colsDepthPT, "Potential depth")
par(mai = c(1.1,2.25,0.3,1.9))
image(1:100, 1, matrix(1:100, ncol = 1), col = rev(heat.colors(130)[1:100]), axes = FALSE, ann = FALSE)
box()
axis(3, labels = c("Legend"), at = 50, tick = FALSE, font = 2)
axis(1, labels = c("High", "Low"), at = c(5,95), tick = FALSE, line = -0.5)
axis(1, labels = c("Low", "High"), at = c(5,95), tick = FALSE, line = 1)
axis(1, labels = c("High", "Low"), at = c(5,95), tick = FALSE, line = 2.5)
axis(1, labels = c("High", "Low"), at = c(5,95), tick = FALSE, line = 4)
axis(1, labels = "Probability density", at = 50, tick = FALSE, line = -0.5, font = 3)
axis(1, labels = "Distance to stable state", at = 50, tick = FALSE, line = 1, font = 3)
axis(1, labels = "Distance to tipping point", at = 50, tick = FALSE, line = 2.5, font = 3)
axis(1, labels = "Potential depth", at = 50, tick = FALSE, line = 4, font = 3)
addL(c(-0.01,0.53,-0.01,0.53), c(1.005,1.005,0.575,0.575), D = TRUE)

#_