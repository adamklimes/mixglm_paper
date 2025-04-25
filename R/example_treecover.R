## example_treecover
# Recreation of Hirota et al. 2011 analyses of tree cover stable states
library(mixglm)
library(raster)

# data_loading
datHirotaSA <- read.csv("data/datHirotaSA.csv")
precVar <- raster("data/wc2.1_2.5m_bio_15.tif")

# data_preparation
squeeze <- function(x) x * 0.98 + 0.01
st <- function(x, y = x) (x - mean(y)) / sd(y)
precVarSA <- crop(precVar, extent(c(-85, -30, -40, 20)))
datHirotaSA$PV <- extract(precVarSA, datHirotaSA[, c("x", "y")], method = "bilinear")
  # subsetting
set.seed(10)
datHirotaSAsub <- datHirotaSA[datHirotaSA$MAP < 3700, ]
datHirotaSAsub <- datHirotaSAsub[sample(1:nrow(datHirotaSAsub), 6000), ]
datHirotaSAsub$Tree_proportion <- squeeze(datHirotaSAsub$Tree_proportion)
datHirotaSAsub$MAPst <- st(datHirotaSAsub$MAP)
datHirotaSAsub$PVst <- st(datHirotaSAsub$PV)

# precipitation only model
for (nStates in 5:8){
  set.seed(21)
  modP <- mixglm(
    stateValModels = Tree_proportion ~ MAPst,
    stateProbModels = ~ MAPst,
    statePrecModels = ~ MAPst,
    stateValError = "beta",
    inputData = datHirotaSAsub,
    numStates = nStates,
    mcmcChains = 2,
    setInit = list(intercept_stateVal = c(-1, rep(if (nStates %in% 5:6) 0.5 else 0.2, nStates - 1)),
                   intercept_statePrec = rep(2, nStates),
                   intercept_stateProb = c(0, rep(0.01, nStates - 1)),
                   MAPst_stateProb = rep(0.01, nStates),
                   MAPst_stateVal = rep(0.01, nStates),
                   MAPst_statePrec = rep(0.01, nStates)),
    setPriors = list(stateVal = list(int1 = "dnorm(0, 0.1)", pred = "dnorm(0, 0.1)"),
                     statePrec = list(int = "dnorm(0, 0.1)", pred = "dnorm(0, 0.1)"))
  )
  save(modP, file = paste0("data/analyses/modP", nStates, ".RData"))
  modP$mcmcSamples$WAIC$WAIC
}
load("data/analyses/modP7.RData")

# figures
pred <- predict(modP, threshold = 0.1)
assignColor <- function(dat, col){
  sq <- seq(min(dat, na.rm = TRUE), max(dat, na.rm = TRUE), length.out = length(col)+1)
  col[findInterval(dat, sq, all.inside = TRUE)]
}
colsDens <- assignColor(- pred$obsDat$potentEn + 1, heat.colors(130)[1:100])
colsState <- assignColor(pred$obsDat$distToState, rev(heat.colors(130)[1:100]))
colsTip <- assignColor(pred$obsDat$distToTip, heat.colors(130)[1:100])
colsDepth <- assignColor(pred$obsDat$potentialDepth, heat.colors(130)[1:100])
axesLas <- function(xat = -8:(-4)*10, yat = -3:1*10, drawX = TRUE, drawY = TRUE) {
  if (drawX) axis(1, labels = paste0(xat, "°"), at = xat)
  if (drawY) axis(2, labels = paste0(yat, "°"), at = yat, las = 2)
  box()
}
assignColDir <- function(tab, dat, cols = c("orange", "darkgreen", "brown", "yellow3")){#c("darkgreen", "lightgreen", "orange", "brown")
  tipResp <- tab$resp[tab$state < 0.5 & tab$catSt == 1]
  tipResp <- tipResp[!is.na(tipResp)]
  if (length(tipResp) < 1) NA else
    if (all(tipResp < dat) & all(tipResp > 0.3)) cols[1] else
      if (all(tipResp < dat) & all(tipResp < 0.3)) cols[3] else
        if (all(tipResp > dat) & all(tipResp > 0.3)) cols[2] else
          if (all(tipResp > dat) & all(tipResp < 0.3)) cols[4] else stop("An observation does not fall into any of predefined categories.")
}
assignColState <- function(tab, dat, cols = c("darkgreen", "lightgreen", "brown")){
  stateResp <- tab$resp[tab$state > 0.5 & tab$catSt == 1]
  tipResp <- tab$resp[tab$state < 0.5 & tab$catSt == 1]
  stateResp <- stateResp[!is.na(stateResp)]
  tipResp <- tipResp[!is.na(tipResp)]
  minE <- function(x) if (length(x) == 0) Inf else min(x)
  stateAttr <- c(stateResp[stateResp <= dat & dat - stateResp < minE((dat - tipResp)[dat - tipResp > 0])],
                 stateResp[stateResp > dat & stateResp - dat < minE((tipResp - dat)[tipResp - dat > 0])])
  if (length(stateAttr) < 1) NA else
    if (stateAttr > 0.6) cols[1] else
      if (stateAttr > 0.1) cols[2] else
        if (stateAttr <= 0.1) cols[3]
}
plotSA <- function(col, main = "", drawX = TRUE, drawY = TRUE){
  plot(datHirotaSAsub$x, datHirotaSAsub$y, type = "n",
       xlab = "", ylab = "", main = main, asp = 1, axes = FALSE)
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = "#81B5FF")
  terra::polys(crop(countriesSA, usr), col = "white", border = NA)
  points(datHirotaSAsub$x, datHirotaSAsub$y, col = col, pch = 15, cex = 0.5)
  terra::polys(crop(countriesSA, usr))
  axesLas(drawX = drawX, drawY = drawY)
}
colDir <- unlist(Map(assignColDir, pred$tipStable, datHirotaSAsub$Tree_proportion))
colState <- unlist(Map(assignColState, pred$tipStable, datHirotaSAsub$Tree_proportion))
countriesSA <- geodata::gadm(c("BRA", "ARG", "PRY", "BOL", "SUR", "GUY", "VEN", "URY", "GUF", "ECU", "CUW", "COL", "CHL", "ABW", "TTO", "PER", "FLK", "PAN", "CRI", "NIC", "GRD"), path = "data/gadm/", level = 0, resolution = 2)
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
  rect(x, -31, x2, -28, col = col, border = NA, xpd = TRUE)
}
precFutureSA <- terra::rast("data/CORDEX South America - Total precipitation (PR) mm_day - Long Term (2081-2100) RCP4.5 - Annual (12 models).nc")
precF <- terra::extract(precFutureSA * 365, cbind(datHirotaSAsub$x, datHirotaSAsub$y))$pr
precF[precF > 3700] <- 3700
colState <- unlist(Map(assignColState, pred$tipStable, datHirotaSAsub$Tree_proportion))
predF <- predict(modP, newdata = data.frame(MAPst = st(precF, datHirotaSAsub$MAP)), threshold = 0.1)
colStateF <- unlist(Map(assignColState, predF$tipStable, datHirotaSAsub$Tree_proportion))
colChange <- vapply(paste0(colState, colStateF), switch,
                    brownbrown = NA_character_, lightgreenlightgreen = NA_character_, darkgreendarkgreen = NA_character_,
                    darkgreenlightgreen = "orange",
                    lightgreendarkgreen = "darkgreen",
                    lightgreenbrown = "brown",
                    brownlightgreen = "yellow3",
                    darkgreenbrown = "black",
                    browndarkgreen = "cyan", FUN.VALUE = "aa")

# figure 2
# png("figures/Fig2_results.png", width = 4800, height = 4800, res = 720)
layout(matrix(c(1,2,3,1,4,5), 3))
par(mai = c(0.6, 0.55, 0.05, 0.01))
landscapeMixglm(modP, Tree_proportion ~ MAPst, 0.1, axes = FALSE,
                xlab = "Precipitation (mm/yr)", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(1, labels = c(0:3*1000), at = st(c(0:3*1000), datHirotaSAsub$MAP))
legend("topleft", lwd = 4, col = c("blue", "red"), legend = c("Stable states", "Tipping points"), bty = "n", cex = 0.9)
par(mai = c(0.3, 0.55, 0.05, 0.85))
plotSA(colState, drawX = FALSE)
legend(-32, -19, col = c("darkgreen", "lightgreen", "brown"), pch = 15, legend = c("Forest", "Savanna", "Treeless"), bty = "n", cex = 1.0, title = "  Domain", title.adj = 0, xpd = TRUE)
plotSA(colsDepth)
text(-30.5, -26, "Potential depth", xpd = TRUE, adj = c(0,0))
xsq <- seq(-30.5, -5, length.out = 101)
invisible(Map(colTr, head(xsq, -1), tail(xsq, -1), heat.colors(130)[1:100]))
text(-30.5, -35, "Low", xpd = TRUE, adj = c(0,0))
text(-10, -35, "High", xpd = TRUE, adj = c(1,0))
par(mai = c(0.3, 0.1, 0.05, 1.3))
plotSA(colDir, drawX = FALSE, drawY = FALSE)
legend(-32, -16, col = c("orange", "darkgreen", "brown", "yellow3"), pch = 15, legend = c("Forest -> Savanna", "Savanna -> Forest", "Savanna -> Treeless", "Treeless -> Savanna"), bty = "n", cex = 1.0, title = "  Alternative states", title.adj = 0, xpd = TRUE)
plotSA(colStateF, drawY = FALSE)
legend(-32, -19, col = c("darkgreen", "lightgreen", "brown"), pch = 15, legend = c("Forest              ", "Savanna", "Treeless"), bty = "n", cex = 1.0, title = "  Predicted domain", title.adj = 0, xpd = TRUE)
addL(c(-0.02,-0.02, 0.48, -0.02, 0.48), c(1.02,0.67, 0.67,0.31,0.31), D = TRUE, E = TRUE)

# figure S1
# png("figures/FigS1_uncertainty.png", width = 4800, height = 4800, res = 720)
set.seed(9)
layout(matrix(c(1,1,1,1,1,1,1,1,2), nrow = 1))
par(mai = c(0.8,0.8,0.1,0.1))
landscapeMixglm(modP, Tree_proportion ~ MAPst, 0.1, axes = FALSE,
                xlab = "Precipitation (mm/yr)", ylab = "Tree cover (%)", randomSample = 20,
                col = terrain.colors(12), addPoints = FALSE)
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(1, labels = c(0:3*1000), at = st(c(0:3*1000), datHirotaSAsub$MAP))
par(mai = c(2.5,0,1.8,0.5))
image(1, 1:12, matrix(1:12, nrow = 1), col = terrain.colors(12), axes = FALSE, ann = FALSE)
box()
axis(4, labels = c("Low", "High"), at = c(1.2,11.8), tick = FALSE, line = -0.5, las = 2)
axis(4, labels = "Uncertainty", at = 6, tick = FALSE, line = -0.5, font = 3)

# figure S2
# png("figures/FigS2_modelplot.png", width = 4800, height = 4800, res = 720)
par(mai = c(0.8,0.8,0.1,0.1))
plot(modP, Tree_proportion ~ MAPst, cex = 0.2, byChains = FALSE, xlab = "Precipitation (mm/yr)", ylab = "Tree cover (%)", drawAxes = FALSE)
axis(1, labels = 0:3*1000, at = st(c(0:3*1000), datHirotaSAsub$MAP))
axis(2, labels = 0:5*20, at = 0:5/5, las = 2)

# figure S3
# png("figures/FigS3_landscape.png", width = 4800, height = 4800, res = 720)
layout(matrix(c(1,2,1,3), 2))
landscapeMixglm(modP, Tree_proportion ~ MAPst, 0.1, axes = FALSE,
  xlab = "Precipitation (mm/yr)", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(1, labels = c(0:3*1000), at = st(c(0:3*1000), datHirotaSAsub$MAP))
legend("topleft", lwd = 4, col = c("blue", "red"), legend = c("Stable states", "Tipping points"), bty = "n", cex = 0.9)
lines(rep(st(1000, datHirotaSAsub$MAP), 2), c(1.05, -0.01), lwd = 2, lty = 2, xpd = TRUE)
lines(rep(st(2000, datHirotaSAsub$MAP), 2), c(1.05, -0.01), lwd = 2, lty = 2, xpd = TRUE)
sc1 <- sliceMixglm(modP, Tree_proportion ~ MAPst, byChains = FALSE, plotEst = FALSE,
  value = st(1000, datHirotaSAsub$MAP), doPlot = FALSE)
plot(sc1$resp, sc1[[1]]$mean, type = "l", axes = FALSE, xlab = "Tree cover (%)",
  ylab = "Probability density", ylim = rev(range(sc1[[1]]$mean)))
box(bty = "l")
axis(2, las = 2)
axis(1, label = 0:5*20, at = 0:5/5)
axis(3, label = "Precipitation 1000 mm/yr", at = mean(par("usr")[1:2]), line = 1,
  cex = 2, tick = FALSE, font = 2)
sc2 <- sliceMixglm(modP, Tree_proportion ~ MAPst, byChains = FALSE, plotEst = FALSE,
  value = st(2000, datHirotaSAsub$MAP), doPlot = FALSE)
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
landscapeMixglm(modP, Tree_proportion ~ MAPst, 0, axes = FALSE,
                xlab = "", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(3, labels = "Threshold = 0.0", at = mean(par("usr")[1:2]), tick = FALSE, line = -0.5, cex.axis = 1.3, font = 2)
par(mai = c(0.2,0.6,0.4,0.1))
landscapeMixglm(modP, Tree_proportion ~ MAPst, 0.1, axes = FALSE,
                xlab = "", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(3, labels = "Threshold = 0.1", at = mean(par("usr")[1:2]), tick = FALSE, line = -0.5, cex.axis = 1.3, font = 2)
par(mai = c(0.4,0.6,0.2,0.1))
landscapeMixglm(modP, Tree_proportion ~ MAPst, 0.2, axes = FALSE,
                xlab = "", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(3, labels = "Threshold = 0.2", at = mean(par("usr")[1:2]), tick = FALSE, line = -0.5, cex.axis = 1.3, font = 2)
par(mai = c(0.6,0.6,0,0.1))
landscapeMixglm(modP, Tree_proportion ~ MAPst, 0.3, axes = FALSE,
                xlab = "Precipitation (mm/yr)", ylab = "Tree cover (%)", ylim = c(-0.01,1), col = grey.colors(12, rev = TRUE))
axis(2, labels = 0:5*20, at = squeeze(0:5/5), las = 2)
axis(3, labels = "Threshold = 0.3", at = mean(par("usr")[1:2]), tick = FALSE, line = -0.5, cex.axis = 1.3, font = 2)
axis(1, labels = c(0:3*1000), at = st(c(0:3*1000), datHirotaSAsub$MAP))

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
# precipitation and prec. variation model
set.seed(21)
modPV <- mixglm(
  stateValModels = Tree_proportion ~ MAPst + PVst,
  stateProbModels = ~ MAPst + PVst,
  statePrecModels = ~ MAPst + PVst,
  stateValError = "beta",
  inputData = datHirotaSAsub,
  numStates = 5,
  mcmcChains = 2,
  setInit = list(intercept_stateVal = c(-1, 0.5, 0.5, 0.5, 0.5),
                 intercept_statePrec = c(2, 2, 2, 2, 2),
                 intercept_stateProb = c(0, 0.01, 0.01, 0.01, 0.01),
                 MAPst_stateProb = c(0.01, 0.01, 0.01, 0.01, 0.01),
                 MAPst_stateVal = c(0.01, 0.01, 0.01, 0.01, 0.01),
                 MAPst_statePrec = c(0.01, 0.01, 0.01, 0.01, 0.01),
                 PVst_stateProb = c(0.01, 0.01, 0.01, 0.01, 0.01),
                 PVst_stateVal = c(0.01, 0.01, 0.01, 0.01, 0.01),
                 PVst_statePrec = c(0.01, 0.01, 0.01, 0.01, 0.01)),
  setPriors = list(stateVal = list(int1 = "dnorm(0, 0.1)", pred = "dnorm(0, 0.1)"),
                   statePrec = list(int = "dnorm(0, 0.1)", pred = "dnorm(0, 0.1)"))
)
# save(modPV, file = "data/analyses/modPV.RData")
load("data/analyses/modPV.RData")

pred <- predict(modPV, threshold = 0.1)
colsDens <- assignColor(- pred$obsDat$potentEn + 1, heat.colors(130)[1:100])
colsState <- assignColor(pred$obsDat$distToState, rev(heat.colors(130)[1:100]))
colsTip <- assignColor(pred$obsDat$distToTip, heat.colors(130)[1:100])
colsDepth <- assignColor(pred$obsDat$potentialDepth, heat.colors(130)[1:100])

# png("figures/FigS8_PVpredictions.png", width = 4800, height = 4800*1.2, res = 720)
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
#_
