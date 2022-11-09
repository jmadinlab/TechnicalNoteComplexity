
# SIMULATE SURFACES AND ESTIMATE FRACTAL DIMENSION USING DIFFERENT METHODS
source('R/MidPointFM2D.R')
source('R/est_D_methods.R')
source('R/functions.R')
library(fishualize)
library(matlib)

# the scale of the map
maxlevel = 8
L = 2^maxlevel

# the fractal dimensions to use to generate the map
all_Ds = rep(seq(2.09, 2.99, 0.001), each=100)

# vectors to store the results
all_Ds_est_box_counting = numeric(length(all_Ds))
all_Ds_est_box_counting_int_scales = numeric(length(all_Ds))
all_Ds_est_variation = numeric(length(all_Ds))
all_Ds_est_variation_int_scales = numeric(length(all_Ds))
all_Ds_est_Torres_Pulliza_method = numeric(length(all_Ds))
R_theory <- numeric(length(all_Ds))
R <- numeric(length(all_Ds))
Height <- numeric(length(all_Ds))

# run different D methods on simulated surfaces w/ known D
t0 = proc.time()
for ( i in 1:length(all_Ds) ) {
  cat('**** i = ',i,'\n',sep='')
  
  # create the map
  D = all_Ds[i]
  H = 3-D
  sigma0 = 2
  sigma = sigma0/(0.5^((H/2)*maxlevel))
  #set.seed(19)
  fM = MidPointFM2D(maxlevel,H,sigma)
  
  # calculate fractal dimension  
  D_est_box_counting = est_D_box_counting_method(fM)
  D_est_box_counting_int_scales = est_D_box_counting_method_int_scales(fM)
  D_est_variation = est_D_variation_method(fM,1,L)
  D_est_variation_int_scales = est_D_variation_method_int_scales(fM,1,L)
  Ds_est_Torres_Pulliza_method = est_D_Torres_Pulliza_method(fM, 1, L)
  
  all_Ds_est_box_counting[i] = D_est_box_counting
  all_Ds_est_box_counting_int_scales[i] = D_est_box_counting_int_scales
  all_Ds_est_variation[i] = D_est_variation
  all_Ds_est_variation_int_scales[i] = D_est_variation_int_scales
  all_Ds_est_Torres_Pulliza_method[i] = Ds_est_Torres_Pulliza_method$D_est
  
  R_theory[i] <- R_func(Ds_est_Torres_Pulliza_method$mean_H_min, 1)
  
  R[i] <- surfaceArea(fM) / L^2 

  Height[i] <- Ds_est_Torres_Pulliza_method$mean_H_max  }
t1 = proc.time()
print(t1-t0)

# store output in dataframe and save
tog <- data.frame("true" = all_Ds,
                  "box" = all_Ds_est_box_counting,
                  "box_int" = all_Ds_est_box_counting_int_scales,
                  "variation" = all_Ds_est_variation,
                  "variation_int" = all_Ds_est_variation_int_scales,
                  "torres_pulliza" = all_Ds_est_Torres_Pulliza_method,
                  "R" = R,
                  "R_theory" = R_theory,
                  "H" = Height)
write.csv(tog, "data/DfromSimulatedSurfaces.csv")

# least squared method to define plane
head(tog)
L0 <- 1

Rs <- c("R","R_theory")
Ds <- c("true","box","box_int","variation","variation_int","torres")
grid.lines = 100

for (r in 1:length(Rs)) {
  for (d in 1:length(Ds)) {
    
    # organize data
    mdat <- tog[c(paste0(Rs[r]), "H")]
    colnames(mdat) <- c("R","H")
    mdat$R <- log10(mdat$R ^2 - 1)
    mdat$H <- log10(mdat$H / (L0 * sqrt(2)))
    mdat$m <- 1
    
    # best fit orthogonal plane 
    A <- data.matrix(mdat)
    A <- matrix(c (sum(mdat$R * mdat$R) , sum(mdat$R * mdat$H) , sum(mdat$R) , 
                   sum(mdat$R * mdat$H) , sum(mdat$H * mdat$H) , sum(mdat$H) ,
                   sum(mdat$R) , sum(mdat$H) , nrow(mdat)), nrow = 3, ncol = 3)
    A.i <- inv(A)
    B <- matrix(c (sum(mdat$R * tog[,Ds[d]]) , 
                   sum(mdat$H * tog[,Ds[d]]) , 
                   sum(tog[,Ds[d]])), nrow = 3, ncol = 1)
    fin <- A.i %*% B
    
    f <- function(x, y) { z <- x * fin[1,1] + y * fin[2,1] + fin[3,1] }
    
    # plane for plotting 
    assign(paste0("R.pred_", Rs[[r]]), seq(min(mdat$R), max(mdat$R), length=grid.lines))
    H.pred <- seq(min(mdat$H), max(mdat$H), length=grid.lines)
    mat <- expand.grid(R=seq(min(mdat$R), max(mdat$R), length=grid.lines), 
                       H=H.pred)
    
    assign(paste0("D.mat_", Rs[[r]], "_", Ds[[d]]), 
           matrix(f(mat$R, mat$H), nrow = grid.lines, ncol = grid.lines))
    
    # d from plane vs data d 
    comp <- data.frame("TrueD" = tog[,Ds[d]],
                       "PlaneD" = f(mdat$R, mdat$H))
    comp$R <- comp$TrueD - comp$PlaneD
    assign(paste0("R2_",Rs[[r]], "_", Ds[[d]]), 
           round(1 - ( sum((comp$TrueD - comp$PlaneD)^2) / (sum((comp$TrueD - mean(comp$TrueD))^2))),4)) }
  }

# Lizard island data 
dta_rdh <- read.csv("data/megaplot.csv")
head(dta_rdh)

for (r in 1:length(Rs)) {
 # organize data
  mdat <- dta_rdh[c(paste0(Rs[r]), "H")]
  colnames(mdat) <- c("R","H")
  mdat$R <- log10(mdat$R ^2 - 1)
  mdat$H <- log10(mdat$H / (L0 * sqrt(2)))
  mdat$m <- 1
    
  # best fit orthogonal plane 
  A <- data.matrix(mdat)
  A <- matrix(c (sum(mdat$R * mdat$R) , sum(mdat$R * mdat$H) , sum(mdat$R) , 
                 sum(mdat$R * mdat$H) , sum(mdat$H * mdat$H) , sum(mdat$H) ,
                 sum(mdat$R) , sum(mdat$H) , nrow(mdat)), nrow = 3, ncol = 3)
  A.i <- inv(A)
  B <- matrix(c (sum(mdat$R * dta_rdh$D) , 
                 sum(mdat$H * dta_rdh$D) , 
                 sum(dta_rdh$D)), nrow = 3, ncol = 1)
  fin <- A.i %*% B
    
  f <- function(x, y) { z <- x * fin[1,1] + y * fin[2,1] + fin[3,1] }
  
  # plane for plotting 
  assign(paste0("R.pred_", Rs[[r]], "_lizard"), seq(min(mdat$R), max(mdat$R), length=grid.lines))
  H.pred_lizard <- seq(min(mdat$H), max(mdat$H), length=grid.lines)
  mat <- expand.grid(R=seq(min(mdat$R), max(mdat$R), length=grid.lines), 
                     H=H.pred_lizard)
    
  assign(paste0("D.mat_", Rs[[r]], "_lizard"), 
         matrix(f(mat$R, mat$H), nrow = grid.lines, ncol = grid.lines))
    
  # d from plane vs data d 
  comp <- data.frame("TrueD" = dta_rdh$D,
                     "PlaneD" = f(mdat$R, mdat$H))
  comp$R <- comp$TrueD - comp$PlaneD
  assign(paste0("R2_",Rs[[r]], "_lizard"), 
         round(1 - ( sum((comp$TrueD - comp$PlaneD)^2) / (sum((comp$TrueD - mean(comp$TrueD))^2))),4)) }

# FIGURE 2
pal <- fish(100, option = "Trimma_lantana")

png("output/figure2.png", width = 4.5, height = 6.5, units = "in", res = 300)

par(mar=c(.2, 3.5, .6, .4), mfrow = c(3,2), ps = 10)
scatter3D(log10(dta_rdh$R ^2 - 1), log10(dta_rdh$H / (L0 * sqrt(2))), dta_rdh$D,
          pch = 20, 
          cex = .5,
          col = pal,
          xlab = "Rugosity",
          ylab = "Height range",
          zlab = "Fractal dimension",
          surf = list(x = R.pred_R_lizard, y = H.pred_lizard, z = D.mat_R_lizard, facets = NA, 
                      col = rgb(0,0,0,0.05), fitpoits = dta_rdh$D),
          theta=215, 
          phi=0,
          colkey=FALSE)
print(R2_R_lizard)
text3D(-0.8, 0.8, 2.48, labels = expression(italic(r)^2 == 0.8714), surf = NULL, add = TRUE)
mtext(expression(bold("DEM Rugosity")), side = 3, line = -1, cex = 1)
mtext(expression(bold("Torres-Pulliza")), side = 2, line= 1.5, cex = 1)
mtext("A", side = 3, at = -.5, line = -.5)

scatter3D(log10(dta_rdh$R_theory ^2 - 1), log10(dta_rdh$H / (L0 * sqrt(2))), dta_rdh$D,
          pch = 20, 
          cex = .5,
          col = pal,
          xlab = "Rugosity",
          ylab = "Height range",
          zlab = "Fractal dimension",
          surf = list(x = R.pred_R_theory_lizard, y = H.pred_lizard, z = D.mat_R_theory_lizard, facets = NA, col = rgb(0,0,0,0.05), fitpoits = dta_rdh$D),
          theta=215, 
          phi=0,
          colkey=FALSE)
print(R2_R_theory_lizard)
text3D(-1, 0.45, 2.55, labels = expression(italic(r)^2 == 0.9860), surf = NULL, add = TRUE)
mtext(expression(bold("Height Range Rugosity")), side = 3, line = -1, cex = 1)
mtext("B", side = 3, at = -.5, line = -.5)

scatter3D(log10(tog$R ^2 - 1), log10(tog$H / (L0 * sqrt(2))), tog$box_int,
          pch = 20, 
          cex = .5,
          col = pal,
          xlab = "Rugosity",
          ylab = "Height range",
          zlab = "Fractal dimension",
          surf = list(x = R.pred_R, y = H.pred, z = D.mat_R_box_int, facets = NA, col = rgb(0,0,0,0.1), fitpoits = tog$box_int),
          theta=215, 
          phi=0,
          colkey=FALSE)
print(R2_R_box_int)
text3D(-0.5, 2.3, 2.75, labels = expression(italic(r)^2 == 0.8958), surf = NULL, add = TRUE)
mtext(expression(bold("Intermediate Box")), side = 2, line= 1.5, cex = 1)
mtext("C", side = 3, at = -.5, line = -.5)

scatter3D(log10(tog$R_theory ^2 - 1), log10(tog$H / (L0 * sqrt(2))), tog$box_int,
          pch = 20, 
          cex = .5,
          col = pal,
          xlab = "Rugosity",
          ylab = "Height range",
          zlab = "Fractal dimension",
          surf = list(x = R.pred_R_theory, y = H.pred, z = D.mat_R_theory_box_int, facets = NA, col = rgb(0,0,0,0.1), fitpoits = tog$box_int),
          theta=215, 
          phi=0,
          colkey=FALSE)
print(R2_R_theory_box_int)
text3D(-0.8, 2.3, 2.75, labels = expression(italic(r)^2 == 0.8956), surf = NULL, add = TRUE)
mtext("D", side = 3, at = -.5, line = -.5)

scatter3D(log10(tog$R ^2 - 1), log10(tog$H / (L0 * sqrt(2))), tog$true,
          pch = 20, 
          cex = .5,
          col = pal,
          xlab = "Rugosity",
          ylab = "Height range",
          zlab = "Fractal dimension",
          surf = list(x = R.pred_R, y = H.pred, z = D.mat_R_true, facets = NA, col = rgb(0,0,0,0.1), fitpoits = tog$true),
          theta=215, 
          phi=0,
          colkey=FALSE)
print(R2_R_true)
text3D(-0.5, 2.3, 2.8, labels = expression(italic(r)^2 == 0.9958), surf = NULL, add = TRUE)
mtext(expression(bold("Actual")), side = 2, line= 1.5, cex = 1)
mtext("E", side = 3, at = -.5, line = -.5)

scatter3D(log10(tog$R_theory ^2 - 1), log10(tog$H / (L0 * sqrt(2))), tog$true,
          pch = 20, 
          cex = .5,
          col = pal,
          xlab = "Rugosity",
          ylab = "Height range",
          zlab = "Fractal dimension",
          surf = list(x = R.pred_R_theory, y = H.pred, z = D.mat_R_theory_true, facets = NA, col = rgb(0,0,0,0.1), fitpoits = tog$true),
          theta=215, 
          phi=0,
          colkey=FALSE)
print(R2_R_theory_true)
text3D(-0.8, 2.3, 2.8, labels = expression(italic(r)^2 == 0.9959), surf = NULL, add = TRUE)
mtext("F", side = 3, at = -.5, line = -.5)

dev.off()

