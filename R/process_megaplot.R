# Process the megaplot (Trimodal large plot) data
source("R/functions.R")
library(dplyr)

# Load geotif
data <- raster("data/megaplot/trimodal3set_photoscan_2chunk_DEM_7mm_proj_clip.tif")
# These are midpoints for 2x2m squares, which are our "patch"-level samples
mids <- readOGR("data/megaplot/trimodal_patch_grid_midpnts.shp")
# These are the annotated coral colonies, ID'ed to species level
anno <- readOGR("data/megaplot/trimodal_ann_aligned_cleaned.shp")

anno$coordx <- coordinates(anno)[,1]
anno$coordy <- coordinates(anno)[,2]


# Scope (extent), scales of variation, and resolution (grain)
L <- 2 # Scope

store <- data.frame()

for (mid in 1:length(mids)) {
  rep <- mid
  print(mid)
  # Get lower corner of 2x2m bounding box; alternatively, prescribe these values some other way.
  x0 <- mids@data$x[mid] - L/2
  y0 <- mids@data$y[mid] - L/2
  
  tax <- raster::crop(anno, extent(x0, x0 + 2, y0, y0 + 2))

 dta <-  as.data.frame(table(tax$Species))
  colnames(dta) <- c("species", "abu")
  dta$rep <- mid

  # Calculate rugosit, fractal dimension and height range (rdh function)
  store <- rbind(store, dta)
  
}

for (mid in 24:length(mids)) {
  rep <- mid
  print(mid)
  # Get lower corner of 2x2m bounding box; alternatively, prescribe these values some other way.
  x0 <- mids@data$x[mid] - L/2
  y0 <- mids@data$y[mid] - L/2
  
  tax <- crop(anno, extent(x0, x0 + 2, y0, y0 + 2))
  # points(tax, pch=20, cex=0.3)
  dta <-  as.data.frame(table(tax$Species))
  colnames(dta) <- c("species", "abu")
  dta$rep <- mid
  
  # Calculate rugosit, fractal dimension and height range (rdh function)
  store <- rbind(store, dta)
  
}

for (mid in 80:length(mids)) {
  rep <- mid
  print(mid)
  # Get lower corner of 2x2m bounding box; alternatively, prescribe these values some other way.
  x0 <- mids@data$x[mid] - L/2
  y0 <- mids@data$y[mid] - L/2
  
  tax <- crop(anno, extent(x0, x0 + 2, y0, y0 + 2))
  # points(tax, pch=20, cex=0.3)
  dta <-  as.data.frame(table(tax$Species))
  colnames(dta) <- c("species", "abu")
  dta$rep <- mid
  
  # Calculate rugosit, fractal dimension and height range (rdh function)
  store <- rbind(store, dta)
}

for (mid in 103:length(mids)) {
  rep <- mid
  print(mid)
  # Get lower corner of 2x2m bounding box; alternatively, prescribe these values some other way.
  x0 <- mids@data$x[mid] - L/2
  y0 <- mids@data$y[mid] - L/2
  
  tax <- crop(anno, extent(x0, x0 + 2, y0, y0 + 2))
  # points(tax, pch=20, cex=0.3)
  dta <-  as.data.frame(table(tax$Species))
  colnames(dta) <- c("species", "abu")
  dta$rep <- mid
  
  # Calculate rugosit, fractal dimension and height range (rdh function)
  store <- rbind(store, dta)
}

groups <- data.frame(
  rep = 1:100,
  rep2 = rep(1:50, each = 2),
  rep3 = rep(1:25, each = 4),
  rep4 = rep(1:10, each = 10)
)

library(tidyverse)
out <- inner_join(store, groups)

a1 <- out  %>%
  group_by(rep) %>%
  summarize(species = length(unique(species)), abu = sum(abu)) %>%
  mutate(pa = 4, id = paste0("s1_", rep))

a2 <- out  %>%
  group_by(rep2) %>%
  summarize(species = length(unique(species)), abu = sum(abu)) %>%
  mutate(pa = 8, id = paste0("s2_", rep2))

a3 <- out  %>%
  group_by(rep3) %>%
  summarize(species = length(unique(species)), abu = sum(abu)) %>%
  mutate(pa = 16, id = paste0("s3_", rep3))

a4 <- out  %>%
  group_by(rep4) %>%
  summarize(species = length(unique(species)), abu = sum(abu)) %>%
  mutate(pa = 40, id = paste0("s4_", rep4))

a5 <- out  %>%
  summarize(species = length(unique(species)), abu = sum(abu)) %>%
  mutate(id = "s5_1", pa = 400)

avar <- bind_rows(list(a1, a2, a3, a4, a5))

dta_rdh <- read.csv("data/megaplot.csv") %>% inner_join(groups)



c1 <- dta_rdh  %>%
  group_by(rep) %>%
  summarize(D = mean(D), R = mean(R)) %>%
  mutate(id = paste0("s1_", rep))

c2 <- dta_rdh  %>%
  group_by(rep2) %>%
  summarize(D = mean(D), R = mean(R)) %>%
  mutate(id = paste0("s2_", rep2))

c3 <- dta_rdh  %>%
  group_by(rep3) %>%
  summarize(D = mean(D), R = mean(R)) %>%
  mutate(id = paste0("s3_", rep3))

c4 <- dta_rdh  %>%
  group_by(rep4) %>%
  summarize(D = mean(D), R = mean(R)) %>%
  mutate( id = paste0("s4_", rep4))

c5 <- dta_rdh  %>%
  summarize(D = mean(D), R = mean(R)) %>%
  mutate(id = "s5_1")

cvar <- bind_rows(list(c1, c2, c3, c4))

comb <- left_join(avar, cvar) %>%
  mutate(sa = pa*R) 

ggplot(comb) +
  geom_point(aes(x = log(sa), y = log(species), color = R))


ggplot(comb) +
  geom_point(aes(x = log(pa), y = log(species)))

ggplot(comb) +
  geom_point(aes(x = (R), y = D, color = log(species/sa))) +
  scale_color_gradient(low = "red", high = "blue")


mod <- lm(log(species) ~ R + I(R^2) + log(sa) + D, data = comb)
mod


nd <- expand.grid(R = seq(1, 2.5, 0.05), 
                  D = seq(2, 2.6, 0.01), 
                  sa = 4:100) %>%
  as.data.frame() 

pred <- nd %>%
  mutate(logS = predict(mod, newdata = nd)) %>%
  as.data.frame() 

library(fishualize)

a <- ggplot(pred[pred$D == 2.5,]) +
  geom_raster(aes(x = (sa), y = R,  fill = exp(logS))) +
  scale_fill_fish(direction = 1, option = "Trimma_lantana") +
  geom_hline(yintercept = 1.2, linetype = 2, size = 1) +
  geom_hline(yintercept = 2, linetype = 3, size = 1) +
  theme_classic() +
  guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0,
                                barwidth = unit(1, 'lines'), 
                                barheight = unit(20, 'lines')))  +
  labs(fill = "S", x = "Surface area") +
  theme(text = element_text(size = 16), legend.position = "right")
a
b <- ggplot(pred[pred$D == 2.5 & pred$R ==1.2,]) +
  geom_line(aes(x = (sa), y = exp(logS)), size = 1, linetype = 2) +
  scale_fill_fish(direction = 1, option = "Trimma_lantana") +
  scale_y_continuous(limits = c(5,200)) +
  theme_classic() +
  labs(fill = "S", x = "Surface area", y = "S", title = "R = 1.2") +
  theme(text = element_text(size = 16))

b

c <- ggplot(pred[pred$D == 2.5 & pred$R == 2,]) +
  geom_line(aes(x = (sa), y = exp(logS)), size = 1, linetype = 3) +
  scale_fill_fish(direction = 1, option = "Trimma_lantana") +
  scale_y_continuous(limits = c(5,200)) +
  theme_classic() +
  labs(fill = "S", x = "Surface area", y = "S", title = "R = 2") +
  theme(text = element_text(size = 16))

c

library(patchwork)

layout <- 
  "AAB
   AAC"

a + c + b + patchwork::plot_layout(design = layout) +
  plot_annotation(tag_levels = "A" ) & theme(text = element_text(size = 18))

ggsave("output/Ra_species_fixedD.png")


