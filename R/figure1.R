library(tidyverse)
library(fishualize)
library(patchwork)
library(readr)

# load data
dta_coral <- read_csv("data/coralspecies.csv")
dta_rdh <- read_csv("data/megaplot.csv") %>%
  dplyr::select(rep, R, D)

dta <- left_join(dta_coral, dta_rdh)

min_x <- min(dta_coral$mid_x) - 1
min_y <- min(dta_coral$mid_y) - 1

max_x <- max(dta_coral$mid_x) + 1
max_y <- max(dta_coral$mid_y) + 1

### 2x2m boxes 
comb2 <- dta_coral %>%
  dplyr::group_by(rep, pa) %>%
  dplyr::summarise(richness = n()) %>%
  dplyr::left_join(dta_rdh) %>%
  mutate(id = paste0("2_", rep)) %>%
  ungroup() %>%
  dplyr::select(id, richness, D, R, pa)

### 4x4m boxes
L = 4
box4 <- expand.grid(
  x0 = seq(min_x, max_x - L, by = L),
  y0 = seq(min_y, max_y - L, by = L)
) %>%
  as.data.frame() %>%
  mutate(id = paste0(L, "_", 1:n()))

comb4 <- lapply(1:nrow(box4), function(i) {
  print(i)
  box <- box4
  dta %>%
    filter(mid_x > box$x0[i], mid_x < (box$x0[i] + L),
           mid_y > box$y0[i], mid_y < (box$y0[i] + L)) %>%
    mutate(id = box$id[i]) %>%
    group_by(id) %>%
    mutate(check = length(unique(rep)) == 4)
}) %>% bind_rows() %>%
  filter(check == T) %>%
  group_by(id) %>%
  summarise(richness = length(unique(species)), D = mean(D), R = mean(R)) %>%
  mutate(pa = 16)


### 6x6m boxes
L = 6
box6 <- expand.grid(
  x0 = seq(min_x, max_x - L, by = L),
  y0 = seq(min_y, max_y - L, by = L)
) %>%
  as.data.frame() %>%
  mutate(id = paste0(L, "_", 1:n()))

comb6 <- lapply(1:nrow(box6), function(i) {
  print(i)
  box <- box6
  dta %>%
    filter(mid_x > box$x0[i], mid_x < (box$x0[i] + L),
           mid_y > box$y0[i], mid_y < (box$y0[i] + L)) %>%
    mutate(id = box$id[i]) %>%
    group_by(id) %>%
    mutate(check = length(unique(rep)) == 9)
}) %>% bind_rows() %>%
  filter(check == T) %>%
  group_by(id) %>%
  summarise(richness = length(unique(species)), D = mean(D), R = mean(R)) %>%
  mutate(pa = 36)

### 8x8m boxes
L = 8
box8 <- expand.grid(
  x0 = seq(min_x, max_x - L, by = L),
  y0 = seq(min_y, max_y - L, by = L)
) %>%
  as.data.frame() %>%
  mutate(id = paste0(L, "_", 1:n()))

comb8 <- lapply(1:nrow(box8), function(i) {
  print(i)
  box <- box6
  dta %>%
    filter(mid_x > box$x0[i], mid_x < (box$x0[i] + L),
           mid_y > box$y0[i], mid_y < (box$y0[i] + L)) %>%
    mutate(id = box$id[i]) %>%
    group_by(id) %>%
    mutate(check = length(unique(rep)) == 16)
}) %>% bind_rows() %>%
  filter(check == T) %>%
  group_by(id) %>%
  summarise(richness = length(unique(species)), D = mean(D), R = mean(R)) %>%
  mutate(pa = 64)

# combine data
dta_boxes <- bind_rows(comb2, comb4, comb6, comb8)

coraldata <- dta_boxes %>%
  mutate(sa = pa*R) # add surface area

#### Analyze
# linear regression model
mod <- lm(log(richness) ~ R + I(R^2) + log(sa) + D, data = coraldata)
summary(mod)


#### figure 
nd <- expand.grid(R = seq(1, 2.5, 0.05), 
                  D = seq(2, 2.6, 0.01), 
                  sa = 4:100) %>%
  as.data.frame() 

pred <- nd %>%
  mutate(logS = predict(mod, newdata = nd)) %>%
  as.data.frame() 


a <- ggplot(pred[pred$D == 2.5,]) +
  geom_raster(aes(x = (sa), y = R,  fill = exp(logS))) +
  scale_fill_fish(direction = 1, option = "Trimma_lantana") +
  geom_hline(yintercept = 1.1, linetype = 2, size = 1) +
  geom_hline(yintercept = 2.2, linetype = 3, size = 1) +
  theme_classic() +
  guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0,
                               barwidth = unit(1, 'lines'), 
                               barheight = unit(20, 'lines')))  +
  labs(fill = "S", x = "Surface area") +
  theme(text = element_text(size = 16), legend.position = "right")
a
b <- ggplot(pred[pred$D == 2.5 & pred$R ==1.1,]) +
  geom_line(aes(x = (sa), y = exp(logS)), size = 1, linetype = 2) +
  scale_fill_fish(direction = 1, option = "Trimma_lantana") +
  scale_y_continuous(limits = c(5,100)) +
  theme_classic() +
  labs(fill = "S", x = "Surface area", y = "S", title = "R = 1.1") +
  theme(text = element_text(size = 16))

b

c <- ggplot(pred[pred$D == 2.5 & pred$R == 2.2,]) +
  geom_line(aes(x = (sa), y = exp(logS)), size = 1, linetype = 3) +
  scale_fill_fish(direction = 1, option = "Trimma_lantana") +
  scale_y_continuous(limits = c(5,100)) +
  theme_classic() +
  labs(fill = "S", x = "Surface area", y = "S", title = "R = 2.2") +
  theme(text = element_text(size = 16))

c


layout <- 
  "AAB
   AAC"

a + c + b + patchwork::plot_layout(design = layout) +
  plot_annotation(tag_levels = "A" ) & theme(text = element_text(size = 18))

ggsave("output/figure1_Ra_richness_fixedD.png", width = 12, height = 8)






