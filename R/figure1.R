library(tidyverse)
library(fishualize)
library(patchwork)
library(readr)
library(brms)

##### load data #####
# data source: Torres-Pulliza et al. (2020)
# coral species per 2x2m box
dta_coral <- read_csv("data/coralspecies.csv")
# Complexity measures per 2x2m box
dta_rdh <- read_csv("data/megaplot.csv") %>%
  dplyr::select(rep, R, D)
# Combine datasets
dta <- left_join(dta_coral, dta_rdh)

# bounding box of coral data
min_x <- min(dta_coral$mid_x) - 1
min_y <- min(dta_coral$mid_y) - 1
max_x <- max(dta_coral$mid_x) + 1
max_y <- max(dta_coral$mid_y) + 1

##### Species diversity for boxes of varying sizes #####
### 2x2m boxes 
comb2 <- dta_coral %>%
  dplyr::group_by(rep, pa) %>%
  dplyr::summarise(richness = n(), abundance = sum(abu)) %>%
  dplyr::left_join(dta_rdh) %>%
  mutate(id = paste0("2_", rep)) %>%
  ungroup() %>%
  mutate(Rsd = NA, Dsd = NA) %>%
  dplyr::select(id, richness, abundance, Dm = D, Rm = R, Dsd, Rsd, pa) 

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
  summarise(richness = length(unique(species)), 
            abundance = sum(abu), 
            Dm = mean(unique(D)), Rm = mean(unique(R)), 
            Dsd = sd(unique(D)), Rsd = sd(unique(R))) %>%
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
  summarise(richness = length(unique(species)), 
            abundance = sum(abu), 
            Dm = mean(unique(D)), Rm = mean(unique(R)), 
            Dsd = sd(unique(D)), Rsd = sd(unique(R))) %>%
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
  summarise(richness = length(unique(species)), 
            abundance = sum(abu), 
            Dm = mean(unique(D)), Rm = mean(unique(R)), 
            Dsd = sd(unique(D)), Rsd = sd(unique(R))) %>%
  mutate(pa = 64)

# combine data
dta_boxes <- bind_rows(comb2, comb4, comb6, comb8)

coraldata <- dta_boxes %>%
  mutate(sa = pa*Rm, R2 = Rm^2) # add surface area

##### Analyze #####

# Bayesian linear regression model
mod <- brm(log(richness) ~ Rm + R2 + log(sa) + Dm, 
            data = coraldata,
            backend = "cmdstanr", cores = 4)

# posterior predictive check
pp_check(mod)
# model summary
summary(mod)
# bayesian R2
bayes_R2(mod)
# Extract priors
prior_summary(mod)

##### figure #####

# new data frame to predict to
nd <- expand.grid(Rm = seq(1, 3, 0.05), 
                  Dm = 2.5, 
                  sa = 4:100) %>%
  as.data.frame() %>%
  mutate(R2 = Rm^2)

# predict
pred <- nd %>%
  mutate(logS = fitted(mod, newdata = nd)[,1]) %>%
  as.data.frame() 

# plots
a <- ggplot(pred) +
  geom_raster(aes(x = (sa), y = Rm,  fill = exp(logS))) +
  scale_fill_fish(direction = 1, option = "Trimma_lantana") +
  geom_hline(yintercept = 1.1, linetype = 2, size = 1) +
  geom_hline(yintercept = 2.2, linetype = 3, size = 1) +
  theme_classic() +
  guides(fill = guide_colorbar(title.position = 'top', title.hjust = 0,
                               barwidth = unit(1, 'lines'), 
                               barheight = unit(20, 'lines')))  +
  labs(fill = "S", x = "Surface area", y = "R") +
  theme(text = element_text(size = 16), legend.position = "right")
a
b <- ggplot(pred[pred$Rm ==1.1,]) +
  geom_line(aes(x = (sa), y = exp(logS)), size = 1, linetype = 2) +
  scale_fill_fish(direction = 1, option = "Trimma_lantana") +
  scale_y_continuous(limits = c(5,100)) +
  theme_classic() +
  labs(fill = "S", x = "Surface area", y = "S", title = "R = 1.1") +
  theme(text = element_text(size = 16))

b

c <- ggplot(pred[pred$Rm == 2.2,]) +
  geom_line(aes(x = (sa), y = exp(logS)), size = 1, linetype = 3) +
  scale_fill_fish(direction = 1, option = "Trimma_lantana") +
  scale_y_continuous(limits = c(5,105)) +
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






