# Load packages ####

library(amt)
library(tidyverse)
library(terra)
library(lubridate)

# Set seed ####

set.seed(22)

# Load data ####

mig_spring_amt <- readRDS("output/mig_spring_100rsteps_with-dyn-covs_pranges-dist_2024-01-08.rds")
mig_fall_amt <- readRDS("output/mig_fall_100rsteps_with-dyn-covs_pranges-dist.rds")

# Load model ####

issa_spring <- readRDS("output/iSSA_spring_2024-01-08.rds")
issa_fall <- readRDS("output/iSSA_fall_2023-12-07.rds")

# Recreate static raster stack ####

elev <- rast("output/processed_layers/elevation_utm.tiff")
cliffs <- rast("output/processed_layers/cliffs.tif")
dist_to_roads <- rast("output/processed_layers/distance_to_roads_crop.tif")
road_poly <- rast("output/processed_layers/road_polygons.tif")
land_cover <- rast("output/processed_layers/land_cover_simple_utm.tiff") %>%
  project(elev)
summ_rast <- rast("output/processed_layers/summer-residencies-by-capture-unit.tiff")
wint_rast <- rast("output/processed_layers/winter-residencies-by-capture-unit.tiff")
summ_dist <- distance(summ_rast)
wint_dist <- distance(wint_rast)

rasts <- rast(list(land_cover, elev, cliffs, dist_to_roads, road_poly))
names(rasts) <- c("land_cover", "elev", "cliffs", "dist_to_roads", "road_poly")

# NDVI ####

ndvi <- rast("output/processed_layers/NDVI.tiff")

# Loads means and sds to scale and center ####

means <- readRDS("output/means.rds")
sds <- readRDS("output/sds.rds")

# Simulate from individual models ####

source("FUN_scale-and-center.R")
source("FUN_dyn.R")

who <- 1

ind <- issa_spring$deploy_ID[who]

# Format model output
issa_spring <- issa_spring %>%
  mutate(coef = map(issa, function(x) broom::tidy(x$model)))

# Get individual coefficients
coefs <- issa_spring %>%
  select(-issa) %>%
  unnest(cols = coef) %>%
  filter(deploy_ID == ind)

coefs_vec <- coefs$estimate
names(coefs_vec) <- coefs$term

## Simulate ####

# First generate a redistribution kernel

start <- make_start(mig_spring_amt$trk[[who]],
                    ta_ = 0,
                    dt = hours(2)) # Starting location

k1_ind <- redistribution_kernel(x = issa_spring$issa[[who]],
                            map = rasts,
                            start = start,
                            fun = function (xy, map) {
                              xy %>%
                                extract_covariates(map, where = "both") %>%
                                attach_dyn(dyn = ndvi, col_name = "ndvi") %>%
                                mutate(log_sl_ = log(sl_),
                                       cos_ta_ = cos(ta_),
                                       dist_to_roads_end_log = log(dist_to_roads_end + 0.001),
                                       ndvi_start = case_when(
                                         land_cover_start == "Open Water" ~ 0,
                                         TRUE ~ ndvi_start
                                       ),
                                       ndvi_end = case_when(
                                         land_cover_end == "Open Water" ~ 0,
                                         TRUE ~ ndvi_end
                                       )) %>%
                                scale_data(means = means, sds = sds)
                            })

# How many tracks to simulate?
n <- 5

# Simulate
sim <- lapply(1:n, function(i){
  cat("Iteration", i, "of", n, "           \n")
  return(simulate_path(k1_ind, n = nrow(mig_spring_amt$trk[[who]])))
  })

## Plot ####

# Unlist results
sim <- bind_rows(sim, .id = "iter")

# Get real data for n steps
dat <- mig_spring_amt$trk[[who]]

# xs
xs <- c(sim$x_, dat$x_)
#ys
ys <- c(sim$y_, dat$y_)

plot(elev,
     xlim = c(min(xs) - 2000, max(xs) + 2000),
     ylim = c(min(ys) - 2000, max(ys) + 2000))
for (i in unique(sim$iter)) {
  lines(sim[sim$iter == i, ]$x_, sim[sim$iter == i, ]$y_)
  points(sim[sim$iter == i, ]$x_, sim[sim$iter == i, ]$y_, pch = i)
}
for (y in unique(year(dat$t_))) {
  lines(dat[year(dat$t_) == y, ]$x_,
        dat[year(dat$t_) == y, ]$y_,
        col = "red", lwd = 2)
  points(dat[year(dat$t_) == y, ]$x_,
         dat[year(dat$t_) == y, ]$y_,
         col = "red", pch = 16)
}
