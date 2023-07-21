# Load packages ####

library(amt)
library(tidyverse)
library(terra)

# Set seed ####

set.seed(22)

# Load data ####

mig_spring_amt <- readRDS("output/mig_spring_100rsteps.rds")

# Load model ####

issa_spring <- readRDS("output/iSSA_spring_2023-07-20.rds")

# Recreate static raster stack ####

elev <- rast("output/processed_layers/elevation_utm.tiff")
cliffs <- rast("output/processed_layers/cliffs.tif")
dist_to_roads <- rast("output/processed_layers/distance_to_roads_crop.tif")
road_poly <- rast("output/processed_layers/road_polygons.tif")
land_cover <- rast("output/processed_layers/land_cover_utm.tiff") %>%
  project(elev)

rasts <- rast(list(land_cover, elev, cliffs, dist_to_roads, road_poly))
names(rasts) <- c("land_cover", "elev", "cliffs", "dist_to_roads", "road_poly")

# NDVI ####

ndvi <- rast("output/processed_layers/NDVI.tiff")

# Loads means and sds to scale and center ####

means <- readRDS("output/means.rds")
sds <- readRDS("output/sds.rds")

# Simulate from average model coefficients ####

## Average model parameters ####

# Simulate from a model that uses the mean value of each beta across individuals
# IT MIGHT BE BETTER TO SIMULATE FROM THE INDIVIDUAL MODELS INSTEAD AND THEN
# AVERAGE THE SIMULATIONS AFTER THE FACT.

# Format model output
issa_spring <- issa_spring %>%
  mutate(coef = map(issa, function(x) broom::tidy(x$model)))

# Unnest coefficients
coefs <- issa_spring %>%
  select(-issa) %>%
  unnest(cols = coef) %>%
  group_by(term) %>%
  summarize(mean_beta = mean(estimate))

coefs_vec <- coefs$mean_beta
names(coefs_vec) <- coefs$term
# For some reason, when I fit the redistribution kernel below, the name of this
# coefficient gets swapped from A interacts B to B interacts A. Then it breaks
# everything because it can't match the coefficient names. Fix it manually here
# and later I'll submit an issue on GitHub.
names(coefs_vec)[6] <- "elev_start_sc:I(elev_end_sc - elev_start_sc)"

## Create simulation model ####

sim_mod_spring <- make_issf_model(coefs = coefs_vec,
                                  sl = issa_spring$issa[[1]]$sl_,
                                  ta = issa_spring$issa[[1]]$ta_)

## Simulate ####

source("FUN_scale-and-center.R")

# First generate a redistribution kernel

start <- make_start(mig_spring_amt$trk[[1]],
                    ta_ = 0,
                    dt = hours(2)) # Starting location

k1 <- redistribution_kernel(x = sim_mod_spring,
                            map = rasts,
                            start = start,
                            fun = function (xy, map) {
                              xy %>%
                                extract_covariates(map, where = "both") %>%
                                mutate(log_sl_ = log(sl_),
                                       cos_ta_ = cos(ta_),
                                       dist_to_roads_end_log = log(dist_to_roads_end + 0.001)) %>%
                                scale_data(means = means, sds = sds)
                            })

# How many tracks to simulate?
n <- 10

# Simulate
sim <- replicate(n, simulate_path(k1, n = 30), simplify = FALSE)

## Plot ####

# Layer
elev_df <- as.data.frame(elev, xy = TRUE)

# Unlist results
sim <- bind_rows(sim, .id = "iter")

# Get real data for n steps
dat <- mig_spring_amt$trk[[1]][1:30, ]

plot(elev,
     xlim = c(min(sim$x_) - 5000, max(sim$x_) + 5000),
     ylim = c(min(sim$y_) - 5000, max(sim$y_) + 5000))
for (i in unique(sim$iter)) {
  lines(sim[sim$iter == i, ]$x_, sim[sim$iter == i, ]$y_)
  points(sim[sim$iter == i, ]$x_, sim[sim$iter == i, ]$y_, pch = i)
}
lines(dat$x_, dat$y_, col = "red", lwd = 2)
points(dat$x_, dat$y_, col = "red", pch = 16)

# Simulate from individual models ####

source("FUN_scale-and-center.R")
source("FUN_NDVI.R")

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
                                attach_ndvi(ndvi = ndvi) %>%
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
lines(dat$x_, dat$y_, col = "red", lwd = 2)
points(dat$x_, dat$y_, col = "red", pch = 16)
