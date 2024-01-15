# Load packages ####

library(amt)
library(tidyverse)
library(terra)
library(lubridate)

# Set seed ####

set.seed(22)

# Load data ####

mig_spring_amt <- readRDS("output/mig_spring_100rsteps_with-dyn-covs_pranges-dist_2024-01-08.rds")
#mig_fall_amt <- readRDS("output/mig_fall_100rsteps_with-dyn-covs_pranges-dist.rds")

# Load model ####

issa_spring <- readRDS("output/iSSA_spring_2024-01-08.rds")
#issa_fall <- readRDS("output/iSSA_fall_2023-12-07.rds")

# Recreate static raster stack ####

elev <- rast("output/processed_layers/elevation_utm.tiff")
cliffs <- rast("output/processed_layers/cliffs.tif")
dist_to_roads <- rast("output/processed_layers/distance_to_roads_crop.tif")
road_poly <- rast("output/processed_layers/road_polygons.tif")
land_cover <- rast("output/processed_layers/land_cover_simple_utm.tiff") %>%
  project(elev)

rasts <- rast(list(land_cover, elev, cliffs, dist_to_roads, road_poly))
names(rasts) <- c("land_cover", "elev", "cliffs", "dist_to_roads", "road_poly")

# NDVI ####

ndvi <- rast("output/processed_layers/NDVI.tiff")

# Choose individual ####

who <- 9
# 4

ind <- issa_spring$deploy_ID[who]

# Residency ranges ####

csu <- mig_spring_amt %>%
  select(deploy_ID, csu) %>%
  filter(deploy_ID == ind) %>%
  pull(csu) %>%
  unique() %>%
  unlist()

summ_rast <- rast("output/processed_layers/summer-residencies-by-capture-unit.tiff")
#wint_rast <- rast("output/processed_layers/winter-residencies-by-capture-unit.tiff")

focal_summ <- summ_rast[[which(names(summ_rast) == csu)]]
#focal_wint <- wint_rast[[which(names(wint_rast) == csu)]]

summ_dist <- distance(focal_summ)
#wint_dist <- distance(focal_wint)

names(summ_dist) <- "dist_to_range"
#names(wint_dist) <- "dist_to_range"

rasts_spring <- c(rasts, summ_dist)
#rasts_fall <- c(rasts, wint_dist)

# Loads means and sds to scale and center ####

means <- readRDS("output/means_2023-12-07.rds")
sds <- readRDS("output/sds_2023-12-07.rds")

# Simulate from individual models ####

source("FUN_scale-and-center.R")
source("FUN_dyn.R")

# ## Re-assign step length and turn angle distributions ####
#
# issa_spring$issa[[which(issa_spring$deploy_ID == ind)]]$sl_ <- attr(mig_spring_amt$rsteps[[1]], "sl_")
# issa_spring$issa[[which(issa_spring$deploy_ID == ind)]]$ta_ <- attr(mig_spring_amt$rsteps[[1]], "ta_")

## Simulate ####

source("ssf_formula_fix2.R")

# First generate a redistribution kernel

start <- make_start(mig_spring_amt$trk[[which(mig_spring_amt$deploy_ID == ind)]],
                    ta_ = 0,
                    dt = hours(2)) # Starting location

k1_ind <- redistribution_kernel(x = issa_spring$issa[[which(issa_spring$deploy_ID == ind)]],
                                map = rasts_spring,
                                start = start,
                                fun = function (xy, map) {
                                  xy %>%
                                    extract_covariates(map, where = "both") %>%
                                    attach_dyn(dyn = ndvi, col_name = "ndvi") %>%
                                    mutate(log_sl_ = log(sl_),
                                           cos_ta_ = cos(ta_),
                                           dist_to_roads_end_log = log(dist_to_roads_end + 0.001),
                                           dist_to_range_end_log = log(dist_to_range_end + 0.001),
                                           land_cover_end = factor(land_cover_end, levels = c("Forest",
                                                                                              "Agricultural",
                                                                                              "Developed",
                                                                                              "Grassland",
                                                                                              "Shrubland",
                                                                                              "Unsuitable")),
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
n <- 100

# Simulate
sim <- lapply(1:n, function(i){
  cat("Iteration", i, "of", n, "           \n")
  return(simulate_path(k1_ind, n = nrow(mig_spring_amt$trk[[which(mig_spring_amt$deploy_ID == ind)]])))
})

## Plot ####

# Unlist results
sim <- bind_rows(sim, .id = "iter")

# Get real data for n steps
dat <- mig_spring_amt$trk[[which(mig_spring_amt$deploy_ID == ind)]]

# xs
xs <- c(sim$x_, dat$x_)
#ys
ys <- c(sim$y_, dat$y_)

library(ggplot2)

minx <- min(dat$x_ - 10000)
maxx <- max(dat$x_ + 10000)
miny <- min(dat$y_ - 10000)
maxy <- max(dat$y_ + 10000)

r1 <- as.data.frame(summ_dist, xy = TRUE)
r2 <- as.data.frame(elev, xy = TRUE)

p1 <- ggplot(dat, aes(x_, y_)) +
  geom_raster(mapping = aes(x, y, fill = LC20_Elev_220),
              data = r2) +
  geom_path(color = "springgreen2") +
  geom_point(color = "springgreen2") +
  scale_fill_viridis_c(option = "cividis") +
  theme_bw() +
  labs(x = "UTM Easting", y = "UTM Northing", fill = "Elevation (m)") +
  coord_cartesian(xlim = c(minx, maxx), ylim = c(miny, maxy))

ggsave(paste0("output/real-track_", ind, "_",
              csu, "_spring_elev_", n, "-sims.tiff"),
       plot = p1, width = 6, height = 4, dpi = 300, compression = "lzw")

p2 <- ggplot(sim, aes(x_, y_)) +
  geom_raster(mapping = aes(x, y, fill = LC20_Elev_220),
              data = r2) +
  geom_path(aes(group = iter), alpha = 0.3, color = "white") +
  geom_point(aes(group = iter), alpha = 0.2, color = "white") +
  scale_fill_viridis_c(option = "cividis") +
  theme_bw() +
  labs(x = "UTM Easting", y = "UTM Northing", fill = "Elevation (m)") +
  coord_cartesian(xlim = c(minx, maxx), ylim = c(miny, maxy))

ggsave(paste0("output/sim-tracks_", ind, "_",
              csu, "_spring_elev_", n, "-sims.tiff"),
       plot = p2, width = 6, height = 4, dpi = 300, compression = "lzw")

p3 <- ggplot(sim, aes(x_, y_)) +
  geom_raster(mapping = aes(x, y, fill = LC20_Elev_220),
              data = r2) +
  geom_path(aes(group = iter), alpha = 0.1, color = "white") +
  geom_point(aes(group = iter), alpha = 0.1, color = "white") +
  geom_path(data = dat, color = "springgreen2", size = 1) +
  geom_point(data = dat, color = "springgreen2", size = 2) +
  scale_fill_viridis_c(option = "cividis") +
  theme_bw() +
  labs(x = "UTM Easting", y = "UTM Northing", fill = "Elevation (m)") +
  coord_cartesian(xlim = c(minx, maxx), ylim = c(miny, maxy))

ggsave(paste0("output/real-and-sim-tracks_", ind, "_",
              csu, "_spring_elev_", n, "-sims.tiff"),
       plot = p3, width = 6, height = 4, dpi = 300, compression = "lzw")

p4 <- ggplot(sim, aes(x_, y_)) +
  geom_path(aes(group = iter)) +
  geom_point(aes(group = iter)) +
  facet_wrap(~ iter, nrow = 2) +
  theme_bw() +
  labs(x = "UTM Easting", y = "UTM Northing", fill = "Elevation (m)") +
  coord_cartesian(xlim = c(minx + 8000, maxx - 8000),
                  ylim = c(miny + 8000, maxy - 8000))

ggsave(paste0("output/sim-tracks_", ind, "_",
              csu, "_spring_no-bg_", n, "-sims.tiff"),
       plot = p4, width = 8, height = 6, dpi = 300, compression = "lzw")

p5 <- ggplot(sim, aes(x_, y_)) +
  geom_path(data = dat, color = "springgreen2", size = 1) +
  geom_point(data = dat, color = "springgreen2", size = 2) +
  geom_path(aes(group = iter)) +
  geom_point(aes(group = iter)) +
  facet_wrap(~ iter, nrow = 2) +
  theme_bw() +
  labs(x = "UTM Easting", y = "UTM Northing", fill = "Elevation (m)") +
  coord_cartesian(xlim = c(minx + 8000, maxx - 8000),
                  ylim = c(miny + 8000, maxy - 8000))

ggsave(paste0("output/real-and-sim-tracks_", ind, "_",
              csu, "_spring_no-bg_", n, "-sims.tiff"),
       plot = p5, width = 8, height = 6, dpi = 300, compression = "lzw")
