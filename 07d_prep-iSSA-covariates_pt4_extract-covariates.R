# Load packages ####

library(tidyverse)
library(terra)
library(amt)
library(lubridate)
library(sf)

# Load data ####

mig_spring_amt <- readRDS("output/mig_spring_100rsteps_2024-06-19.rds")
mig_fall_amt <- readRDS("output/mig_fall_100rsteps_2024-06-19.rds")

# Intersect static covariates ####

elev <- rast("output/processed_layers/elevation_utm.tiff")
cliffs <- rast("output/processed_layers/cliffs.tif")
dist_to_roads <- rast("output/processed_layers/distance_to_roads_crop.tif")
road_poly <- rast("output/processed_layers/road_polygons.tif")
land_cover <- rast("output/processed_layers/land_cover_simple_utm.tiff") %>%
  project(elev)

rasts <- rast(list(land_cover, elev, cliffs, dist_to_roads, road_poly))
names(rasts) <- c("land_cover", "elev", "cliffs", "dist_to_roads", "road_poly")

mig_spring_amt <- mig_spring_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function (x) {
    x %>%
      extract_covariates(rasts, where = "both")
  }))

mig_fall_amt <- mig_fall_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function (x) {
    x %>%
      extract_covariates(rasts, where = "both")
  }))

# Intersect dynamic covariates ####

ndvi <- rast("output/processed_layers/NDVI.tiff")
snow <- rast("output/processed_layers/snow.tiff")

source("FUN_dyn.R")

mig_spring_amt <- mig_spring_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function (x) {
    x %>%
      # There is something weird happening with extract_covariates_var_time.
      # It attaches NA values when the values are not NA.
      # See troubleshooting example in FUN_ndvi.R
      # Which is why I'm using a custom function instead
      # extract_covariates_var_time(ndvi,
      #                             where = "both",
      #                             when = "before",
      #                             max_time = days(17)) %>%
      # rename(ndvi_start = time_var_covar_start,
      #        ndvi_end = time_var_covar_end)
      attach_dyn(dyn = ndvi, col_name = "ndvi") %>%
      # The same NDVI function works for snow as well
      attach_dyn(dyn = snow, col_name = "snow") %>%
      mutate(ndvi_start = case_when(
        land_cover_start == "Open Water" ~ 0,
        TRUE ~ ndvi_start
      ),
      ndvi_end = case_when(
        land_cover_end == "Open Water" ~ 0,
        TRUE ~ ndvi_end
      ))
  }))

mig_fall_amt <- mig_fall_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function (x) {
    x %>%
      attach_dyn(dyn = ndvi, col_name = "ndvi") %>%
      attach_dyn(dyn = snow, col_name = "snow") %>%
      mutate(ndvi_start = case_when(
        land_cover_start == "Open Water" ~ 0,
        TRUE ~ ndvi_start
      ),
      ndvi_end = case_when(
        land_cover_end == "Open Water" ~ 0,
        TRUE ~ ndvi_end
      ))
  }))

saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_with-dyn-covs_2024-06-19.rds")
saveRDS(mig_fall_amt, "output/mig_fall_100rsteps_with-dyn-covs_2024-06-19.rds")
