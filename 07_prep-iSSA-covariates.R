# Prepare covariate layers for iSSA

# I fit the HMM on data subsampled at a 12h resolution.
# However, for the iSSA, I want to fit models to the 2-hour data.
# Load the 2h data and subset migrations by intersecting the dates
# from the HMM output.

# Load packages ####

library(tidyverse)
library(terra)
library(patchwork)
library(amt)
library(lubridate)
library(sf)

# Set seed ####

set.seed(22)

# Load data ####

mule2 <- readRDS("output/mule-deer_regularized-2h_with-NAs.rds")

# Load HMM output ####

mule12 <- readRDS("output/mule-deer_regularized-12h_CTMM_with-viterbi.rds")

# Get migration intervals ####

mule_mig_bouts <- mule12 %>%
  select(ID, t_, vit_3states) %>%
  arrange(ID, t_) %>%
  mutate(start_and_end = case_when(
    # Single-location bouts
    # Within-individual
    (lead(vit_3states) != 3 &
      lag(vit_3states) != 3 &
      vit_3states == 3) |
      # At the beginning
      (lag(ID) != ID &
         lead(vit_3states) != 3 &
         vit_3states == 3) |
      # At the end
      (lead(ID) != ID &
         lag(vit_3states) != 3 &
         vit_3states == 3) ~ "single",
    # First location of one ID
    lag(ID) != ID & vit_3states == 3 ~ "start",
    # First location of a new migration bout
    lag(ID) == ID &
      lag(vit_3states) != 3 & vit_3states == 3 ~ "start",
    # Last location of a migration bout
    lag(ID) == ID & lead(ID) == ID &
      lead(vit_3states) != 3 & vit_3states == 3 ~ "end",
    # Last location of one ID
    lead(ID) != ID & vit_3states == 3 ~ "end",
    # Else, NA
    TRUE ~ NA_character_
  ))

mule_mig_starts <- mule_mig_bouts %>%
  filter(start_and_end == "start") %>%
  select(ID, start_ = t_)

mule_mig_ends <- mule_mig_bouts %>%
  filter(start_and_end == "end") %>%
  select(ID, end_ = t_)

mig_intervals <- mule_mig_starts %>%
  select(-ID) %>%
  bind_cols(mule_mig_ends) %>%
  select(ID, start_, end_)

table(month(mig_intervals$start_))
# April-May and December are peak migrations
# Could restrict analysis to just these
# Or could do it for all/most
# Split:
# Spring migration = April-June
# Fall migration = October-December

spring_mig <- mig_intervals %>%
  filter(month(start_) %in% c(4:6)) %>%
  mutate(mig_int = start_ %--% end_)

fall_mig <- mig_intervals %>%
  filter(month(start_) %in% c(10:12)) %>%
  mutate(mig_int = start_ %--% end_)

# Filter data within migration intervals ####

# Spring

# Test
# i <- spring_mig$ID[352]

spring_mig_data <- data.frame()

for (i in unique(spring_mig$ID)) {

  ints <- spring_mig %>%
    filter(ID == i)

  depID <- word(i, 1, 2, "_")

  all_res <- data.frame()

  for (r in 1:nrow(ints)) {

    res <- mule2 %>%
      filter(deploy_ID == depID &
               t_ %within% ints$mig_int[r])

    all_res <- bind_rows(all_res, res)

  }

  spring_mig_data <- bind_rows(spring_mig_data, all_res)

}

# Fall

fall_mig_data <- data.frame()

for (i in unique(fall_mig$ID)) {

  ints <- fall_mig %>%
    filter(ID == i)

  depID <- word(i, 1, 2, "_")

  all_res <- data.frame()

  for (r in 1:nrow(ints)) {

    res <- mule2 %>%
      filter(deploy_ID == depID &
               t_ %within% ints$mig_int[r])

    all_res <- bind_rows(all_res, res)

  }

  fall_mig_data <- bind_rows(fall_mig_data, all_res)

}

# Format for use in amt ####

mig_spring_amt <- spring_mig_data %>%
  filter(!is.na(x_) & !is.na(y_)) %>%
  nest(cols = -deploy_ID) %>%
  mutate(trk = lapply(cols, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      make_track(.x = x_, .y = y_, .t = t_) %>%
      track_resample(rate = hours(2))
  })) %>%
  mutate(steps = lapply(trk, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      steps_by_burst(keep_cols = "end")
  }))

# Check that all steps are at 2h resolution (FALSE means they all are)
any(unlist(lapply(mig_spring_amt$steps, function(x) {unique(x$dt_)})) != 2)

mig_fall_amt <- fall_mig_data %>%
  filter(!is.na(x_) & !is.na(y_)) %>%
  nest(cols = -deploy_ID) %>%
  mutate(trk = lapply(cols, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      make_track(.x = x_, .y = y_, .t = t_) %>%
      track_resample(rate = hours(2))
  })) %>%
  mutate(steps = lapply(trk, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      steps_by_burst(keep_cols = "end")
  }))

# Check that all steps are at 2h resolution (FALSE means they all are)
any(unlist(lapply(mig_fall_amt$steps, function(x) {unique(x$dt_)})) != 2)

# Remove IDs with few steps and steps with no turn angle ####

mig_spring_amt <- mig_spring_amt %>%
  mutate(steps = lapply(steps, function(x) {
    x %>%
      filter(!is.na(ta_))
  }),
  n_steps = map_dbl(steps, nrow)) %>%
  filter(n_steps >= 20)

mig_fall_amt <- mig_fall_amt %>%
  mutate(steps = lapply(steps, function(x) {
    x %>%
      filter(!is.na(ta_))
  }),
  n_steps = map_dbl(steps, nrow)) %>%
  filter(n_steps >= 20)

# Fit steps and angle distributions ####

# Steps

# Spring
sl_spring <- mig_spring_amt %>%
  dplyr::select(deploy_ID, steps) %>%
  unnest(cols = steps) %>%
  filter(!is.na(sl_) & sl_ > 0) %>%
  pull(sl_)

plot(density(sl_spring))

sl_distr_spring <- fit_distr(sl_spring, dist_name = "gamma")

# Fall
sl_fall <- mig_fall_amt %>%
  dplyr::select(deploy_ID, steps) %>%
  unnest(cols = steps) %>%
  filter(!is.na(sl_) & sl_ > 0) %>%
  pull(sl_)

plot(density(sl_fall))

sl_distr_fall <- fit_distr(sl_fall, dist_name = "gamma")

# Angles

# Spring
ta_spring <- mig_spring_amt %>%
  dplyr::select(deploy_ID, steps) %>%
  unnest(cols = steps) %>%
  filter(!is.na(ta_)) %>%
  pull(ta_)

plot(density(ta_spring))

ta_distr_spring <- fit_distr(ta_spring, dist_name = "vonmises")

# Fall
ta_fall <- mig_fall_amt %>%
  dplyr::select(deploy_ID, steps) %>%
  unnest(cols = steps) %>%
  filter(!is.na(ta_)) %>%
  pull(ta_)

plot(density(ta_fall))

ta_distr_fall <- fit_distr(ta_fall, dist_name = "vonmises")

# Sample random steps ####

mig_spring_amt <- mig_spring_amt %>%
  mutate(rsteps = lapply(steps, FUN = function (x) {
    x %>% random_steps(n_control = 100,
                       sl_distr = sl_distr_spring,
                       ta_distr = ta_distr_spring)
  }))

mig_fall_amt <- mig_fall_amt %>%
  mutate(rsteps = lapply(steps, FUN = function (x) {
    x %>% random_steps(n_control = 100,
                       sl_distr = sl_distr_fall,
                       ta_distr = ta_distr_fall)
  }))


saveRDS(mig_spring_amt, "output/mig_spring_100rsteps.rds")
saveRDS(mig_fall_amt, "output/mig_fall_100rsteps.rds")

# Create extent ####

# Compute minimum and maximum coordinate values

minxs <- min(c(
  unlist(lapply(mig_spring_amt$rsteps, FUN = function(x) {
  xs <- c(x$x1_, x$x2_)
  minx <- min(xs)
  return(minx)
})),
unlist(lapply(mig_fall_amt$rsteps, FUN = function(x) {
  xs <- c(x$x1_, x$x2_)
  minx <- min(xs)
  return(minx)
}))
))

maxxs <- max(c(
  unlist(lapply(mig_spring_amt$rsteps, FUN = function(x) {
    xs <- c(x$x1_, x$x2_)
    maxx <- max(xs)
    return(maxx)
  })),
  unlist(lapply(mig_fall_amt$rsteps, FUN = function(x) {
    xs <- c(x$x1_, x$x2_)
    maxx <- max(xs)
    return(maxx)
  }))
))

minys <- min(c(
  unlist(lapply(mig_spring_amt$rsteps, FUN = function(x) {
    ys <- c(x$y1_, x$y2_)
    miny <- min(ys)
    return(miny)
  })),
  unlist(lapply(mig_fall_amt$rsteps, FUN = function(x) {
    ys <- c(x$y1_, x$y2_)
    miny <- min(ys)
    return(miny)
  }))
))

maxys <- max(c(
  unlist(lapply(mig_spring_amt$rsteps, FUN = function(x) {
    ys <- c(x$y1_, x$y2_)
    maxy <- max(ys)
    return(maxy)
  })),
  unlist(lapply(mig_fall_amt$rsteps, FUN = function(x) {
    ys <- c(x$y1_, x$y2_)
    maxy <- max(ys)
    return(maxy)
  }))
))

ext_mig <- ext(c(minxs, maxxs, minys, maxys))
sa <- as.polygons(ext_mig, crs = "epsg:32612")

# Take a look ####

# mig_spring_amt %>%
#   unnest(rsteps) %>%
#   ggplot() +
#   geom_path(aes(x = x2_, y = y2_, group = deploy_ID, color = deploy_ID)) +
#   facet_wrap(~ case_) +
#   theme(legend.position = "none")

# Looks good but takes forever to plot

# Load and process static covariates ####

## Land cover ####

# Load EVT
evt_layer <- rast("input/LF2020_EVT_220_CONUS/Tif/LC20_EVT_220.tif")

# Crop to larger area to avoid tilt when reprojecting
larger_sa <- rast(extent = ext(sa) + 200000,
                  resolution = res(evt_layer),
                  crs = crs(sa))
larger_sa <- larger_sa %>%
  project(crs(evt_layer)) %>%
  ext()

evt_layer <- crop(evt_layer, larger_sa)

# Load attribute table
evt <- read.csv("input/LF20_EVT_220.csv")

# Filter only categories that exist within the study area
evt_df <- as.data.frame(evt_layer)

sa_evt <- unique(evt_df$EVT_NAME)

evt %>%
  filter(EVT_NAME %in% sa_evt) %>%
  pull(EVT_PHYS) %>%
  unique() %>%
  sort()

new_classif <- evt %>%
  mutate(land_cover = case_when(
    grepl("Blackbrush", EVT_NAME) ~ "Blackbrush",
    grepl("Sagebrush", EVT_NAME) ~ "Sagebrush",
    EVT_PHYS %in% c("Conifer",
                    "Conifer-Hardwood",
                    "Exotic Tree-Shrub",
                    "Hardwood") ~ "Forest",
    EVT_PHYS %in% c("Grassland", "Exotic Herbaceous") ~ "Grassland",
    grepl("Developed", EVT_PHYS) | grepl("Quarries", EVT_PHYS) ~ "Developed",
    # CHECK THIS GROUP BELOW. I MADE AN EXECUTIVE DECISION TO GROUP THESE TWO
    EVT_PHYS %in% c("Snow-Ice", "Sparsely Vegetated") ~ "Bare",
    TRUE ~ EVT_PHYS
  )) %>%
  select(EVT_NAME, land_cover)

evt_df <- evt_df %>%
  left_join(new_classif)

values(evt_layer) <- evt_df$land_cover

evt_layer <- evt_layer %>%
  project(crs(sa)) %>%
  crop(sa)

writeRaster(evt_layer, "output/processed_layers/land_cover_utm.tiff",
            overwrite = T)

evt_freq <- as.data.frame(evt_layer) %>%
  pull(EVT_NAME) %>%
  table() %>%
  as.data.frame()

names(evt_freq) <- c("Cover", "Freq")

p_bar <- ggplot(evt_freq, aes(x = Cover, y = Freq, fill = Cover)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(y = "Frequency", fill = "Land Cover") +
  scale_fill_manual(values = rev(terrain.colors(11))[2:11])

p_map <- ggplot(as.data.frame(evt_layer, xy = TRUE),
                aes(x = x, y = y, fill = EVT_NAME)) +
  geom_raster() +
  scale_fill_manual(values = rev(terrain.colors(11))[2:11]) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_blank())

p_map | p_bar

ggsave("output/land-cover-map_full-extent.tiff", width = 12, height = 6,
       dpi = 300, compression = "lzw")

## Elevation ####

elev <- rast("input/LF2020_Elev_220_CONUS/Tif/LC20_Elev_220.tif") %>%
  crop(larger_sa) %>%
  project(crs(sa)) %>%
  crop(sa)

writeRaster(elev, "output/processed_layers/elevation_utm.tiff", overwrite = T)

## Cliffs ####

# Code from Kezia, extract_cliffs.R

# Extract slopes from DEM
slope_rast <- terra::terrain(elev,
                             v = "slope",
                             unit = "degrees",
                             neighbors = 8)

# Skip this part because it creates problem to have different resolutions
# # Aggregate from 30x30 resolution to 120x120 (factor = 4)
# slope_raster_aggregate <- aggregate(slope_rast, fact = 4, fun = mean)

m <- c(0, 30, 0,
       30, 180, 1)

rclmat <- matrix(m, ncol = 3, byrow = TRUE)

cliffs <- classify(slope_rast,
                                rclmat,
                                include.lowest = TRUE)

writeRaster(cliffs, "output/processed_layers/cliffs.tif",
            overwrite = TRUE)

## Roads ####

# Need a slightly different extent than the actual study area to avoid
# edge effect when calculating distance to roads

larger_sa <- ext(sa) + 100000

# Load roads for UT, CO, WY
roads_ut <- read_sf("input/TIGER Roads/UT/tl_2022_49_prisecroads.shp") %>%
  st_transform(st_crs(32612)) %>%
  st_crop(larger_sa) %>%
  st_union()
roads_co <- read_sf("input/TIGER Roads/CO/tl_2022_08_prisecroads.shp") %>%
  st_transform(st_crs(32612)) %>%
  st_crop(larger_sa) %>%
  st_union()
roads_wy <- read_sf("input/TIGER Roads/WY/tl_2022_56_prisecroads.shp") %>%
  st_transform(st_crs(32612)) %>%
  st_crop(larger_sa) %>%
  st_union()

# Merge into one object
roads <- st_union(roads_ut, roads_co) %>%
  st_union(roads_wy)

st_write(roads, "output/processed_layers/roads.shp")

# Rasterize (to calculate distance in QGIS)
roads_vect <- vect(roads)
roads_rast <- rasterize(roads_vect, cliffs, field = 1)

writeRaster(roads_rast, "output/processed_layers/roads_raster.tif",
            overwrite = TRUE)

# After calculating distance in QGIS, load layer and crop it to study area again
dist_to_roads <- rast("output/processed_layers/distance_to_roads.tif") %>%
  project(elev)

writeRaster(dist_to_roads, "output/processed_layers/distance_to_roads_crop.tif",
            overwrite = TRUE)

## Road-delimited polygons ####

# To see whether a step crossed a road or not, I am going to create a layer
# of polygons delimited by roads. If the starting polygon and the end polygon
# are different, the step crossed a road. If they are the same, the step did not
# cross a road.

# Get roads object and crop to the desired extent
rdp <- roads %>%
  st_crop(elev)

# Get the bounding box as a vector object
bbox <- st_as_sfc(st_bbox(rdp)) %>%
  st_cast("MULTILINESTRING")

# Combine the roads and their bounding box
rdp <- st_union(rdp, bbox)

# Cast the object as a multilinestring
rdp_mls <- st_cast(rdp, "MULTILINESTRING")

# Find intersection points between different lines
rdp_int <- st_intersection(rdp_mls)

# Use the intersection points as vertices of polygons
poly <- rdp_int %>%
  st_union() %>%
  st_polygonize() %>%
  st_collection_extract()

# Transform into terra-compatible object
poly_terra <- vect(poly)

# Assign unique ID to polygons
poly_terra$poly_ID <- 1:nrow(poly_terra)

# Transform to raster
poly_rast <- rasterize(poly_terra, elev, field = "poly_ID")

writeRaster(poly_rast, "output/processed_layers/road_polygons.tif",
            overwrite = TRUE)

# Intersect static covariates ####

elev <- rast("output/processed_layers/elevation_utm.tiff")
cliffs <- rast("output/processed_layers/cliffs.tif")
dist_to_roads <- rast("output/processed_layers/distance_to_roads_crop.tif")
road_poly <- rast("output/processed_layers/road_polygons.tif")
land_cover <- rast("output/processed_layers/land_cover_utm.tiff") %>%
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

# Save dataset ####

saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_with-covs.rds")
saveRDS(mig_fall_amt, "output/mig_fall_100rsteps_with-covs.rds")

# Load and process dynamic covariates ####

## NDVI ####

# Downloaded using script MODIStsp.R

ndvi_files <- list.files("output/processed_layers/MODIS/VI_16Days_250m_v6/NDVI",
                         full.names = TRUE)

ndvi <- rast(ndvi_files)

# plot(ndvi$MOD13A1_NDVI_2016_001)

# test <- ndvi$MOD13A1_NDVI_2016_129
# values(test) <- ifelse(values(ndvi$MOD13A1_NDVI_2016_145) < 0,
#                        NA, values(ndvi$MOD13A1_NDVI_2016_145))
#
# plot(elev, col = rev(grDevices::gray.colors(50)))
# plot(test, alpha = 0.3, add = T)
#
# hist(values(ndvi$MOD13A1_NDVI_2016_145))

# It looks like negative values are invalid (snow or open water).
# I'll replace every pixel that is < 0 with 0.

ndvi_pos <- ndvi
values(ndvi_pos) <- ifelse(values(ndvi_pos) < 0,
                           0, values(ndvi_pos))

# Add z column (timestamp)

# Extract year and julian day from name of layer
y <- stringr::word(names(ndvi_pos), 3, 3, "_")
j <- as.numeric(stringr::word(names(ndvi_pos), 4, 4, "_"))
# Transform julian day to day of year
doy <- stringr::word(as_date(j, origin = "2015-12-31"), 2, 3, "-")
# Add the correct year
z <- ymd(paste0(y, "-", doy))
# Add to SpatRaster
time(ndvi_pos) <- z

ndvi <- ndvi_pos
writeRaster(ndvi, "output/processed_layers/NDVI.tiff", overwrite = TRUE)

# Intersect dynamic covariates ####

source("FUN_ndvi.R")

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
      attach_ndvi(ndvi = ndvi) %>%
      mutate(ndvi_start = case_when(
        land_cover_start == "Open Water" ~ 0,
        TRUE ~ ndvi_start
      ),
      ndvi_end = case_when(
        land_cover_end == "Open Water" ~ 0,
        TRUE ~ ndvi_end
      ))
  }))

saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_with-dyn-covs.rds")
