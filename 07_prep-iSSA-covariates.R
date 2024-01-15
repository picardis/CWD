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

mule2 <- readRDS("output/mule-deer_regularized-2h_with-NAs_2023-07-22.rds")

# Load HMM output ####

mule12 <- readRDS("output/mule-deer_regularized-12h_CTMM_with-viterbi.rds")

# Load new information on capture sites ####

cap_lasal <- read.csv("input/lasals_animals_with_cap_areas_krmmanual.csv") %>%
  mutate(deploy_ID = paste0(uniqueID, "_", collarID)) %>%
  select(deploy_ID, realCaptur)
cap_sjuan <- read.csv("input/sanjuan_animals_with_cap_areas_krmmanual.csv") %>%
  mutate(deploy_ID = paste0(uniqueID, "_", collarID)) %>%
  select(deploy_ID, realCaptur)

cap_all <- bind_rows(cap_lasal, cap_sjuan)

mule2 <- mule2 %>%
  left_join(cap_all)

mule12 <- mule12 %>%
  left_join(cap_all)

# Now create new column that uses realCaptur when available and the old
# capture subunit when realCaptur is not available
mule2 <- mule2 %>%
  mutate(csu = case_when(
    !is.na(realCaptur) ~ realCaptur,
    is.na(realCaptur) ~ captureSubUnit
  )) %>%
  mutate(csu = case_when(
    csu == "Abajo Mtns" ~ "Abajo",
    TRUE ~ csu
  ))

mule12 <- mule12 %>%
  mutate(csu = case_when(
    !is.na(realCaptur) ~ realCaptur,
    is.na(realCaptur) ~ captureSubUnit
  )) %>%
  mutate(csu = case_when(
    csu == "Abajo Mtns" ~ "Abajo",
    TRUE ~ csu
  ))

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

# This is one spot where I lose a bunch of individuals. I go from 596 to
# 387 in the spring and 511 in the fall. Others have migration steps in
# different months.

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
      make_track(.x = x_, .y = y_, .t = t_,
                 all_cols = TRUE) %>%
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
      make_track(.x = x_, .y = y_, .t = t_,
                 all_cols = TRUE) %>%
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

saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_2023-11-21.rds")
saveRDS(mig_fall_amt, "output/mig_fall_100rsteps_2023-11-21.rds")

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

# Further simplify land cover
land_cover <- rast("output/processed_layers/land_cover_utm.tiff")

land_cover_simple <- as.data.frame(land_cover) %>%
  mutate(new_lc = case_when(
    EVT_NAME %in% c("Sagebrush", "Shrubland", "Blackbrush") ~ "Shrubland",
    EVT_NAME %in% c("Bare", "Open Water") ~ "Unsuitable",
    EVT_NAME %in% c("Grassland", "Riparian") ~ "Grassland",
    TRUE ~ as.character(EVT_NAME)
  ))

lc_simple <- land_cover
values(lc_simple) <- land_cover_simple$new_lc

unique(land_cover_simple$new_lc)

writeRaster(lc_simple, "output/processed_layers/land_cover_simple_utm.tiff",
            overwrite = T)

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

# Save dataset ####

saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_with-covs.rds")
saveRDS(mig_fall_amt, "output/mig_fall_100rsteps_with-covs.rds")

# Load and process dynamic covariates ####

## NDVI ####

# Downloaded using script MODIStsp.R

ndvi_files <- list.files("output/processed_layers/MODIS/VI_16Days_250m_v6/NDVI",
                         full.names = TRUE)

ndvi <- rast(ndvi_files)

# Negative values are invalid (snow or open water).
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

## Snow ####

# Downloaded using script MODIStsp.R

snow_files <- list.files("output/processed_layers/MODIS/Snow_Cov_8-Day_500m_v6/MAX_SNW",
                         full.names = TRUE)

snow <- rast(snow_files)
# table(values(snow))
#
# 1        25        37        50       100       200
# 637 302024113   4325443  10490864   1474619 102972884

# Values assigned based on the user manual. 1 = snow, 0 = no snow, NA = unknown
values(snow) <- ifelse(values(snow) %in% c(0, 1, 11, 50, 254, 255),
                       NA, ifelse(values(snow) %in% c(25, 37, 39, 100),
                                  0, 1))

# Add z column (timestamp)

# Extract year and julian day from name of layer
y <- stringr::word(names(snow), 5, 5, "_")
j <- as.numeric(stringr::word(names(snow), 6, 6, "_"))
# Transform julian day to day of year
doy <- stringr::word(as_date(j, origin = "2015-12-31"), 2, 3, "-")
# Add the correct year
z <- ymd(paste0(y, "-", doy))
# Add to SpatRaster
time(snow) <- z

writeRaster(snow, "output/processed_layers/snow.tiff", overwrite = TRUE)

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

saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_with-dyn-covs.rds")
saveRDS(mig_fall_amt, "output/mig_fall_100rsteps_with-dyn-covs.rds")

# Population ranges ####

# Use 12h individuals which weren't included in iSSA
pop_data <- mule12 %>%
  filter(!deploy_ID %in% c(mig_spring_amt$deploy_ID,
                           mig_fall_amt$deploy_ID))

# Residency data
summ_res <- pop_data %>%
  filter(month(t_) %in% 7:9 &
           vit_3states %in% 1:2)
wint_res <- pop_data %>%
  filter(month(t_) %in% 1:3 &
           vit_3states %in% 1:2)

## WHY WOULD THIS MATTER?
# # Take only deploy_ID that have data for both winter and summer
all_ids <- unique(c(summ_res$deploy_ID, wint_res$deploy_ID))
# not_in_wint <- setdiff(unique(summ_res$deploy_ID), unique(wint_res$deploy_ID))
# not_in_summ <- setdiff(unique(wint_res$deploy_ID), unique(summ_res$deploy_ID))
# focal_ids <- all_ids[!all_ids %in% c(not_in_wint, not_in_summ)]

n_by_unit_focal <- pop_data %>%
  as.data.frame() %>%
  #filter(deploy_ID %in% focal_ids) %>%
  filter(deploy_ID %in% all_ids & !is.na(csu)) %>%
  dplyr::select(deploy_ID, csu) %>%
  distinct() %>%
  group_by(csu) %>%
  tally() %>%
  rename(n_both_seasons = n)

focal_known_cap <- pop_data %>%
  as.data.frame() %>%
  # replaced focal with all
  filter(deploy_ID %in% all_ids & !is.na(csu)) %>%
  pull(deploy_ID) %>%
  unique()

summ_res <- summ_res %>%
  filter(deploy_ID %in% focal_known_cap)
wint_res <- wint_res %>%
  filter(deploy_ID %in% focal_known_cap)

ggplot(summ_res, aes(x = x, y = y, color = csu)) +
  geom_point()
ggplot(wint_res, aes(x = x, y = y, color = csu)) +
  geom_point()

summ_res %>%
  mutate(season = "Summer") %>%
  bind_rows(wint_res) %>%
  mutate(season = case_when(
    is.na(season) ~ "Winter",
    TRUE ~ season
  )) %>%
  ggplot(aes(x = x, y = y, color = season)) +
  geom_point(alpha = 0.1)

# I need to have residency data for at least one animal in each CSU that's in
# the iSSA data. I need one summer for the spring iSSA CSUs and one winter for
# the fall iSSA CSUs. Check:
need_summer <- unique(unlist(lapply(1:nrow(mig_spring_amt), FUN = function(x) {
  unique(mig_spring_amt$rsteps[[x]]$csu)})))
need_winter <- unique(unlist(lapply(1:nrow(mig_fall_amt), FUN = function(x) {
  unique(mig_fall_amt$rsteps[[x]]$csu)})))
needed_csu <- data.frame(csu = sort(c(unique(need_summer, need_winter)))) %>%
  mutate(needed_summer = csu %in% need_summer,
         needed_winter = csu %in% need_winter) %>%
  # Do I have what I need from the 12h data?
  mutate(got_summer = csu %in% summ_res$csu,
         got_winter = csu %in% wint_res$csu)

# For the CSUs that I don't have info for in the 12h data, I can subset a few
# individuals from the 2h data.
miss_summer <- needed_csu %>%
  filter(needed_summer == TRUE & got_summer == FALSE) %>%
  pull(csu)
miss_winter <- needed_csu %>%
  filter(needed_winter == TRUE & got_winter == FALSE) %>%
  pull(csu)

# Get residency data from 2h data
mig_data_2h <- spring_mig_data %>%
  bind_rows(fall_mig_data) %>%
  mutate(mig = TRUE)

res_data_2h <- mule2 %>%
  left_join(mig_data_2h) %>%
  filter(is.na(mig)) %>%
  dplyr::select(-mig)

# How many individuals do I have for each of those in the 2h data?
n_summ_2h <- res_data_2h %>%
  filter(csu %in% miss_summer,
         month(t_) %in% 7:9) %>%
  select(csu, deploy_ID) %>%
  distinct() %>%
  group_by(csu) %>%
  tally()

n_wint_2h <- res_data_2h %>%
  filter(csu %in% miss_winter,
         month(t_) %in% 1:3) %>%
  select(csu, deploy_ID) %>%
  distinct() %>%
  group_by(csu) %>%
  tally()

# Randomly sample 2 of them per CSU
summ_options <- res_data_2h %>%
  filter(csu %in% miss_summer,
         month(t_) %in% 7:9) %>%
  select(csu, deploy_ID) %>%
  distinct()

summ_samp <- c()
for (u in 1:length(unique(summ_options$csu))) {
  choices <- summ_options[summ_options$csu == unique(summ_options$csu)[u], ]
  if (nrow(choices) > 1) {
  summ_samp <- c(summ_samp,
                 sample(choices$deploy_ID, size = 2, replace = FALSE))
  } else {
    summ_samp[u] <- choices$deploy_ID
    }
}

wint_options <- res_data_2h %>%
  filter(csu %in% miss_winter,
         month(t_) %in% 7:9) %>%
  select(csu, deploy_ID) %>%
  distinct()

wint_samp <- c()
for (u in 1:length(unique(wint_options$csu))) {
  choices <- wint_options[wint_options$csu == unique(wint_options$csu)[u], ]
  if (nrow(choices) > 1) {
    wint_samp <- c(wint_samp,
                   sample(choices$deploy_ID, size = 2, replace = FALSE))
  } else {
    wint_samp[u] <- choices$deploy_ID
  }
}

# Get them
summ_2h <- res_data_2h %>%
  filter(deploy_ID %in% summ_samp,
         month(t_) %in% 7:9)

wint_2h <- res_data_2h %>%
  filter(deploy_ID %in% wint_samp,
         month(t_) %in% 1:3)

# Merge all residency data together now
summ_res <- summ_2h %>%
  rename(x = x_, y = y_) %>%
  select(deploy_ID, x, y, csu) %>%
  bind_rows(summ_res[, c("deploy_ID", "x", "y", "csu")])

wint_res <- wint_2h %>%
  rename(x = x_, y = y_) %>%
  select(deploy_ID, x, y, csu) %>%
  bind_rows(wint_res[, c("deploy_ID", "x", "y", "csu")])

# Make table of how many individuals per CSU
summ_res %>%
  select(deploy_ID, csu) %>%
  distinct() %>%
  group_by(csu) %>%
  tally() %>%
  as.data.frame()

wint_res %>%
  select(deploy_ID, csu) %>%
  distinct() %>%
  group_by(csu) %>%
  tally() %>%
  as.data.frame()

# Rasterize the residency layers by capture unit

summ_rast <- lapply(unique(summ_res$csu), function(x) {
  rasterize(as.matrix(summ_res[summ_res$csu == x,
                               c("x", "y")]),
            elev,
            values = 1)
}) %>%
  rast()
names(summ_rast) <- unique(summ_res$csu)

writeRaster(summ_rast, "output/processed_layers/summer-residencies-by-capture-unit.tiff",
            overwrite = TRUE, NAflag = NA)

gc()

wint_rast <- lapply(unique(wint_res$csu), function(x) {
  rasterize(as.matrix(wint_res[wint_res$csu == x,
                               c("x", "y")]),
            elev,
            values = 1)
}) %>%
  rast()
names(wint_rast) <- unique(wint_res$csu)

writeRaster(wint_rast, "output/processed_layers/winter-residencies-by-capture-unit.tiff",
            overwrite = TRUE, NAflag = NA)

gc()

summ_dist <- distance(summ_rast)
gc()
writeRaster(summ_dist, "output/processed_layers/distance_to_summer_ranges.tiff",
            overwrite = TRUE, NAflag = NA)
rm(summ_dist)

wint_dist <- distance(wint_rast)
gc()
writeRaster(wint_dist, "output/processed_layers/distance_to_winter_ranges.tiff",
            overwrite = TRUE, NAflag = NA)

rm(list=ls()[! ls() %in% c("mig_spring_amt",
                           "mig_fall_amt",
                           "summ_dist",
                           "wint_dist",
                           "summ_2h",
                           "wint_2h")])

gc()

### ALTERNATIVE
summ_dist <- rast("output/processed_layers/distance_to_summer_ranges.tiff")

mig_spring_amt <- mig_spring_amt %>%
  mutate(csu = lapply(rsteps, FUN = function(x) {x$csu[1]}))

mig_spring_amt <- mig_spring_amt %>%
  filter(!deploy_ID %in% summ_2h$deploy_ID &
           !is.na(csu)) %>%
  mutate(rsteps = lapply(rsteps, FUN = function (x) {
    cap <- x$csu[1]
    r <- summ_dist[[cap]]
    names(r) <- "dist_to_summ"
    x %>%
      extract_covariates(r, where = "both")
  }))
### END ALTERNATIVE

# ### OLD
# # Extract distance values
# mig_spring_amt <- mig_spring_amt %>%
#   filter(!deploy_ID %in% summ_2h$deploy_ID) %>%
#   mutate(rsteps = lapply(rsteps, FUN = function (x) {
#     x %>%
#       extract_covariates(summ_dist, where = "both")
#   }))
#
# gc()
#
# mig_fall_amt <- mig_fall_amt %>%
#   filter(!deploy_ID %in% wint_2h$deploy_ID) %>%
#   mutate(rsteps = lapply(rsteps, FUN = function (x) {
#     x %>%
#       extract_covariates(wint_dist, where = "both")
#   }))
#
# gc()
#
# saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_with-pop-ranges.rds")
# saveRDS(mig_fall_amt, "output/mig_fall_100rsteps_with-pop-ranges.rds")
#
# mig_spring_amt$rsteps <- lapply(mig_spring_amt$rsteps, function(x) {
#
#   cap <- unique(x$csu)
#   if(is.na(cap)) {
#     x$dist_to_summ_end <- NA
#   } else {
#   vals <- x[, grepl(paste0(cap, "_end"),
#                     colnames(x))]
#   names(vals) <- "dist_to_summ_end"
#   x <- cbind(x, vals)}
#   return(x)
#
# })
#
# mig_fall_amt$rsteps <- lapply(mig_fall_amt$rsteps, function(x) {
#
#   cap <- unique(x$csu)
#   if(is.na(cap)) {
#     x$dist_to_wint_end <- NA
#   } else {
#     vals <- x[, grepl(paste0(cap, "_end"),
#                       colnames(x))]
#     names(vals) <- "dist_to_wint_end"
#     x <- cbind(x, vals)}
#   return(x)
#
# })
#
# # Get rid of distance to irrelevant CSUs
# mig_spring_amt <- mig_spring_amt %>%
#   mutate(rsteps = lapply(rsteps, function(x) {
#     x %>%
#       select(x1_:snow_end, dist_to_summ_end)
#   }))
#
# mig_fall_amt <- mig_fall_amt %>%
#   mutate(rsteps = lapply(rsteps, function(x) {
#     x %>%
#       select(x1_:snow_end, dist_to_wint_end)
#   }))
# ### END OLD

saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_with-dyn-covs_pranges-dist_2024-01-08.rds")
saveRDS(mig_fall_amt, "output/mig_fall_100rsteps_with-dyn-covs_pranges-dist.rds")
