# Load packages ####

library(tidyverse)
library(terra)
library(amt)
library(lubridate)
library(sf)

# Load study area polygon ####

sa <- vect("output/study_area_polygon.shp")

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
