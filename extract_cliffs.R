library(terra)
### i. Bring in DEM ----
# read in elevation raster
#elev_rast <- terra::rast("Data/GISLayers/DEMs/USGS_13_n39w110_20220510.tif")
elev_rast <- terra::rast("input/USGS_13_n39w110_20220510.tif")

### ii. Extract slopes
slope_rast <- terra::terrain(elev_rast, v = "slope",
                             unit = "degrees",
                             neighbors = 8)

plot(slope_rast)

### v. Distance to escape terrain
# classify the values into three blocks
# all values >= 0 and <= 0.25 become 1, etc.
# m <- c(0, 30, 0,
#        30, 180, 1)
# rclmat <- matrix(m, ncol=3, byrow=TRUE)
# cliffs_class <- classify(slope_rast, rclmat, include.lowest=TRUE)
# plot(cliffs_class)

#aggregate from 40x40 resolution to 120x120 (factor = 3)
slope_raster_aggregate <- aggregate(slope_rast, fact = 20, fun=mean)
res(slope_raster_aggregate)
plot(slope_raster_aggregate)
m <- c(0, 30, 0,
       30, 180, 1)
rclmat <- matrix(m, ncol = 3, byrow = TRUE)
cliffs_coarse_class <- classify(slope_raster_aggregate, rclmat, include.lowest=TRUE)
plot(cliffs_coarse_class)

# FOR TESTING
library(sf)
library(stars)
# clif <- crop(cliffs_coarse_class, ext(-109.6, -109.4, 38.4, 38.6))
clif <- cliffs_coarse_class

# clif2 <- focal(clif, w = 3, fun = function(x) {
#   return(as.numeric(sum(x) > 3))
# })

clif2 <- resample(clif)

x <- st_as_stars(clif2) %>%
  st_as_sf(merge = TRUE) %>% # this is the raster to polygons part
  st_cast("MULTILINESTRING")

plot(x)

# Find intersection points between different lines
x_int <- st_intersection(x)

# Use the intersection points as vertices of polygons
poly <- x %>%
  st_union() %>%
  st_polygonize() %>%
  st_collection_extract()

# Transform into terra-compatible object
poly_terra <- vect(poly)

# Assign unique ID to polygons
poly_terra$poly_ID <- 1:nrow(poly_terra)

# Transform to raster
poly_rast <- rasterize(poly_terra, clif2, field = "poly_ID")

plot(poly_rast)
