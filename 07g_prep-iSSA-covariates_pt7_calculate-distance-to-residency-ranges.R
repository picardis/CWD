# Load packages ####

library(terra)

# Load rasters ####

summ_rast <- rast("output/processed_layers/summer-residencies-by-capture-unit.tiff")
wint_rast <- rast("output/processed_layers/winter-residencies-by-capture-unit.tiff")

# Calculate distance ####

summ_dist <- distance(summ_rast)
gc()
writeRaster(summ_dist, "output/processed_layers/distance_to_summer_ranges.tiff",
            overwrite = TRUE, NAflag = NA)
rm(summ_dist)

wint_dist <- distance(wint_rast)
gc()
writeRaster(wint_dist, "output/processed_layers/distance_to_winter_ranges.tiff",
            overwrite = TRUE, NAflag = NA)
