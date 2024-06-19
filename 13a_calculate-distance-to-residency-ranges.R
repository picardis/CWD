# I already calculated distances to summer ranges, but for the purposes of
# parallelizing the simulation, it's better if I can just load one population
# at a time. Create distance rasters and save them as separate files.

library(terra)

summ_rast <- rast("output/processed_layers/summer-residencies-by-capture-unit.tiff")

for (i in 1:nlyr(summ_rast)) {

  focal_summ <- summ_rast[[i]]
  summ_dist <- distance(focal_summ)

  writeRaster(summ_dist,
              paste0("output/processed_layers/dist-to-summer-range_",
                     names(summ_dist),
                     ".tiff"), overwrite = TRUE)

}

wint_rast <- rast("output/processed_layers/winter-residencies-by-capture-unit.tiff")

for (i in 1:nlyr(wint_rast)) {

  focal_wint <- wint_rast[[i]]
  wint_dist <- distance(focal_wint)

  writeRaster(wint_dist,
              paste0("output/processed_layers/dist-to-winter-range_",
                     names(wint_dist),
                     ".tiff"), overwrite = TRUE)

}
