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
