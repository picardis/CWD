# Load packages ####

library(tidyverse)
library(terra)

# Load data ####

summ_res <- readRDS("output/residency_data_summer.rds")
wint_res <- readRDS("output/residency_data_winter.rds")

elev <- rast("output/processed_layers/elevation_utm.tiff")

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

rm(list=ls()[! ls() %in% c("wint_res", "elev")])

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
