###############################################################################
### Troubleshooting ###########################################################
## In case of error, figure out which rows of the data failed
# Load deploy_IDs from object
# inds <- readRDS("output/iSSA_spring_2024-01-08.rds")[[1]]
# # Load files from output
# op <- list.files("output/simulations/")
# # Get individuals from output
# op_ind_1 <- sapply(strsplit(op, split = "_", fixed = TRUE),
#                    getElement, 2)
# op_ind_2 <- sapply(strsplit(op, split = "_", fixed = TRUE),
#                    getElement, 3)
#
# op_inds <- paste(op_ind_1, op_ind_2, sep = "_")
# # Find out which inds are missing from output
# (missing <- inds[which(!(inds %in% op_inds))])
# # Indices (should be row numbers)
# (missing_inds <- which(!(inds %in% op_inds)))
# # Or find a specific one
# x <- which(inds == "MD19F0176_45035")
###############################################################################X

# To run in batches:
# # Select only individuals that have not been done yet
# done <- stringr::word(list.files("output/simulations"), 2, 3, "_")
# issa_spring <- readRDS("output/iSSA_spring_2024-06-19.rds")
# sub <- issa_spring[!issa_spring$deploy_ID %in% done, ]
# saveRDS(sub, "output/iSSA_spring_2024-06-19_SUBSET_2024-06-19.rds")
# # Restart R!

# Load packages ####

library(parallel)

# Parallelize ####

# Uses ~ 5GB RAM/cluster. Limit to 10 clusters for RAM considerations.
clust <- makeCluster(10)

#system.time({parLapply(clust, 1:107, function(x) {
  system.time({parLapply(clust, 1:2, function(x) {

  library(amt)
  library(tidyverse)
  library(terra)
  library(lubridate)

  # Set seed ####

  set.seed(22 + x)

  # Source functions ####

  source("FUN_scale-and-center.R")
  source("FUN_dyn.R")
  source("ssf_formula_fix2.R")

  # Load data ####

  mig_spring_amt <- readRDS("output/mig_spring_100rsteps_with-dyn-covs_pranges-dist_2024-06-19.rds")

  # Load model ####

  issa_spring <- readRDS("output/iSSA_spring_2024-06-19.rds") # Replace with subset if running in batches

  # Recreate static raster stack ####
  lyrs <- c("output/processed_layers/elevation_utm.tiff",
            "output/processed_layers/cliffs.tif",
            "output/processed_layers/distance_to_roads_crop.tif",
            "output/processed_layers/road_polygons.tif",
            "output/processed_layers/land_cover_simple_utm.tiff")

  rasts <- rast(lyrs)

  names(rasts) <- c("elev", "cliffs", "dist_to_roads", "road_poly", "land_cover")

  # NDVI ####

  ndvi <- rast("output/processed_layers/NDVI.tiff")

  # Loads means and sds to scale and center ####

  means <- readRDS("output/means_2024-06-19.rds")
  sds <- readRDS("output/sds_2024-06-19.rds")

  ind <- issa_spring$deploy_ID[x]

  # Residency ranges ####

  csu <- mig_spring_amt %>%
    select(deploy_ID, csu) %>%
    filter(deploy_ID == ind) %>%
    pull(csu) %>%
    unique() %>%
    unlist()

  summ_dist <- rast(paste0(
    "output/processed_layers/dist-to-summer-range_",
    csu,
    ".tiff"))

  names(summ_dist) <- "dist_to_range"

  # Garbage cleanup after processing all the rasters
  gc()

  # Simulate from individual models ####

  ## Simulate ####

  # First generate a redistribution kernel

  start <- make_start(mig_spring_amt$trk[[which(mig_spring_amt$deploy_ID == ind)]],
                      ta_ = 0,
                      dt = hours(2)) # Starting location

  # Store coefficients and set any NA coefficients to 0
  B <- coef(issa_spring$issa[[which(issa_spring$deploy_ID == ind)]])
  B[is.na(B)] <- 0

  # Make iSSF model with substituted coefficients
  mm <- make_issf_model(coefs = B,
                        sl = sl_distr(issa_spring$issa[[which(issa_spring$deploy_ID == ind)]]),
                        ta = ta_distr(issa_spring$issa[[which(issa_spring$deploy_ID == ind)]]))

  k1_ind <- redistribution_kernel(x = mm,
                                  map = rasts,
                                  start = start,
                                  fun = function (xy, map) {
                                    xy %>%
                                      extract_covariates(map, where = "both") %>%
                                      extract_covariates(summ_dist, where = "both") %>%
                                      mutate(land_cover_endAgricultural = as.numeric(
                                        land_cover_end == "Agricultural"),
                                        land_cover_endDeveloped = as.numeric(
                                          land_cover_end == "Developed"),
                                        land_cover_endGrassland = as.numeric(
                                          land_cover_end == "Grassland"),
                                        land_cover_endShrubland = as.numeric(
                                          land_cover_end == "Shrubland"),
                                        land_cover_endUnsuitable = as.numeric(
                                          land_cover_end == "Unsuitable")) %>%
                                      attach_dyn(dyn = ndvi, col_name = "ndvi") %>%
                                      mutate(log_sl_ = log(sl_),
                                             cos_ta_ = cos(ta_),
                                             dist_to_roads_end_log = log(dist_to_roads_end + 0.001),
                                             dist_to_range_end_log = log(dist_to_range_end + 0.001),
                                             # land_cover_end = factor(land_cover_end, levels = c("Forest",
                                             #                                                    "Agricultural",
                                             #                                                    "Developed",
                                             #                                                    "Grassland",
                                             #                                                    "Shrubland",
                                             #                                                    "Unsuitable")),
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
  sim <- NULL

  try({
    sim <- lapply(1:n, function(i){
      cat("Iteration", i, "of", n, "           \n")
      return(simulate_path(k1_ind,
                           #n = nrow(mig_spring_amt$trk[[which(mig_spring_amt$deploy_ID == ind)]])
                           n = 360 #max(sapply(mig_spring_amt$trk, nrow))
      ))
    })
  })

  if (is.null(sim)) {
    # Delete anything saved
    unlink(paste0("output/simulations/sim_", ind, "_spring_", Sys.Date(), ".rds"))

  } else {
    saveRDS(sim, paste0("output/simulations/sim_", ind, "_spring_", Sys.Date(), ".rds"))
  }


})
})
