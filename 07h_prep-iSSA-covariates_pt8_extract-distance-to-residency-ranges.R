# Load packages ####

library(tidyverse)
library(terra)
library(amt)

# Load data ####

mig_spring_amt <- readRDS("output/mig_spring_100rsteps_with-dyn-covs_pranges-sample_2024-06-19.rds")
mig_fall_amt <- readRDS("output/mig_fall_100rsteps_with-dyn-covs_pranges-sample_2024-06-19.rds")

# Extract distance to summer ranges ####
summ_dist <- rast("output/processed_layers/distance_to_summer_ranges.tiff")

mig_spring_amt <- mig_spring_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function (x) {
    cap <- x$csu[1]
    r <- summ_dist[[cap]]
    names(r) <- "dist_to_summ"
    x %>%
      extract_covariates(r, where = "both")
  }))

saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_with-dyn-covs_pranges-dist_2024-06-19.rds")

rm(mig_spring_amt, summ_dist)

# Extract distance to winter ranges
wint_dist <- rast("output/processed_layers/distance_to_winter_ranges.tiff")

mig_fall_amt <- mig_fall_amt %>%
  filter(!csu == "East Canyon") %>% # No residency data available for this CSU
  mutate(rsteps = lapply(rsteps, FUN = function (x) {
    cap <- x$csu[1]
    r <- wint_dist[[cap]]
    names(r) <- "dist_to_wint"
    x %>%
      extract_covariates(r, where = "both")
  }))

saveRDS(mig_fall_amt, "output/mig_fall_100rsteps_with-dyn-covs_pranges-dist_2024-06-19.rds")
