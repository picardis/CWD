# Make a map of migration corridors based on the state assignments from the HMM

# Load packages ####

library(momentuHMM)
library(raster)
library(tidyverse)
library(patchwork)

# Load HMM ####

hmm2 <- readRDS("output/HMM2_12h_2023-04-25_NAs.rds")

# Get Viterbi assignments ####

vit <- viterbi(hmm2)

# Load Utah boundary ####

# Load UT boundary
utah <- sf::read_sf("../../../Other People/Ben Crabb/From Ben/clean_data/boundaries/Utah_NAD83_12N.shp")

# Create template raster
template <- raster(ext = extent(utah), res = 1000, crs = crs(utah))

# Separate states ####

mig <- hmm2$data %>%
  cbind(vit) %>%
  rename(state = vit) %>%
  filter(state == 2)

res <- hmm2$data %>%
  cbind(vit) %>%
  rename(state = vit) %>%
  filter(state == 1)

# Rasterize points ####

mig_raster <- rasterize(mig[, c("x", "y")],
                        template, field = 1, fun = "sum")

res_raster <- rasterize(res[, c("x", "y")],
                        template, field = 1, fun = "sum")

res_yn_raster <- res_raster
values(res_yn_raster) <- ifelse(is.na(values(res_raster)),
                                  NA, 1)

mig_only_raster <- mig_raster
values(mig_only_raster) <- ifelse(is.na(values(res_yn_raster)),
                                  values(mig_raster), NA)

mor_df <- as.data.frame(mig_only_raster, xy = TRUE)
res_df <- as.data.frame(res_raster, xy = TRUE)
mig_df <- as.data.frame(mig_raster, xy = TRUE)
res_yn_df <- as.data.frame(res_yn_raster, xy = TRUE)

p_mor <- ggplot() +
  geom_raster(data = mor_df,
              mapping = aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(na.value = NA) +
  theme_void()

p_mig <- ggplot() +
  geom_raster(data = mig_df,
              mapping = aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(na.value = NA) +
  theme_void()

p_res <- ggplot() +
  geom_raster(data = res_df,
              mapping = aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(na.value = NA) +
  theme_void()

p_res_yn <- ggplot() +
  geom_raster(data = res_yn_df,
              mapping = aes(x = x, y = y, fill = layer)) +
  theme_void()

p_res + p_res_yn + p_mig + p_mor + plot_layout(guides = "collect")

ggplot() +
  geom_sf(data = utah) +
  geom_raster(data = mor_df,
              mapping = aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(na.value = NA) +
  theme_void()

# Brownian bridges ####
