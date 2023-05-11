# Make a map of migration corridors based on the state assignments from the HMM

# Load packages ####

library(momentuHMM)
library(raster)
library(tidyverse)
library(patchwork)
library(sf)

# Load data ####

mule <- readRDS("output/mule-deer_regularized-12h_CTMM_with-viterbi.rds")

# Load Utah boundary ####

# Load UT boundary
utah <- sf::read_sf("../../../Other People/Ben Crabb/From Ben/clean_data/boundaries/Utah_NAD83_12N.shp")

# Create template raster
template <- raster(ext = extent(utah), res = 1000, crs = crs(utah))

# Separate states ####

mig <- mule %>%
  filter(vit_3states == 3) %>%
  as.data.frame()

res <- mule %>%
  filter(vit_3states %in% c(1:2)) %>%
  as.data.frame()
#
# # Rasterize points ####
#
# mig_raster <- rasterize(mig[, c("x", "y")],
#                         template, field = 1, fun = "sum")
#
# res_raster <- rasterize(res[, c("x", "y")],
#                         template, field = 1, fun = "sum")
#
# res_yn_raster <- res_raster
# values(res_yn_raster) <- ifelse(is.na(values(res_raster)),
#                                   NA, 1)
#
# mig_only_raster <- mig_raster
# values(mig_only_raster) <- ifelse(is.na(values(res_yn_raster)),
#                                   values(mig_raster), NA)
#
# mig_res_raster <- res_yn_raster
# # 1 is migration, 0 is residence, NA is not used
# values(mig_res_raster) <- case_when(
#   values(res_yn_raster) == 1 ~ 0,
#   !is.na(values(mig_only_raster)) ~ 1,
#   TRUE ~ NA_real_
# )
#
# mor_df <- as.data.frame(mig_only_raster, xy = TRUE)
# res_df <- as.data.frame(res_raster, xy = TRUE)
# mig_df <- as.data.frame(mig_raster, xy = TRUE)
# res_yn_df <- as.data.frame(res_yn_raster, xy = TRUE)
# mig_res_df <- as.data.frame(mig_res_raster, xy = TRUE)
#
# p_mor <- ggplot() +
#   geom_raster(data = mor_df,
#               mapping = aes(x = x, y = y, fill = layer)) +
#   scale_fill_viridis_c(na.value = NA) +
#   theme_void()
#
# p_mig <- ggplot() +
#   geom_raster(data = mig_df,
#               mapping = aes(x = x, y = y, fill = layer)) +
#   scale_fill_viridis_c(na.value = NA) +
#   theme_void()
#
# p_res <- ggplot() +
#   geom_raster(data = res_df,
#               mapping = aes(x = x, y = y, fill = layer)) +
#   scale_fill_viridis_c(na.value = NA) +
#   theme_void()
#
# p_res_yn <- ggplot() +
#   geom_raster(data = res_yn_df,
#               mapping = aes(x = x, y = y, fill = factor(layer))) +
#   scale_fill_discrete(na.value = NA) +
#   theme_void()
#
# p_res | p_res_yn | p_mig | p_mor | plot_layout(guides = "collect")
#
# ggplot() +
#   geom_sf(data = utah) +
#   geom_raster(data = mig_res_df,
#               mapping = aes(x = x, y = y, fill = factor(layer))) +
#   scale_fill_viridis_d(na.value = NA) +
#   theme_void()

# Lines ####

## Migration steps ####
# Separate out each migratory burst before connecting the lines

mig_no_nas <- mig %>%
  filter(!is.na(x))

mig_no_nas$burst <- NA
mig_no_nas$burst[1] <- 1

for (i in 2:nrow(mig_no_nas)) {

  cat(i, "of", nrow(mig_no_nas), "\n")

  gap <- as.numeric(mig_no_nas$t_[i] - mig_no_nas$t_[i - 1])
  mig_no_nas$burst[i] <- ifelse(gap == 12,
                                mig_no_nas$burst[i - 1],
                                mig_no_nas$burst[i - 1] + 1)

}

# Get rid of bursts with only one point
toss <- mig_no_nas %>%
  group_by(burst) %>%
  tally() %>%
  filter(n == 1) %>%
  pull(burst)

mig_no_nas <- mig_no_nas %>%
  filter(!burst %in% toss)

mig_lines <- mig_no_nas %>%
  st_as_sf(coords = c("x", "y"), crs = 32612) %>%
  group_by(burst) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

ggplot() +
  geom_sf(data = utah) +
  geom_sf(data = mig_lines, aes(color = burst)) +
  scale_color_viridis_c(na.value = NA) +
  theme_void() +
  theme(legend.position = "none")

mig_lines_rast <- rasterize(mig_lines, template, fun = "count")

mig_lines_df <- as.data.frame(mig_lines_rast, xy = TRUE)

ggplot() +
  geom_sf(data = utah) +
  geom_raster(data = mig_lines_df,
              mapping = aes(x = x, y = y, fill = layer_burst)) +
  scale_fill_viridis_c(na.value = NA) +
  theme_void()

## Residency steps ####
# Separate out each migratory burst before connecting the lines

res_no_nas <- res %>%
  filter(!is.na(x))

res_no_nas$burst <- NA
res_no_nas$burst[1] <- 1

for (i in 2:nrow(res_no_nas)) {

  cat(i, "of", nrow(res_no_nas), "\n")

  gap <- as.numeric(res_no_nas$t_[i] - res_no_nas$t_[i - 1])
  res_no_nas$burst[i] <- ifelse(gap == 12,
                                res_no_nas$burst[i - 1],
                                res_no_nas$burst[i - 1] + 1)

}

# Get rid of bursts with only one point
toss <- res_no_nas %>%
  group_by(burst) %>%
  tally() %>%
  filter(n == 1) %>%
  pull(burst)

res_no_nas <- res_no_nas %>%
  filter(!burst %in% toss)

res_lines <- res_no_nas %>%
  st_as_sf(coords = c("x", "y"), crs = 32612) %>%
  group_by(burst) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

ggplot() +
  geom_sf(data = utah) +
  geom_sf(data = res_lines, aes(color = burst)) +
  scale_color_viridis_c(na.value = NA) +
  theme_void() +
  theme(legend.position = "none")

res_lines_rast <- rasterize(res_lines, template, fun = "count")

res_lines_df <- as.data.frame(res_lines_rast, xy = TRUE)

ggplot() +
  geom_sf(data = utah) +
  geom_raster(data = res_lines_df,
              mapping = aes(x = x, y = y, fill = layer_burst)) +
  scale_fill_viridis_c(na.value = NA) +
  theme_void()

## Binary migration raster ####

mig_bin_ras <- mig_lines_rast
values(mig_bin_ras) <- ifelse(is.na(values(mig_lines_rast)), NA, 1)

## Binary residency raster ####

res_bin_ras <- res_lines_rast
values(res_bin_ras) <- ifelse(is.na(values(res_lines_rast)), NA, 1)

## Combined raster ####

bin_ras <- res_bin_ras
# 1 is migration, 0 is residence, NA is not used
values(bin_ras) <- ifelse(values(bin_ras) == 1, 0, NA)
values(bin_ras) <- ifelse(values(mig_bin_ras) == 1 &
                            is.na(values(bin_ras)), 1, values(bin_ras))

bin_ras <- setValues(raster(bin_ras), bin_ras[])

bin_ras_df <- as.data.frame(bin_ras, xy = TRUE)

ggplot() +
  geom_sf(data = utah) +
  geom_raster(data = bin_ras_df,
              mapping = aes(x = x, y = y, fill = factor(layer))) +
  scale_fill_viridis_d(na.value = NA) +
  theme_void() +
  labs(fill = "Migration")

ggsave("output/map_mig-vs-res_1km_HMM3_2023-05-10.tiff",
       height = 8, width = 7, dpi = 300, compression = "lzw")

writeRaster(bin_ras, "output/mig-vs-res_1km_HMM3_2023-05-10.tiff")
