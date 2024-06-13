# Load packages ####

library(terra)
library(tidyverse)
library(patchwork)

# Load data on real migrations ####

mig_spring_amt <- readRDS("output/mig_spring_100rsteps_with-dyn-covs_pranges-dist_2024-01-08.rds")

mig_spring <- mig_spring_amt %>%
  select(-csu) %>%
  unnest(cols = trk)

rm(mig_spring_amt)

# Table of migration durations
mig_spring %>%
  group_by(burst) %>%
  mutate(track_dur = difftime(max(t_), min(t_), unit = "days"),
         start = min(t_),
         end = max(t_)) %>%
  select(burst, track_dur, start, end) %>%
  distinct()

mig_durations <- mig_spring %>%
  mutate(ind_year = paste0(deploy_ID, "_", year(t_))) %>%
  group_by(ind_year) %>%
  mutate(track_dur = difftime(max(t_), min(t_), unit = "days"),
         start = min(t_),
         end = max(t_)) %>%
  select(ind_year, track_dur, start, end) %>%
  distinct()

saveRDS(mig_durations, "output/spring_mig_durations.rds")

# Load validation data ####

vali <- readRDS("output/validation_data.rds")

vali <- vali %>%
  filter(vit_3states == 3 &
           month(t_) %in% 4:6) %>%
  mutate(csu = case_when(
             !is.na(realCaptur) ~ realCaptur,
             is.na(realCaptur) ~ captureSubUnit
           ))

# Load simulation results ####

# Spring
files <- list.files("output/simulations", pattern = "spring", full.names = TRUE)

ind_sims <- data.frame()

for (f in 1:length(files)) {

  print(f)

  ind_sim <- readRDS(files[f])

  for (s in 1:length(ind_sim)) {

    ind_sim[[s]]$sim <- s

  }

  ind_sim <- dplyr::bind_rows(ind_sim)

  ind_sim$id <- stringr::word(files[f], 2, 3, "_")

  ind_sims <- dplyr::bind_rows(ind_sims, ind_sim)

}

# Load template raster ####

template <- rast("output/processed_layers/elevation_utm.tiff")

# Subset one region at a time? ####

# If I calculate the percent points per pixel, imbalances in different efforts
# across region are going to complicate the interpretation of results
# Might be better to calculate percentages relative to one region (one or more
# CSUs)

# csus <- unique(mig_spring[mig_spring$x_ >= 630000 &
#              mig_spring$x_ <= 680000 &
#              mig_spring$y_ >= 4222000 &
#              mig_spring$y_ <= 4290000, ]$csu)
#
# mig_spring <- mig_spring %>%
#   filter(csu %in% csus) %>%
#   nrow()
#
# vali <- vali %>%
#   filter(csu %in% csus) %>%
#   nrow()
#
# ind_sims <- ind_sims %>%
#   filter(id %in% unique(mig_spring$deploy_ID))

# Or crop to an extent and calculate everything relative to that extent
# r <- rast()
# ext(r) <- c(560000, 700000,
#             4130000, 4300000)
# template <- crop(template, r)

# Make points ####

## Simulated points ####

pts <- vect(ind_sims, geom = c("x_", "y_"), crs = crs(template))

## Real points ####

real_pts <- vect(mig_spring[, c("x_", "y_")],
                 geom = c("x_", "y_"), crs = crs(template))

## Validation points ####

vali_pts <- vect(as.data.frame(vali[, c("x", "y")]),
                 geom = c("x", "y"), crs = crs(template))

# Rasterize points ####

## Simulated points ####

sim_ras <- rasterize(x = pts, y = template, fun = "count")

# Log of the raw values
sim_ras_log <- sim_ras
values(sim_ras_log) <- log(values(sim_ras))
values(sim_ras_log) <- ifelse(is.na(values(sim_ras_log)), 0, values(sim_ras_log))

plot(sim_ras_log, col = viridis::viridis(n = 10))

# Aggregated log values
sim_ras_log_agg2 <- aggregate(sim_ras_log, fact = 100)
plot(sim_ras_log_agg2, col = viridis::viridis(n = 10))

# Standardized log values
mean_srl <- mean(values(sim_ras_log_agg2), na.rm = TRUE)
sd_srl <- sd(values(sim_ras_log_agg2), na.rm = TRUE)
sim_ras_log_sc <- (sim_ras_log_agg2 - mean_srl)/sd_srl
plot(sim_ras_log_sc, col = viridis::viridis(n = 10))

# Percent values
sim_ras_perc <- sim_ras/nrow(ind_sims)

# Standardized log of percent values
sim_ras_perc_log <- log(sim_ras_perc)
mean_srpl <- mean(values(sim_ras_perc_log), na.rm = TRUE)
sd_srpl <- sd(values(sim_ras_perc_log), na.rm = TRUE)
sim_ras_perc_log_sc <- (sim_ras_perc_log - mean_srpl)/sd_srpl
plot(sim_ras_perc_log_sc, col = viridis::viridis(n = 10))

# THIS IS BAD: FIRST SCALE AND CENTER AND THEN TAKE LOG? THROWS OUT HALF THE VALUES
# sim_ras_perc <- sim_ras/nrow(ind_sims) # Maybe need to divide by number of points?
# mean_srp <- mean(values(sim_ras_perc), na.rm = TRUE)
# sd_srp <- sd(values(sim_ras_perc), na.rm = TRUE)
# sim_ras_perc <- (sim_ras_perc - mean_srp)/sd_srp
# sim_ras_perc_log <- log(sim_ras_perc)
# values(sim_ras_perc_log) <- ifelse(is.na(values(sim_ras_perc_log)),
#                                    min(values(sim_ras_perc_log), na.rm = TRUE) - 1,
#                                    values(sim_ras_perc_log))
# sim_ras_perc_log_agg <- aggregate(sim_ras_perc_log, fact = 100)
# plot(sim_ras_perc_log_agg, col = viridis::viridis(n = 10))

tiff("output/sim_spring.tiff", width = 4, height = 8, units = "in",
     res = 300, compression = "lzw")
plot(sim_ras_log_sc, col = viridis::viridis(n = 10))
dev.off()

## Real points ####

real_ras <- rasterize(x = real_pts, y = template, fun = "count")

# Log of the raw values
real_ras_log <- real_ras
values(real_ras_log) <- log(values(real_ras))
values(real_ras_log) <- ifelse(is.na(values(real_ras_log)), 0, values(real_ras_log))

plot(real_ras_log, col = viridis::viridis(n = 10))

# Aggregated log values
real_ras_log_agg2 <- aggregate(real_ras_log, fact = 100)
plot(real_ras_log_agg2, col = viridis::viridis(n = 10))

# Standardized log values
mean_srl <- mean(values(real_ras_log_agg2), na.rm = TRUE)
sd_srl <- sd(values(real_ras_log_agg2), na.rm = TRUE)
real_ras_log_sc <- (real_ras_log_agg2 - mean_srl)/sd_srl
plot(real_ras_log_sc, col = viridis::viridis(n = 10))

# Percent values
real_ras_perc <- real_ras/nrow(mig_spring)

# Standardized log of percent values
real_ras_perc_log <- log(real_ras_perc)
mean_srpl <- mean(values(real_ras_perc_log), na.rm = TRUE)
sd_srpl <- sd(values(real_ras_perc_log), na.rm = TRUE)
real_ras_perc_log_sc <- (real_ras_perc_log - mean_srpl)/sd_srpl
plot(real_ras_perc_log_sc, col = viridis::viridis(n = 10))

# THIS IS BAD: FIRST SCALE AND CENTER AND THEN TAKE LOG? THROWS OUT HALF THE VALUES
# real_ras_perc <- real_ras/nrow(mig_spring)
# mean_srp <- mean(values(real_ras_perc), na.rm = TRUE)
# sd_srp <- sd(values(real_ras_perc), na.rm = TRUE)
# real_ras_perc <- (real_ras_perc - mean_srp)/sd_srp
# real_ras_perc_log <- log(real_ras_perc)
# values(real_ras_perc_log) <- ifelse(is.na(values(real_ras_perc_log)),
#                                     min(values(real_ras_perc_log), na.rm = TRUE) - 1,
#                                     values(real_ras_perc_log))
# real_ras_perc_log_agg <- aggregate(real_ras_perc_log, fact = 100)
# plot(real_ras_perc_log_agg, col = viridis::viridis(n = 10))

tiff("output/real_spring.tiff", width = 4, height = 8, units = "in",
     res = 300, compression = "lzw")
plot(real_ras_perc_log_agg, col = viridis::viridis(n = 10))
dev.off()

## Validation points ####

vali_ras <- rasterize(x = vali_pts, y = template, fun = "count")

vali_ras_log <- vali_ras
values(vali_ras_log) <- log(values(vali_ras))
values(vali_ras_log) <- ifelse(is.na(values(vali_ras_log)), 0, values(vali_ras_log))

plot(vali_ras_log, col = viridis::viridis(n = 10))

vali_ras_log_agg2 <- aggregate(vali_ras_log, fact = 100)
plot(vali_ras_log_agg2, col = viridis::viridis(n = 10))

vali_ras_perc <- vali_ras/nrow(vali)
mean_srp <- mean(values(vali_ras_perc), na.rm = TRUE)
sd_srp <- sd(values(vali_ras_perc), na.rm = TRUE)
vali_ras_perc <- (vali_ras_perc - mean_srp)/sd_srp
vali_ras_perc_log <- log(vali_ras_perc)
values(vali_ras_perc_log) <- ifelse(is.na(values(vali_ras_perc_log)),
                                    min(values(vali_ras_perc_log), na.rm = TRUE) - 1,
                                    values(vali_ras_perc_log))
vali_ras_perc_log_agg <- aggregate(vali_ras_perc_log, fact = 100)
plot(vali_ras_perc_log_agg, col = viridis::viridis(n = 10))

tiff("output/vali_spring.tiff", width = 4, height = 8, units = "in",
     res = 300, compression = "lzw")
plot(vali_ras_perc_log_agg, col = viridis::viridis(n = 10))
dev.off()

tiff("output/real_sim_vali_spring.tiff", width = 13, height = 8, units = "in",
     res = 300, compression = "lzw")
par(mfrow = c(1, 3))
plot(real_ras_perc_log_agg, col = viridis::viridis(n = 10), main = "Observed")
plot(sim_ras_perc_log_agg, col = viridis::viridis(n = 10), main = "Simulated")
plot(vali_ras_perc_log_agg, col = viridis::viridis(n = 10), main = "Validation")
dev.off()

# Plot cropped maps ####

## Load UT boundary
utah <- sf::read_sf("../../../../Other People/Ben Crabb/From Ben/clean_data/boundaries/Utah_NAD83_12N.shp")

## Example 1 ####
r <- rast()
ext(r) <- c(560000, 700000,
           4130000, 4300000)
values(r) <- 1

# Who is in this example?

ex1 <- as.data.frame(pts) %>%
  cbind.data.frame(geom(pts)) %>%
  filter(x >= 560000 & x <= 700000 &
           y >= 4130000 & y <= 4300000) %>%
  group_by(id) %>%
  tally() %>%
  arrange(desc(n))

ex1$id

# Map inset
map <- ggplot() +
  geom_sf(data = utah) +
  geom_raster(data = as.data.frame(r, xy = TRUE),
              mapping = aes(x = x, y = y), fill = "tomato") +
  theme_void() +
  theme(legend.position = "none")

crop1 <- crop(sim_ras_log_agg2, r)
crop1_df <- as.data.frame(crop1, xy = TRUE)

ggplot() +
  geom_raster(data = crop1_df, aes(x = x, y = y, fill = count)) +
  geom_path(data = mig_spring, aes(x = x_, y = y_, group = burst),
            color = "white", alpha = 0.1) +
  geom_point(data = mig_spring, aes(x = x_, y = y_),
            color = "white", alpha = 0.1) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(630000, 680000),
                  ylim = c(4222000, 4290000)) +
  theme_void() +
  theme(legend.position = "none")

crop2sim <- crop(sim_ras_perc_log_agg, r)
crop2sim_df <- as.data.frame(crop2sim, xy = TRUE)

crop2real <- crop(real_ras_perc_log_agg, r)
crop2real_df <- as.data.frame(crop2real, xy = TRUE)

crop2vali <- crop(vali_ras_perc_log_agg, r)
crop2vali_df <- as.data.frame(crop2vali, xy = TRUE)

sim2 <- ggplot() +
  geom_raster(data = crop2sim_df, aes(x = x, y = y, fill = count)) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(630000, 680000),
                  ylim = c(4222000, 4290000)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Density of simulated points",
       x = "UTM Easting", y = "UTM Northing")

real2 <- ggplot() +
  geom_raster(data = crop2real_df, aes(x = x, y = y, fill = count)) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(630000, 680000),
                  ylim = c(4222000, 4290000)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Density of observed points",
       x = "UTM Easting", y = "UTM Northing")

real_trk2 <- ggplot() +
  geom_path(data = mig_spring, aes(x = x_, y = y_, group = burst),
            color = "black", alpha = 0.1) +
  geom_point(data = mig_spring, aes(x = x_, y = y_),
             color = "black", fill = "black", alpha = 0.1) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(630000, 680000),
                  ylim = c(4222000, 4290000)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Observed tracks",
       x = "UTM Easting", y = "UTM Northing")

vali2 <- ggplot() +
  geom_raster(data = crop2vali_df, aes(x = x, y = y, fill = count)) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(630000, 680000),
                  ylim = c(4222000, 4290000)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Density of validation points",
       x = "UTM Easting", y = "UTM Northing")

vali_trk2 <- ggplot() +
  geom_path(data = vali, aes(x = x, y = y, group = ID),
            color = "black", alpha = 0.1) +
  geom_point(data = vali, aes(x = x, y = y),
             color = "black", fill = "black", alpha = 0.1) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(630000, 680000),
                  ylim = c(4222000, 4290000)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Validation tracks",
       x = "UTM Easting", y = "UTM Northing")

vali_combo <- ggplot() +
  geom_raster(data = crop2sim_df, aes(x = x, y = y, fill = count)) +
  geom_path(data = vali, aes(x = x, y = y, group = ID),
            color = "white", alpha = 0.3) +
  geom_point(data = vali, aes(x = x, y = y),
             color = "white", fill = "white", alpha = 0.6) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(630000, 680000),
                  ylim = c(4222000, 4290000)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Density of simulated points (independent data overlaid)",
       x = "UTM Easting", y = "UTM Northing")

(real_trk2 + real2) / (sim2 + vali_combo)

real_trk2 + real2 + sim2

ggsave("output/real-vs-sim_spring-mig_example1.tiff",
       width = 12, height = 4, dpi = 300, compression = "lzw")

vali_combo

ggsave(plot = vali_combo, "output/sim-with-valid_spring-mig_example1.tiff",
       width = 5, height = 5, dpi = 300, compression = "lzw")

map | real_trk2 | real2 | sim2

ggsave("output/real-vs-sim_spring-mig_example1_with-map.tiff",
       width = 16, height = 4, dpi = 300, compression = "lzw")

map + vali_combo

ggsave("output/sim-with-valid_spring-mig_example1_with-map.tiff",
       width = 10, height = 6, dpi = 300, compression = "lzw")

## Example 2 ####
r <- rast()
ext(r) <- c(520000, 680000,
            4450000, 4580000)
values(r) <- 1

# Who is in this example?

ex2 <- as.data.frame(pts) %>%
  cbind.data.frame(geom(pts)) %>%
  filter(x >= 520000 & x <= 680000 &
           y >= 4450000 & y <= 4580000) %>%
  group_by(id) %>%
  tally() %>%
  arrange(desc(n))

ex2$id

# Map inset
map <- ggplot() +
  geom_sf(data = utah) +
  geom_raster(data = as.data.frame(r, xy = TRUE),
              mapping = aes(x = x, y = y), fill = "tomato") +
  theme_void() +
  theme(legend.position = "none")

crop2sim <- crop(sim_ras_perc_log_agg, r)
crop2sim_df <- as.data.frame(crop2sim, xy = TRUE)

crop2real <- crop(real_ras_perc_log_agg, r)
crop2real_df <- as.data.frame(crop2real, xy = TRUE)

crop2vali <- crop(vali_ras_perc_log_agg, r)
crop2vali_df <- as.data.frame(crop2vali, xy = TRUE)

sim2 <- ggplot() +
  geom_raster(data = crop2sim_df, aes(x = x, y = y, fill = count)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Density of simulated points",
       x = "UTM Easting", y = "UTM Northing")

real2 <- ggplot() +
  geom_raster(data = crop2real_df, aes(x = x, y = y, fill = count)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Density of observed points",
       x = "UTM Easting", y = "UTM Northing")

real_trk2 <- ggplot() +
  geom_path(data = mig_spring, aes(x = x_, y = y_, group = burst),
            color = "black", alpha = 0.1) +
  geom_point(data = mig_spring, aes(x = x_, y = y_),
             color = "black", fill = "black", alpha = 0.1) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(520000, 710000),
                  ylim = c(4450000, 4580000)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Observed tracks",
       x = "UTM Easting", y = "UTM Northing")

vali2 <- ggplot() +
  geom_raster(data = crop2vali_df, aes(x = x, y = y, fill = count)) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(520000, 680000),
                  ylim = c(4450000, 4580000)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Density of validation points",
       x = "UTM Easting", y = "UTM Northing")

vali_trk2 <- ggplot() +
  geom_path(data = vali, aes(x = x, y = y, group = ID),
            color = "black", alpha = 0.1) +
  geom_point(data = vali, aes(x = x, y = y),
             color = "black", fill = "black", alpha = 0.1) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(520000, 680000),
                  ylim = c(4450000, 4580000)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Validation tracks",
       x = "UTM Easting", y = "UTM Northing")

vali_combo <- ggplot() +
  geom_raster(data = crop2sim_df, aes(x = x, y = y, fill = count)) +
  geom_path(data = vali, aes(x = x, y = y, group = ID),
            color = "white", alpha = 0.1) +
  geom_point(data = vali, aes(x = x, y = y),
             color = "white", fill = "white", alpha = 0.1) +
  scale_fill_viridis_c() +
  coord_cartesian(xlim = c(530000, 680000),
                  ylim = c(4450000, 4580000)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Density of simulated points (independent data overlaid)",
       x = "UTM Easting", y = "UTM Northing")

(real_trk2 + real2) / (sim2 + vali_combo)

real_trk2 + real2 + sim2

ggsave("output/real-vs-sim_spring-mig_example2.tiff",
       width = 20, height = 4, dpi = 300, compression = "lzw")

vali_combo

ggsave(plot = vali_combo, "output/sim-with-valid_spring-mig_example2.tiff",
       width = 8, height = 5, dpi = 300, compression = "lzw")

map | real_trk2 | real2 | sim2

ggsave("output/real-vs-sim_spring-mig_example2_with-map.tiff",
       width = 25, height = 5, dpi = 300, compression = "lzw")

map + vali_combo

ggsave("output/sim-with-valid_spring-mig_example2_with-map.tiff",
       width = 12, height = 5, dpi = 300, compression = "lzw")

# Plot tracks ####

## All ####

p1 <- ggplot() +
  geom_path(data = mig_spring, aes(x = x_, y = y_, group = burst),
            color = "black", alpha = 0.1) +
  geom_point(data = mig_spring, aes(x = x_, y = y_),
             color = "black", fill = "black", alpha = 0.1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Observed tracks",
       x = "UTM Easting", y = "UTM Northing")

p2 <- ggplot() +
  geom_path(data = ind_sims, aes(x = x_, y = y_, group = id),
            color = "black", alpha = 0.1) +
  geom_point(data = ind_sims, aes(x = x_, y = y_),
             color = "black", fill = "black", alpha = 0.1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Simulated tracks",
       x = "UTM Easting", y = "UTM Northing")

p3 <- ggplot() +
  geom_path(data = vali, aes(x = x, y = y, group = ID),
            color = "black", alpha = 0.1) +
  geom_point(data = vali, aes(x = x, y = y),
             color = "black", fill = "black", alpha = 0.1) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Independent observed tracks",
       x = "UTM Easting", y = "UTM Northing")

p1 + p2 + p3

ggsave("output/real-sim-vali_spring_tracks.tiff",
       width = 12, height = 5, dpi = 300, compression = "lzw")

## Individual ####

inds <- unique(mig_spring$deploy_ID)
csus <- mig_spring %>%
  select(deploy_ID, csu) %>%
  distinct()

elev <- rast("output/processed_layers/elevation_utm.tiff")

for (i in 1:length(inds)) {

  print(i)

  csu <- csus[csus$deploy_ID == inds[i], ]$csu

  sim <- ind_sims %>%
    filter(id == inds[i])

  real <- mig_spring %>%
    filter(deploy_ID == inds[i])

  x_lims <- c(min(c(sim$x_, real$x_)) - 1000,
              max(c(sim$x_, real$x_)) + 1000)

  y_lims <- c(min(c(sim$y_, real$y_)) - 1000,
              max(c(sim$y_, real$y_)) + 1000)

  r <- rast()
  ext(r) <- c(x_lims, y_lims)
  el <- crop(elev, r)
  el_df <- as.data.frame(el, xy = TRUE)

  ggplot() +
    geom_raster(data = el_df, aes(x = x, y = y, fill = LC20_Elev_220),
                alpha = 0.3) +
    geom_path(data = sim, aes(x = x_, y = y_, group = sim),
              color = "#148484", alpha = 0.01) +
    geom_point(data = sim, aes(x = x_, y = y_),
               color = "#148484", fill = "#148484", alpha = 0.01) +
    geom_path(data = real, aes(x = x_, y = y_),
              color = "#F8AF67") +
    geom_point(data = real, aes(x = x_, y = y_),
               color = "#F8AF67", fill = "#F8AF67") +
    scale_fill_distiller(type = "seq",
                         direction = -1,
                         palette = "Greys") +
    coord_cartesian(xlim = x_lims, ylim = y_lims) +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    labs(x = "UTM Easting", y = "UTM Northing",
         title = paste0(inds[i], " ", csu))

  ggsave(paste0("output/real-vs-sim_spring_", inds[i], "_with-elev_with-ID.tiff"),
                height = 5, width = 5, dpi = 300, compression = "lzw")

}
