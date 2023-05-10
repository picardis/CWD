# Make a map of migration corridors based on the state assignments from the HMM

# Load packages ####

library(momentuHMM)
library(raster)
library(tidyverse)
library(patchwork)
library(sf)

# Load HMM ####

hmm2 <- readRDS("output/HMM2_12h_2023-04-25.rds")
hmm3 <- readRDS("output/HMM3_12h_2023-05-09.rds")

# Get Viterbi assignments ####

vit_2states <- viterbi(hmm2)
vit_3states <- viterbi(hmm3)

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
              mapping = aes(x = x, y = y, fill = factor(layer))) +
  scale_fill_discrete(na.value = NA) +
  theme_void()

p_res | p_res_yn | p_mig | p_mor | plot_layout(guides = "collect")

ggplot() +
  geom_sf(data = utah) +
  geom_raster(data = mor_df,
              mapping = aes(x = x, y = y, fill = layer)) +
  scale_fill_viridis_c(na.value = NA) +
  theme_void()

# Lines ####

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
  filter(n == 1)

mig_no_nas <- mig_no_nas %>%
  filter(!burst %in% toss)

mig_no_nas %>%
  group_by(burst) %>%
  tally() %>%
  arrange(n)

# mig_pts <- split(mig_no_nas[, c("x", "y")],
#                  f = mig_no_nas$burst)
# mig_pts <- lapply(mig_pts, as.matrix)
# mig_lines <- st_multilinestring(mig_pts, dim = "XY")

mig_lines <- mig_no_nas %>%
  st_as_sf(coords = c("x", "y"), crs = 32612) %>%
  group_by(burst) %>%
  summarize(do_union = FALSE) %>%
  st_cast("LINESTRING")

ggplot() +
  geom_sf(data = mig_lines, aes(color = burst)) +
  theme(legend.position = "none")

mig_lines_rast <- rasterize(mig_lines, template, fun = "count")

# EXPERIMENTING WITH CODE BELOW, STOPPED HERE

rast1 <- template
values(rast1) <- 1:ncell(rast1)

sp <- raster::extract(rast1, mig_lines)

count <- as.data.frame(table(unlist(sp)))

newraster <- rast1
newraster <- subs(rast1, count, by = "test")

ggplot() +
  geom_sf(data = utah) +
  geom_sf(data = mig_lines) +
  scale_color_viridis_c(na.value = NA) +
  theme_void()

# Brownian bridges ####

library(move)

mig_move <- move::move(x = mig$x, y = mig$y, animal = mig$ID,
                 time = mig$t_, proj = CRS("+init=epsg:32612"))

bb <- brownian.bridge.dyn(mig_move,
                          raster = template,
                          location.error = 15)
