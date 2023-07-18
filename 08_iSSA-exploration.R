# Load packages ####

library(amt)
library(tidyverse)

# Load data ####

mig_spring_amt <- readRDS("output/mig_spring_100rsteps_with-covs.rds")
mig_fall_amt <- readRDS("output/mig_fall_100rsteps_with-covs.rds")

# Explore ####

mig_spring_amt <- mig_spring_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function (x) {
  x %>%
    mutate(elev_gain = elev_end - elev_start,
           road_cross = road_poly_start != road_poly_end)
})) %>%
  select(deploy_ID, rsteps)

mig_fall_amt <- mig_fall_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function (x) {
    x %>%
      mutate(elev_gain = elev_end - elev_start,
             road_cross = road_poly_start != road_poly_end)
  })) %>%
  select(deploy_ID, rsteps)

mig_spring_explo <- unnest(mig_spring_amt, cols = rsteps)
mig_fall_explo <- unnest(mig_fall_amt, cols = rsteps)

mig_explo <- mig_spring_explo %>%
  mutate(season = "Spring") %>%
  bind_rows(mig_fall_explo) %>%
  mutate(season = case_when(
    season == "Spring" ~ "Spring",
    TRUE ~ "Fall"
  ))

## Land cover ####

n_by_case <- mig_explo %>%
  group_by(land_cover_end, season, case_) %>%
  tally() %>%
  rename(n_by_case = n)

n_tot <- mig_explo %>%
  group_by(case_, season) %>%
  tally()

land_cover_end_explo <- left_join(n_by_case, n_tot) %>%
  mutate(freq = n_by_case/n) %>%
  mutate(step_type = case_when(
    case_ ~ "Used",
    TRUE ~ "Available"
  )) %>%
  filter(!is.na(land_cover_end))

ggplot(land_cover_end_explo, aes(x = season, y = freq, fill = step_type)) +
  geom_bar(stat = "identity", position = position_dodge(0.8)) +
  facet_wrap(~ land_cover_end, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(fill = " ", x = "Land cover", y = "Rel. frequency")

ggsave("output/used-vs-available_land-cover.tiff", width = 10, height = 6,
       dpi = 300, compression = "lzw")

## Elevation ####

mig_explo %>%
  mutate(step_type = case_when(
    case_ ~ "Used",
    TRUE ~ "Available"
  )) %>%
  ggplot(aes(x = step_type, y = elev_end, fill = step_type)) +
  geom_boxplot() +
  facet_wrap(~ season) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(fill = " ", x = " ", y = "Elevation")

ggsave("output/used-vs-available_elevation.tiff", width = 5, height = 4,
       dpi = 300, compression = "lzw")

mig_explo %>%
  mutate(step_type = case_when(
    case_ ~ "Used",
    TRUE ~ "Available"
  )) %>%
  ggplot(aes(x = elev_end, color = step_type)) +
  geom_density() +
  facet_wrap(~ season) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(fill = " ", x = " ", y = "Elevation")

ggsave("output/used-vs-available_density_elevation.tiff", width = 5, height = 4,
       dpi = 300, compression = "lzw")

## Elevation gain ####

mig_explo %>%
  mutate(step_type = case_when(
    case_ ~ "Used",
    TRUE ~ "Available"
  )) %>%
  ggplot(aes(x = step_type, y = elev_gain, fill = step_type)) +
  geom_violin() +
  facet_wrap(~ season) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(fill = " ", x = " ", y = "Elevation gain")

ggsave("output/used-vs-available_elevation-gain_v2.tiff", width = 5, height = 4,
       dpi = 300, compression = "lzw")

mig_explo %>%
  mutate(step_type = case_when(
    case_ ~ "Used",
    TRUE ~ "Available"
  )) %>%
  ggplot(aes(x = elev_gain, color = step_type)) +
  geom_density() +
  facet_wrap(~ season) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(fill = " ", x = " ", y = "Elevation gain")

ggsave("output/used-vs-available_density_elevation-gain.tiff", width = 5, height = 4,
       dpi = 300, compression = "lzw")

## Distance to roads ####

mig_explo %>%
  mutate(step_type = case_when(
    case_ ~ "Used",
    TRUE ~ "Available"
  )) %>%
  ggplot(aes(x = step_type, y = dist_to_roads_end, fill = step_type)) +
  geom_boxplot() +
  facet_wrap(~ season) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(fill = " ", x = " ", y = "Distance to roads")

ggsave("output/used-vs-available_distance-to-roads.tiff", width = 5, height = 4,
       dpi = 300, compression = "lzw")

mig_explo %>%
  mutate(step_type = case_when(
    case_ ~ "Used",
    TRUE ~ "Available"
  )) %>%
  ggplot(aes(x = dist_to_roads_end, color = step_type)) +
  geom_density() +
  facet_wrap(~ season) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(fill = " ", x = " ", y = "Distance to roads")

ggsave("output/used-vs-available_density_distance-to-roads.tiff", width = 5, height = 4,
       dpi = 300, compression = "lzw")

## cliffs_end ####

n_by_case <- mig_explo %>%
  group_by(cliffs_end, season, case_) %>%
  tally() %>%
  rename(n_by_case = n)

n_tot <- mig_explo %>%
  group_by(case_, season) %>%
  tally()

cliffs_end_explo <- left_join(n_by_case, n_tot) %>%
  mutate(freq = n_by_case/n) %>%
  mutate(step_type = case_when(
    case_ ~ "Used",
    TRUE ~ "Available"
  )) %>%
  filter(!is.na(cliffs_end)) %>%
  mutate(cliffs_end = factor(cliffs_end))

ggplot(cliffs_end_explo, aes(x = cliffs_end, y = freq, fill = cliffs_end)) +
  geom_bar(stat = "identity") +
  facet_grid(step_type ~ season) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(fill = "Cliff", x = "cliffs_end", y = "Rel. frequency")

ggsave("output/used-vs-available_cliffs_end.tiff", width = 5, height = 4,
       dpi = 300, compression = "lzw")

## Road crossing ####

n_by_case <- mig_explo %>%
  group_by(road_cross, season, case_) %>%
  tally() %>%
  rename(n_by_case = n)

n_tot <- mig_explo %>%
  group_by(case_, season) %>%
  tally()

road_cross_explo <- left_join(n_by_case, n_tot) %>%
  mutate(freq = n_by_case/n) %>%
  mutate(step_type = case_when(
    case_ ~ "Used",
    TRUE ~ "Available"
  )) %>%
  filter(!is.na(road_cross)) %>%
  mutate(road_cross = factor(road_cross))

ggplot(road_cross_explo, aes(x = road_cross, y = freq, fill = road_cross)) +
  geom_bar(stat = "identity") +
  facet_grid(step_type ~ season) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(fill = "Road crossing", x = "Road crossing", y = "Rel. frequency")

ggsave("output/used-vs-available_road-crossings.tiff", width = 5, height = 4,
       dpi = 300, compression = "lzw")
