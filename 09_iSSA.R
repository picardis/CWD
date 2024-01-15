# Load packages ####

library(amt)
library(tidyverse)

# Load data ####

mig_spring_amt <- readRDS("output/mig_spring_100rsteps_with-dyn-covs_pranges-dist_2024-01-08.rds")
mig_fall_amt <- readRDS("output/mig_fall_100rsteps_with-dyn-covs_pranges-dist.rds")

# Calculate movement variables ####

mig_spring_amt <- mig_spring_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function(x) {
    x %>%
      mutate(log_sl_ = log(sl_),
             cos_ta_ = cos(ta_),
             burst_step_id_ = paste0(burst_, "_", step_id_))
  }))

mig_fall_amt <- mig_fall_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function(x) {
    x %>%
      mutate(log_sl_ = log(sl_),
             cos_ta_ = cos(ta_),
             burst_step_id_ = paste0(burst_, "_", step_id_))
  }))

# Scale and center ####

source("FUN_scale-and-center.R")

vars <- c("elev_start", "elev_end", "dist_to_roads_end_log",
          "ndvi_start", "ndvi_end", "dist_to_range_end_log")

spring_issa_dat <- mig_spring_amt %>%
  dplyr::select(deploy_ID, rsteps) %>%
  unnest(cols = rsteps)
fall_issa_dat <- mig_fall_amt %>%
  dplyr::select(deploy_ID, rsteps) %>%
  unnest(cols = rsteps)

means <- spring_issa_dat %>%
  bind_rows(fall_issa_dat) %>%
  mutate(dist_to_roads_end_log = log(dist_to_roads_end + 0.001),
         dist_to_summ_end_log = log(dist_to_summ_end + 0.001),
         dist_to_wint_end_log = log(dist_to_wint_end + 0.001),
         dist_to_range_end_log = case_when(
           is.na(dist_to_summ_end_log) ~ dist_to_wint_end_log,
           is.na(dist_to_wint_end_log) ~ dist_to_summ_end_log
         )) %>%
  summarize(across(.cols = all_of(vars), .fns = ~ mean(.x, na.rm = TRUE)))

saveRDS(means, "output/means_2023-12-07.rds")

sds <- spring_issa_dat %>%
  bind_rows(fall_issa_dat) %>%
  mutate(dist_to_roads_end_log = log(dist_to_roads_end + 0.001),
         dist_to_summ_end_log = log(dist_to_summ_end + 0.001),
         dist_to_wint_end_log = log(dist_to_wint_end + 0.001),
         dist_to_range_end_log = case_when(
           is.na(dist_to_summ_end_log) ~ dist_to_wint_end_log,
           is.na(dist_to_wint_end_log) ~ dist_to_summ_end_log
         )) %>%
  summarize(across(.cols = all_of(vars), .fns = ~ sd(.x, na.rm = TRUE)))

saveRDS(sds, "output/sds_2023-12-07.rds")

mig_spring_amt <- mig_spring_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function(x) {
    x %>%
      mutate(dist_to_roads_end_log = log(dist_to_roads_end + 0.001),
             dist_to_range_end_log = log(dist_to_summ_end + 0.001)) %>%
      mutate(snow_start = factor(snow_start),
             snow_end = factor(snow_end),
             land_cover_end = factor(land_cover_end, levels = c("Forest",
                                                        "Agricultural",
                                                        "Developed",
                                                        "Grassland",
                                                        "Shrubland",
                                                        "Unsuitable")),
             cliffs_end = factor(cliffs_end)) %>%
      scale_data(means, sds)
  }))

mig_fall_amt <- mig_fall_amt %>%
  mutate(rsteps = lapply(rsteps, FUN = function(x) {
    x %>%
      mutate(dist_to_roads_end_log = log(dist_to_roads_end + 0.001),
             dist_to_range_end_log = log(dist_to_wint_end + 0.001)) %>%
      mutate(snow_start = factor(snow_start),
             snow_end = factor(snow_end),
             land_cover_end = factor(land_cover_end, levels = c("Forest",
                                                                "Agricultural",
                                                                "Developed",
                                                                "Grassland",
                                                                "Shrubland",
                                                                "Unsuitable")),
             cliffs_end = factor(cliffs_end)) %>%
      scale_data(means, sds)
  }))

# Fit model ####

## Spring ####

x <- mig_spring_amt$rsteps[[1]]

# Many individuals never have steps that cross roads (not even available)
x %>%
  mutate(road_cross = road_poly_start != road_poly_end) %>%
  pull(road_cross) %>%
  table()

# How many individuals have at least one step that cross a road?
mig_spring_amt %>%
  dplyr::select(deploy_ID, rsteps) %>%
  unnest(cols = rsteps) %>%
  mutate(road_cross = road_poly_start != road_poly_end) %>%
  group_by(deploy_ID, road_cross) %>%
  tally() %>%
  group_by(deploy_ID) %>%
  tally() %>%
  ungroup() %>%
  pull(n) %>%
  table()
# 53 out of 152.

# How many individuals have at least one step on a cliff?
mig_spring_amt %>%
  dplyr::select(deploy_ID, rsteps) %>%
  unnest(cols = rsteps) %>%
  group_by(deploy_ID, cliffs_end) %>%
  tally() %>%
  group_by(deploy_ID) %>%
  tally() %>%
  ungroup() %>%
  pull(n) %>%
  table()
# 136 out of 152

# Both of those are problematic and cause convergence problems.
# Land cover is problematic too, some classes are not represented enough.
# Update: I now simplified land cover.

# Because elev_start is always the same in a given stratum, elev_end_sc
# and I(elev_end_sc - elev_start_sc) are perfectly correlated. Can't have both.

# How many individuals that never encounter snow?
mig_fall_amt %>%
  dplyr::select(deploy_ID, rsteps) %>%
  unnest(cols = rsteps) %>%
  group_by(deploy_ID, snow_end) %>%
  tally() %>%
  group_by(deploy_ID) %>%
  tally() %>%
  ungroup() %>%
  pull(n) %>%
  table()
# 59 out of 152.

issa_spring <- mig_spring_amt %>%
  mutate(issa = lapply(rsteps, FUN = function (x) {
    res <- NA
    try(res <- fit_issf(data = x,
                        formula = case_ ~
                          land_cover_end +
                          I(elev_end_sc - elev_start_sc) +
                          I((elev_end_sc - elev_start_sc)^2) +
                          elev_start_sc : I(elev_end_sc - elev_start_sc) +
                          elev_start_sc : I((elev_end_sc - elev_start_sc)^2) +
                          dist_to_roads_end_log_sc +
                          dist_to_range_end_log_sc +
                          cos_ta_:dist_to_range_end_log_sc +
                          #I(road_poly_start != road_poly_end) +
                          #cliffs_end +
                          I(ndvi_end_sc - ndvi_start_sc) +
                          ndvi_start_sc : I(ndvi_end_sc - ndvi_start_sc) +
                          #snow_end +
                          sl_ +
                          log_sl_ +
                          cos_ta_ +
                          cos_ta_:log_sl_ +
                          strata(burst_step_id_),
                        model = TRUE))
  })) %>%
  select(deploy_ID, issa)

sum(lapply(issa_spring$issa, class) == "try-error")
# 38 out of 146 do not converge because they're the ones that didn't have a
# capture unit and therefore didn't have any values of distance to range
issa_spring <- issa_spring[!lapply(issa_spring$issa, class) == "try-error", ]
# Filter them out for now

saveRDS(issa_spring, "output/iSSA_spring_2024-01-08.rds")

## Fall ####

issa_fall <- mig_fall_amt %>%
  mutate(issa = lapply(rsteps, FUN = function (x) {
    res <- NA
    try(res <- fit_issf(data = x,
                        formula = case_ ~
                          land_cover_end +
                          I(elev_end_sc - elev_start_sc) +
                          I((elev_end_sc - elev_start_sc)^2) +
                          elev_start_sc : I(elev_end_sc - elev_start_sc) +
                          elev_start_sc : I((elev_end_sc - elev_start_sc)^2) +
                          dist_to_roads_end_log_sc +
                          dist_to_range_end_log_sc +
                          #I(road_poly_start != road_poly_end) +
                          #cliffs_end +
                          I(ndvi_end_sc - ndvi_start_sc) +
                          ndvi_start_sc : I(ndvi_end_sc - ndvi_start_sc) +
                          #snow_end +
                          sl_ +
                          log_sl_ +
                          cos_ta_ +
                          strata(burst_step_id_),
                        model = TRUE))
  })) %>%
  select(deploy_ID, issa)

sum(lapply(issa_fall$issa, class) == "try-error")
# 38 out of 146 do not converge because they're the ones that didn't have a
# capture unit and therefore didn't have any values of distance to range
issa_fall <- issa_fall[!lapply(issa_fall$issa, class) == "try-error", ]
# Filter them out for now

saveRDS(issa_fall, "output/iSSA_fall_2023-12-07.rds")
# 2023-07-26 is the last one I ran. Compare the new one to that (only difference
# is the addition of dist_to_range_end_log_sc)

# Explore results ####

# Look at elevation gain parameters

b <- mean(sapply(issa_spring$issa, function (x) {
  x$model$coefficients[["I(elev_end_sc - elev_start_sc)"]]
}))

a <- mean(sapply(issa_spring$issa, function (x) {
  x$model$coefficients[["I((elev_end_sc - elev_start_sc)^2)"]]
}))

-1*b/(2*a)

# Format model output
issa_spring <- issa_spring %>%
  mutate(coef = map(issa, function(x) broom::tidy(x$model)))

# Unnest coefficients
coefs <- issa_spring %>%
  select(-issa) %>%
  unnest(cols = coef)

# Plot coefficients
ggplot(coefs, aes(x = term, y = estimate)) +
  geom_boxplot() +
  # coord_cartesian(ylim = c(-10, 10)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Make prediction scenarios
prediction_data <- function(elev_end = means$elev_end,
                            elev_start = means$elev_start,
                            dist_to_roads_end_log = means$dist_to_roads_end_log,
                            sl_ = 100,
                            log_sl_ = log(100),
                            cos_ta_ = 1){
  df <- expand.grid(elev_end = elev_end,
                    elev_start = elev_start,
                    dist_to_roads_end_log = dist_to_roads_end_log,
                    sl_ = sl_,
                    log_sl_ = log_sl_,
                    cos_ta_ = cos_ta_)
  return(df)
}

# Function to scale values so I can define scenarios on interpretable values
source("FUN_scale-and-center.R")

# Scenario 1: predict selection for a 500m elevation gain starting at 2500m
issa_spring <- issa_spring %>%
  mutate(elev_pred = map(issa, function(m) {
    elev_st = 2500
    # x1
    x1 <- prediction_data(elev_start = elev_st,
                     elev_end = seq(elev_st - 500,
                                    elev_st + 500, length.out = 100)) %>%
      scale_data(means, sds) %>%
      mutate(elev_gain = elev_end - elev_start)
    # x2
    x2 <- prediction_data(elev_start = elev_st,
                     elev_end = elev_st) %>%
      scale_data(means, sds) %>%
      mutate(elev_gain = elev_end - elev_start)
    # log-RSS
    lr <- log_rss(m, x1 = x1, x2 = x2)
    # Return just the data.frame
    return(lr$df)
  }))

elev_data <- issa_spring %>%
  select(deploy_ID, elev_pred) %>%
  unnest(cols = elev_pred)

elev_mean <- elev_data %>%
  group_by(elev_gain_x1) %>%
  summarize(log_rss = median(log_rss))

ggplot() +
  geom_line(data = elev_data, aes(x = elev_gain_x1, y = log_rss,
                                  group = deploy_ID),
            alpha = 0.1) +
  geom_line(data = elev_mean, aes(x = elev_gain_x1, y = log_rss),
            color = "red") +
  coord_cartesian(ylim = c(-100, 100)) +
  theme_bw()
