# Load packages ####

library(tidyverse)
library(terra)
library(amt)
library(lubridate)
library(sf)

# Set seed ####

set.seed(22)

# Load data ####

mig_spring_amt <- readRDS("output/mig-trk_spring_2h-res.rds")
mig_fall_amt <- readRDS("output/mig-trk_fall_2h-res.rds")

# Fit steps and angle distributions ####

# Steps

# Spring
sl_spring <- mig_spring_amt %>%
  dplyr::select(deploy_ID_year, steps) %>%
  unnest(cols = steps) %>%
  filter(!is.na(sl_) & sl_ > 0) %>%
  pull(sl_)

plot(density(sl_spring))

sl_distr_spring <- fit_distr(sl_spring, dist_name = "gamma")

# Fall
sl_fall <- mig_fall_amt %>%
  dplyr::select(deploy_ID_year, steps) %>%
  unnest(cols = steps) %>%
  filter(!is.na(sl_) & sl_ > 0) %>%
  pull(sl_)

plot(density(sl_fall))

sl_distr_fall <- fit_distr(sl_fall, dist_name = "gamma")

# Angles

# Spring
ta_spring <- mig_spring_amt %>%
  dplyr::select(deploy_ID_year, steps) %>%
  unnest(cols = steps) %>%
  filter(!is.na(ta_)) %>%
  pull(ta_)

plot(density(ta_spring))

ta_distr_spring <- fit_distr(ta_spring, dist_name = "vonmises")

# Fall
ta_fall <- mig_fall_amt %>%
  dplyr::select(deploy_ID_year, steps) %>%
  unnest(cols = steps) %>%
  filter(!is.na(ta_)) %>%
  pull(ta_)

plot(density(ta_fall))

ta_distr_fall <- fit_distr(ta_fall, dist_name = "vonmises")

# Sample random steps ####

mig_spring_amt <- mig_spring_amt %>%
  mutate(rsteps = lapply(steps, FUN = function (x) {
    x %>% random_steps(n_control = 100,
                       sl_distr = sl_distr_spring,
                       ta_distr = ta_distr_spring)
  }))

mig_fall_amt <- mig_fall_amt %>%
  mutate(rsteps = lapply(steps, FUN = function (x) {
    x %>% random_steps(n_control = 100,
                       sl_distr = sl_distr_fall,
                       ta_distr = ta_distr_fall)
  }))

saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_2024-06-19.rds")
saveRDS(mig_fall_amt, "output/mig_fall_100rsteps_2024-06-19.rds")

# Create extent ####

# Compute minimum and maximum coordinate values

minxs <- min(c(
  unlist(lapply(mig_spring_amt$rsteps, FUN = function(x) {
    xs <- c(x$x1_, x$x2_)
    minx <- min(xs)
    return(minx)
  })),
  unlist(lapply(mig_fall_amt$rsteps, FUN = function(x) {
    xs <- c(x$x1_, x$x2_)
    minx <- min(xs)
    return(minx)
  }))
))

maxxs <- max(c(
  unlist(lapply(mig_spring_amt$rsteps, FUN = function(x) {
    xs <- c(x$x1_, x$x2_)
    maxx <- max(xs)
    return(maxx)
  })),
  unlist(lapply(mig_fall_amt$rsteps, FUN = function(x) {
    xs <- c(x$x1_, x$x2_)
    maxx <- max(xs)
    return(maxx)
  }))
))

minys <- min(c(
  unlist(lapply(mig_spring_amt$rsteps, FUN = function(x) {
    ys <- c(x$y1_, x$y2_)
    miny <- min(ys)
    return(miny)
  })),
  unlist(lapply(mig_fall_amt$rsteps, FUN = function(x) {
    ys <- c(x$y1_, x$y2_)
    miny <- min(ys)
    return(miny)
  }))
))

maxys <- max(c(
  unlist(lapply(mig_spring_amt$rsteps, FUN = function(x) {
    ys <- c(x$y1_, x$y2_)
    maxy <- max(ys)
    return(maxy)
  })),
  unlist(lapply(mig_fall_amt$rsteps, FUN = function(x) {
    ys <- c(x$y1_, x$y2_)
    maxy <- max(ys)
    return(maxy)
  }))
))

ext_mig <- ext(c(minxs, maxxs, minys, maxys))
sa <- as.polygons(ext_mig, crs = "epsg:32612")

writeVector(sa, "output/study_area_polygon.shp", overwrite = TRUE)

# Take a look ####

# mig_spring_amt %>%
#   unnest(rsteps) %>%
#   ggplot() +
#   geom_path(aes(x = x2_, y = y2_, group = deploy_ID, color = deploy_ID)) +
#   facet_wrap(~ case_) +
#   theme(legend.position = "none")

# Looks good but takes forever to plot
