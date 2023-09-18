# Prepare covariate layers for iSSA

# This script is analogous to 07_prep-iSSA-covariates.R, except that here I am
# using the data at 12h resolution. We decided that we care more about having
# as many individuals as possible than about having the highest possible time
# resolution.

# Load packages ####

library(tidyverse)
library(terra)
library(patchwork)
library(amt)
library(lubridate)
library(sf)

# Set seed ####

set.seed(22)

# Load HMM output ####

mule12 <- readRDS("output/mule-deer_regularized-12h_CTMM_with-viterbi.rds") %>%
  as.data.frame()

length(unique(mule12$deploy_ID))

# Subset migration only  and format for amt ####

mule_mig <- mule12 %>%
  filter(vit_3states == 3)

hist(month(mule_mig$t_))

# Spring
mig_spring_amt <- mule_mig %>%
  filter(month(t_) %in% 3:8) %>%
  filter(!is.na(x) & !is.na(y)) %>%
  nest(cols = -deploy_ID) %>%
  mutate(trk = lapply(cols, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      make_track(.x = x, .y = y, .t = t_) %>%
      track_resample(rate = hours(12))
  })) %>%
  mutate(steps = lapply(trk, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      steps_by_burst(keep_cols = "end")
  }))

# Fall
mig_fall_amt <- mule_mig %>%
  filter(month(t_) %in% c(1:2, 9:12)) %>%
  filter(!is.na(x) & !is.na(y)) %>%
  nest(cols = -deploy_ID) %>%
  mutate(trk = lapply(cols, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      make_track(.x = x, .y = y, .t = t_) %>%
      track_resample(rate = hours(12))
  })) %>%
  mutate(steps = lapply(trk, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      steps_by_burst(keep_cols = "end")
  }))

# Remove IDs with few steps and steps with no turn angle ####

mig_spring_amt <- mig_spring_amt %>%
  mutate(steps = lapply(steps, function(x) {
    x %>%
      filter(!is.na(ta_))
  }),
  n_steps = map_dbl(steps, nrow)) %>%
  filter(n_steps >= 20)

mig_fall_amt <- mig_fall_amt %>%
  mutate(steps = lapply(steps, function(x) {
    x %>%
      filter(!is.na(ta_))
  }),
  n_steps = map_dbl(steps, nrow)) %>%
  filter(n_steps >= 20)

# Already here I ended up with 154 individuals in the spring and 177 in the fall.
# That's less than I ended up with when using the 2h data (152 and 203). Stick to
# 2 h then. Unless I include all months of the year (Mar-Aug and Sep-Feb) and
# then I end up with 211 and 259. Worth it? I risk contaminating the information
# with spurious migrations. I think I'm better off using few data points that I
# am confident about.
