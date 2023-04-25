# Explore (plot and inspect) deer movement tracks

# Load packages ####

library(tidyverse)
library(amt)

# Load data ####

# Make sure you adjust this path to match the location of the file on your computer
mule <- readRDS("output/mule-deer_HMM-prep_clean.rds")

# Plot one track at a time ####

# Which track do you want to plot?
who <- 1 # replace this with the numbers from 1 to length(unique(mule$ID))

mule %>%
  # Grab data for the current individual
  filter(ID == unique(mule$ID)[who]) %>%
  ggplot(aes(x = x, y = y)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  labs(x = "UTM X", y = "UTM Y")

# Dynamically inspect one track at a time with amt ####

mule_amt <- mule %>%
  filter(!is.na(longitude)) %>%
  nest(trk = -ID) %>%
  mutate(trk = lapply(trk, function(x) {
    x %>%
      amt::make_track(.x = longitude, .y = latitude, .t = t_, crs = 4326)
  }))

amt::inspect(mule_amt$trk[[1]], cluster = FALSE, popup = mule_amt$trk[[1]]$t_)
