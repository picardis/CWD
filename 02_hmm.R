# Load packages ####

library(tidyverse)
library(momentuHMM)

# Load data ####

mule <- readRDS("output/mule-deer_regularized-2h_with-NAs.rds")

# Format for momentuHMM ####

mule <- mule %>%
  rename(ID = burst) %>%
  arrange(ID, timestamp) %>%
  as.data.frame()

mule <- prepData(mule,
                 coordNames = c("utm_x", "utm_y"),
                 type = "UTM")

# Save because it took a long time to run
saveRDS(mule, "output/mule-deer_HMM-prep.rds")

# Explore steps and angles ####

summary(mule$step)
hist(mule$step)
