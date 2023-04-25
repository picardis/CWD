# Load packages ####

library(tidyverse)
library(momentuHMM)

# Load data ####

mule <- readRDS("output/mule-deer_regularized-12h_CTMM_NA-gaps.rds")

# Format for momentuHMM ####

mule <- mule %>%
  arrange(ID, t_) %>%
  as.data.frame()

mule <- prepData(mule,
                 coordNames = c("x_", "y_"),
                 type = "UTM")

# Save
saveRDS(mule, "output/mule-deer_HMM-prep_12h_NA-gaps.rds")

mule <- readRDS("output/mule-deer_HMM-prep_12h_NA-gaps.rds")

# Explore steps and angles ####

summary(mule$step)
hist(mule$step, breaks = 1000, xlim = c(1, 10000))

# Choose initial parameter values ####

source("FUN_gamma_pars.R")

# Two states
# Short state
test1 <- gamma_pars(mean = 100, sd = 50)
# Long state
test2 <- gamma_pars(mean = 1000, sd = 500)

hist(mule$step, breaks = 1000, xlim = c(1, 10000), freq = FALSE)
lines(dgamma(x = c(0:10000), shape = test1[2], rate = test1[1]), col = "blue")
lines(dgamma(x = c(0:10000), shape = test2[2], rate = test2[1]), col = "green")

# 2-state HMM ####

system.time({
  hmm2 <- fitHMM(data = mule,
                 nbStates = 2,
                 dist = list(step = "gamma", angle = "vm"),
                 estAngleMean = list(angle = TRUE),
                 Par0 = list(step = c(mean_1 = 100,
                                      mean_2 = 1000,
                                      sd_1 = 50,
                                      sd_2 = 500,
                                      zeromass_1 = 0.5,
                                      zeromass_2 = 0.01),
                             angle = c(mean_1 = pi,
                                       mean_2 = 0,
                                       concentration_1 = 0.1,
                                       concentration_2 = 0.99)))
})

saveRDS(hmm2, "output/HMM2_12h_2023-04-25_NAs.rds")
beepr::beep(2)

