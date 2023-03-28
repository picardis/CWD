# Load packages ####

library(tidyverse)
library(momentuHMM)

# Load data ####

mule <- readRDS("output/mule-deer_regularized-2h_with-NAs.rds")

# Format for momentuHMM ####

mule <- mule %>%
  rename(ID = burst) %>%
  arrange(ID, t_) %>%
  as.data.frame()

mule <- prepData(mule,
                 coordNames = c("x_", "y_"),
                 type = "UTM")

# Save because it took a long time to run
saveRDS(mule, "output/mule-deer_HMM-prep.rds")

mule <- readRDS("output/mule-deer_HMM-prep.rds")

# Explore steps and angles ####

summary(mule$step)
hist(mule$step)

# Choose initial parameter values ####

source("FUN_gamma_pars.R")

# Short state
test1 <- gamma_pars(mean = 10, sd = 10)
# Intermediate state
test2 <- gamma_pars(mean = 150, sd = 150) # 150, 200
# Long state
test3 <- gamma_pars(mean = 500, sd = 500)

hist(mule$step, breaks = 100, freq = FALSE)
lines(dgamma(x = c(0:7000), shape = test1[2], rate = test1[1]), col = "red")
lines(dgamma(x = c(0:7000), shape = test2[2], rate = test2[1]), col = "green")
lines(dgamma(x = c(0:7000), shape = test3[2], rate = test3[1]), col = "blue")

# 3-state HMM ####

system.time({
  hmm1 <- fitHMM(data = mule,
               nbStates = 3,
               dist = list(step = "gamma", angle = "vm"),
               Par0 = list(step = c(mean_1 = 10,
                                    mean_2 = 150,
                                    mean_3 = 500,
                                    sd_1 = 10,
                                    sd_2 = 150,
                                    sd_3 = 500,
                                    zeromass_1 = 0.5,
                                    zeromass_2 = 0.01,
                                    zeromass_3 = 0.01),
                           angle = c(concentration_1 = 0.1,
                                     concentration_2 = 0.5,
                                     concentration_3 = 0.99)))
  })

saveRDS(hmm1, "output/HMM1_2023-03-27.rds")
