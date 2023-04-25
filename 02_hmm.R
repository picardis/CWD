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

# Re-filter again based on sample size. Not only x and y, but also angle
# needs to not be NA.
keepers <- mule %>%
  filter(!is.na(longitude) & !is.na(latitude) & !is.na(angle)) %>%
  group_by(ID) %>%
  tally() %>%
  arrange(n) %>%
  filter(n >= 100) %>%
  pull(ID)

mule <- mule %>%
  filter(ID %in% keepers)

# Save again
saveRDS(mule, "output/mule-deer_HMM-prep_clean.rds")

mule <- readRDS("output/mule-deer_HMM-prep_clean.rds")

# Explore steps and angles ####

summary(mule$step)
hist(mule$step)

# Choose initial parameter values ####

source("FUN_gamma_pars.R")

# Three states

# Short state
test1 <- gamma_pars(mean = 10, sd = 5)
# Intermediate state
test2 <- gamma_pars(mean = 150, sd = 75) # 150, 200
# Long state
test3 <- gamma_pars(mean = 1000, sd = 500)

hist(mule$step, breaks = 100, freq = FALSE)
lines(dgamma(x = c(0:7000), shape = test1[2], rate = test1[1]), col = "red")
lines(dgamma(x = c(0:7000), shape = test2[2], rate = test2[1]), col = "green")
lines(dgamma(x = c(0:7000), shape = test3[2], rate = test3[1]), col = "blue")

# Two states
# Short state
test1 <- gamma_pars(mean = 30, sd = 15)
# Long state
test2 <- gamma_pars(mean = 600, sd = 300)

hist(mule$step, breaks = 100, freq = FALSE)
lines(dgamma(x = c(0:7000), shape = test1[2], rate = test1[1]), col = "blue")
lines(dgamma(x = c(0:7000), shape = test2[2], rate = test2[1]), col = "green")

# 3-state HMM ####

system.time({
  hmm1 <- fitHMM(data = mule,
               nbStates = 3,
               dist = list(step = "gamma", angle = "vm"),
               Par0 = list(step = c(mean_1 = 10,
                                    mean_2 = 150,
                                    mean_3 = 1000,
                                    sd_1 = 5,
                                    sd_2 = 75,
                                    sd_3 = 500,
                                    zeromass_1 = 0.5,
                                    zeromass_2 = 0.01,
                                    zeromass_3 = 0.01),
                           angle = c(concentration_1 = 0.1,
                                     concentration_2 = 0.5,
                                     concentration_3 = 0.99)))
  })

saveRDS(hmm1, "output/HMM1_2023-04-06.rds")

mule$state3 <- viterbi(hmm1)

mule %>%
  ggplot(aes(x = factor(lubridate::month(t_)), fill = factor(state3), group = factor(state3))) +
  geom_bar() +
  labs(group = "State", fill = "State", y = "Frequency", x = "Month")

# 2-state HMM ####

# Things I have tried:
# - Fiddle with initial parameter values
# - Add estAngleMean = list(angle = TRUE)

# Things to try next:
# - Exponential distribution for step lengths

system.time({
  hmm2 <- fitHMM(data = mule,
                 nbStates = 2,
                 dist = list(step = "gamma", angle = "vm"),
                 estAngleMean = list(angle = TRUE),
                 Par0 = list(step = c(mean_1 = 30,
                                      mean_2 = 600,
                                      sd_1 = 15,
                                      sd_2 = 300,
                                      zeromass_1 = 0.5,
                                      zeromass_2 = 0.01),
                             angle = c(mean_1 = pi,
                                       mean_2 = 0,
                                       concentration_1 = 0.1,
                                       concentration_2 = 0.99)))
})

saveRDS(hmm2, "output/HMM2_2023-04-07.rds")

mule$state <- viterbi(hmm2)

mule %>%
  ggplot(aes(x = factor(lubridate::month(t_)), fill = factor(state), group = factor(state))) +
  geom_bar() +
  labs(fill = "State", y = "Frequency", x = "Month")

ggsave("states_by_month.tiff", dpi= 300, width = 6, height = 4)

# WHAT THE HELL IS WRONG?

mule %>%
  ggplot(aes(x = step)) +
  geom_histogram() +
  facet_wrap(~ lubridate::month(t_))

mule %>%
  ggplot(aes(x = factor(lubridate::month(t_)), fill = factor(ID), group = factor(ID))) +
  geom_bar() +
  labs(fill = "Burst", y = "Frequency", x = "Month") +
  theme(legend.position = "none")

mule %>%
  mutate(na = is.na(longitude) & is.na(latitude)) %>%
  ggplot(aes(x = factor(lubridate::month(t_)), fill = factor(na), group = factor(na))) +
  geom_bar() +
  labs(fill = "Missing data", y = "Frequency", x = "Month")

mule_amt <- mule %>%
  filter(!is.na(longitude)) %>%
  nest(trk = -ID) %>%
  mutate(trk = lapply(trk, function(x) {
    x %>%
      amt::make_track(.x = longitude, .y = latitude, .t = t_, crs = 4326)
  }))

class(mule_amt$trk[[1]])

?amt::inspect

amt::inspect(mule_amt$trk[[1]], cluster = FALSE, popup = mule_amt$trk[[1]]$t_)

mule %>%
  group_by(ID, lubridate::month(t_)) %>%
  tally() %>%
  arrange(n)

# Try next:
## 1: Coarsen to 12 or 24 h
## 2: Train on known states
## 3: Restrict to ~March-June (try with spring first, can either predict or fit)
## 4: Isolate bursts where migration transitions into residency or vice-versa
# We're going to try #1 (me) and #4 (Emily) and then if that doesn't work we try
# #2. If all else fails, try #3.
# Send Kezia the exact dataset I'm working with and share some basic code to
# plot and to inspect with amt!
