# Load packages ####

library(tidyverse)
library(momentuHMM)
library(lubridate)
library(patchwork)

# Load data ####

mule <- readRDS("output/mule-deer_regularized-12h_CTMM_no-gaps.rds")

# Format for momentuHMM ####

mule <- mule %>%
  arrange(ID, t_) %>%
  as.data.frame()

mule <- prepData(mule,
                 coordNames = c("x_", "y_"),
                 type = "UTM")

# Save
saveRDS(mule, "output/mule-deer_HMM-prep_12h_no-gaps.rds")

mule <- readRDS("output/mule-deer_HMM-prep_12h_no-gaps.rds")

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

# Three states
# Short state
test1 <- gamma_pars(mean = 50, sd = 25)
# Medium state
test2 <- gamma_pars(mean = 300, sd = 150)
# Long state
test3 <- gamma_pars(mean = 1000, sd = 500)

hist(mule$step, breaks = 1000, xlim = c(1, 10000), freq = FALSE)
lines(dgamma(x = c(0:10000), shape = test1[2], rate = test1[1]), col = "blue")
lines(dgamma(x = c(0:10000), shape = test2[2], rate = test2[1]), col = "green")
lines(dgamma(x = c(0:10000), shape = test3[2], rate = test3[1]), col = "red")

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

saveRDS(hmm2, "output/HMM2_12h_2023-04-25.rds")
beepr::beep(2)

# 3-state HMM ####

system.time({
  hmm3 <- fitHMM(data = mule,
                 nbStates = 3,
                 dist = list(step = "gamma", angle = "vm"),
                 estAngleMean = list(angle = TRUE),
                 Par0 = list(step = c(mean_1 = 50,
                                      mean_2 = 300,
                                      mean_3 = 1000,
                                      sd_1 = 25,
                                      sd_2 = 150,
                                      sd_3 = 500,
                                      zeromass_1 = 0.5,
                                      zeromass_2 = 0.1,
                                      zeromass_3 = 0.001),
                             angle = c(mean_1 = pi,
                                       mean_2 = pi,
                                       mean_3 = 0,
                                       concentration_1 = 0.1,
                                       concentration_2 = 0.5,
                                       concentration_3 = 0.99)))
})

saveRDS(hmm3, "output/HMM3_12h_2023-05-09.rds")
beepr::beep(2)

# Get Viterbi assignments ####

mule$vit_2states <- viterbi(hmm2)
mule$vit_3states <- viterbi(hmm3)

# Check results ####

# Check distribution of states throughout the year
# (I expect the long movement state to be mostly in migration times)

data_by_month <- mule %>%
  mutate(month = lubridate::month(t_)) %>%
  group_by(month) %>%
  tally()

s3 <- mule %>%
  mutate(month = lubridate::month(t_)) %>%
  left_join(data_by_month, by = "month") %>%
  group_by(month, vit_3states, n) %>%
  tally() %>%
  mutate(perc_state_by_month = nn/n) %>%
  ggplot(aes(x = factor(month),
             y = perc_state_by_month,
             fill = factor(vit_3states))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Month", y = "Frequency", fill = "State", title = "3-state HMM") +
  theme_bw()

s2 <- mule %>%
  mutate(month = lubridate::month(t_)) %>%
  left_join(data_by_month, by = "month") %>%
  group_by(month, vit_2states, n) %>%
  tally() %>%
  mutate(perc_state_by_month = nn/n) %>%
  ggplot(aes(x = factor(month),
             y = perc_state_by_month,
             fill = factor(vit_2states))) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Month", y = "Frequency", fill = "State", title = "2-state HMM") +
  theme_bw()

s2 | s3

ggsave("output/states-by-month_2-vs-3-states_2023-05-10.tiff",
       width = 12, height = 5, dpi = 300, compression = "lzw")

# Example tracks ####

# MD14F0001_33768_15

i3_1 <- mule %>%
  filter(ID == "MD14F0001_33768_15") %>%
  mutate(vit_3states = factor(vit_3states)) %>%
  ggplot(aes(x = x, y = y, color = vit_3states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_1 <- mule %>%
  filter(ID == "MD14F0001_33768_15") %>%
  mutate(vit_2states = factor(vit_2states)) %>%
  ggplot(aes(x = x, y = y, color = vit_2states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_1 | i3_1

# MD14F0001_33768_4

i3_2 <- mule %>%
  filter(ID == "MD14F0001_33768_4") %>%
  mutate(vit_3states = factor(vit_3states)) %>%
  ggplot(aes(x = x, y = y, color = vit_3states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_2 <- mule %>%
  filter(ID == "MD14F0001_33768_4") %>%
  mutate(vit_2states = factor(vit_2states)) %>%
  ggplot(aes(x = x, y = y, color = vit_2states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_2 | i3_2

# MD14F0002_33763_25

i3_3 <- mule %>%
  filter(ID == "MD14F0002_33763_25") %>%
  mutate(vit_3states = factor(vit_3states)) %>%
  ggplot(aes(x = x, y = y, color = vit_3states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_3 <- mule %>%
  filter(ID == "MD14F0002_33763_25") %>%
  mutate(vit_2states = factor(vit_2states)) %>%
  ggplot(aes(x = x, y = y, color = vit_2states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_3 | i3_3

# MD14F0002_33763_3

i3_4 <- mule %>%
  filter(ID == "MD14F0002_33763_3") %>%
  mutate(vit_3states = factor(vit_3states)) %>%
  ggplot(aes(x = x, y = y, color = vit_3states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_4 <- mule %>%
  filter(ID == "MD14F0002_33763_3") %>%
  mutate(vit_2states = factor(vit_2states)) %>%
  ggplot(aes(x = x, y = y, color = vit_2states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_4 | i3_4

# MD14F0008_33876_10

i3_5 <- mule %>%
  filter(ID == "MD14F0008_33876_10") %>%
  mutate(vit_3states = factor(vit_3states)) %>%
  ggplot(aes(x = x, y = y, color = vit_3states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_5 <- mule %>%
  filter(ID == "MD14F0008_33876_10") %>%
  mutate(vit_2states = factor(vit_2states)) %>%
  ggplot(aes(x = x, y = y, color = vit_2states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_5 | i3_5

# MD14F0008_33876_19

i3_6 <- mule %>%
  filter(ID == "MD14F0008_33876_19") %>%
  mutate(vit_3states = factor(vit_3states)) %>%
  ggplot(aes(x = x, y = y, color = vit_3states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_6 <- mule %>%
  filter(ID == "MD14F0008_33876_19") %>%
  mutate(vit_2states = factor(vit_2states)) %>%
  ggplot(aes(x = x, y = y, color = vit_2states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_6 | i3_6

# MD14F0008_33876_24

i3_7 <- mule %>%
  filter(ID == "MD14F0008_33876_24") %>%
  mutate(vit_3states = factor(vit_3states)) %>%
  ggplot(aes(x = x, y = y, color = vit_3states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_7 <- mule %>%
  filter(ID == "MD14F0008_33876_24") %>%
  mutate(vit_2states = factor(vit_2states)) %>%
  ggplot(aes(x = x, y = y, color = vit_2states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_7 | i3_7

# MD14F0009_33870_103

i3_8 <- mule %>%
  filter(ID == "MD14F0009_33870_103") %>%
  mutate(vit_3states = factor(vit_3states)) %>%
  ggplot(aes(x = x, y = y, color = vit_3states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_8 <- mule %>%
  filter(ID == "MD14F0009_33870_103") %>%
  mutate(vit_2states = factor(vit_2states)) %>%
  ggplot(aes(x = x, y = y, color = vit_2states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_8 | i3_8

# MD14F0009_33870_104

i3_9 <- mule %>%
  filter(ID == "MD14F0009_33870_104") %>%
  mutate(vit_3states = factor(vit_3states)) %>%
  ggplot(aes(x = x, y = y, color = vit_3states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_9 <- mule %>%
  filter(ID == "MD14F0009_33870_104") %>%
  mutate(vit_2states = factor(vit_2states)) %>%
  ggplot(aes(x = x, y = y, color = vit_2states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_9 | i3_9

# MD14F0009_33870_3

i3_10 <- mule %>%
  filter(ID == "MD14F0009_33870_3") %>%
  mutate(vit_3states = factor(vit_3states)) %>%
  ggplot(aes(x = x, y = y, color = vit_3states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_10 <- mule %>%
  filter(ID == "MD14F0009_33870_3") %>%
  mutate(vit_2states = factor(vit_2states)) %>%
  ggplot(aes(x = x, y = y, color = vit_2states, group = NA)) +
  geom_path() +
  geom_point() +
  theme_bw() +
  theme(legend.position = "none")

i2_10 | i3_10

# Patchwork

(i2_1 | i3_1) /
  (i2_2 | i3_2) /
  (i2_3 | i3_3) /
  (i2_4 | i3_4) /
  (i2_5 | i3_5) /
  (i2_6 | i3_6) /
  (i2_7 | i3_7) /
  (i2_8 | i3_8) /
  (i2_9 | i3_9) /
  (i2_10 | i3_10)

ggsave("output/example-tracks_2-vs-3-states_2023-05-10.tiff",
       width = 8, height = 30, dpi = 300, compression = "lzw")

# Save ####

# Reload individual info back in

path <- "../../../MuleDeerUtah/CleanedFiles/"

# Load and format original data

la_sal <- read.csv(paste0(path, "LaSal_Newest_Cleaned.csv")) %>%
  mutate(study_area = "La Sal",
         mortality = case_when(
           mortality == "true" ~ TRUE,
           mortality == "false" ~ FALSE
         )) %>%
  rename(currentCohort = currentCoh,
         projectName = projectNam,
         captureUnit = captureUni,
         captureSubUnit = captureSub,
         realCaptureArea = realCaptur,
         uniqueID_Range = uniqueID_R) %>%
  dplyr::select(-uniqueID_1, -Used, -dateYearAn) %>%
  mutate(dateYearAndJulian = ymd_hms(dateYearAndJulian, tz = "America/Denver"))
north_slope <- read.csv(paste0(path, "NorthSlope_Newest_Cleaned.csv")) %>%
  mutate(study_area = "North Slope",
         mortality = case_when(
           mortality == "true" ~ TRUE,
           mortality == "false" ~ FALSE
         )) %>%
  rename(currentCohort = currentCoh,
         projectName = projectNam,
         captureUnit = captureUni,
         captureSubUnit = captureSub,
         realCaptureArea = realCaptur,
         uniqueID_Range = uniqueID_R) %>%
  dplyr::select(-uniqueID_1, -Used, -lat_new, -long_new, -dateYearAn) %>%
  mutate(dateYearAndJulian = ymd_hms(dateYearAndJulian, tz = "America/Denver"))
san_juan <- read.csv(paste0(path, "SanJuan_Newest_Cleaned.csv")) %>%
  mutate(study_area = "San Juan",
         mortality = case_when(
           mortality == "true" ~ TRUE,
           mortality == "false" ~ FALSE
         )) %>%
  rename(collarObje = collarObjectId) %>%
  dplyr::select(-lat_new, -long_new) %>%
  mutate(dateYearAndJulian = ymd_hms(dateYearAndJulian, tz = "America/Denver"))
south_manti <- read.csv(paste0(path, "SouthManti_Newest_Cleaned.csv")) %>%
  mutate(study_area = "South Manti") %>%
  rename(currentCohort = currentCoh,
         projectName = projectNam,
         captureUnit = captureUni,
         captureSubUnit = captureSub,
         realCaptureArea = realCaptur,
         uniqueID_Range = uniqueID_R) %>%
  dplyr::select(-Used, -dateYearAn) %>%
  mutate(dateYearAndJulian = mdy_hm(dateYearAndJulian, tz = "America/Denver"))
south_slope <- read.csv(paste0(path, "SouthSlope_Newest_Cleaned.csv")) %>%
  mutate(study_area = "South Slope",
         mortality = case_when(
           mortality == "true" ~ TRUE,
           mortality == "false" ~ FALSE
         )) %>%
  rename(collarObje = collarObjectId) %>%
  dplyr::select(-lat_new, -long_new) %>%
  mutate(dateYearAndJulian = ymd_hms(dateYearAndJulian, tz = "America/Denver"))

# Merge

mule_info <- bind_rows(la_sal, north_slope, san_juan, south_manti, south_slope) %>%
  rename(timestamp = dateYearAndJulian) %>%
  dplyr::select(-projectName, -uniqueID_1, -uniqueID_Range_Year, -X, -da) %>%
  # The same individual was sometimes re-collared with a new collar programmed
  # with a different fix rate (e.g., MD15F0016).
  # Treat these as separate deployments
  mutate(deploy_ID = paste0(uniqueID, "_", collarID)) %>%
  dplyr::select(uniqueID, collarID, deploy_ID, sex, captureUnit, captureSubUnit,
                birthYear, study_area) %>%
  distinct()

# Join info back into data
mule <- mule %>%
  left_join(mule_info, by = "deploy_ID") %>%
  relocate(uniqueID:study_area, .after = ID)

# Get rid of unnecessary columns
mule <- mule %>%
  dplyr::select(-c(nu.x:speed))

# Save
saveRDS(mule, "output/mule-deer_regularized-12h_CTMM_with-viterbi.rds")
