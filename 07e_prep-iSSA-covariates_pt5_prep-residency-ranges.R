# Load packages ####

library(tidyverse)
library(terra)
library(amt)
library(lubridate)
library(sf)

# Load data ####

mig_spring_amt <- readRDS("output/mig_spring_100rsteps_with-dyn-covs_2024-06-19.rds")
mig_fall_amt <- readRDS("output/mig_fall_100rsteps_with-dyn-covs_2024-06-19.rds")

mule2 <- readRDS("output/mule-deer_2h_with-csu.rds")
mule12 <- readRDS("output/mule-deer_12h_with-csu.rds")

spring_mig_data <- readRDS("output/mig-data_spring_2h-res.rds")
fall_mig_data <- readRDS("output/mig-data_fall_2h-res.rds")

# Set seed ####

set.seed(22)

# Population ranges ####

# Use 12h individuals which weren't included in iSSA
pop_data <- mule12 %>%
  filter(!deploy_ID %in% c(mig_spring_amt$deploy_ID,
                           mig_fall_amt$deploy_ID))

# Save these as validation data for later
saveRDS(pop_data, "output/validation_data.rds")

# Residency data
summ_res <- pop_data %>%
  filter(month(t_) %in% 7:9 &
           vit_3states %in% 1:2)
wint_res <- pop_data %>%
  filter(month(t_) %in% 1:3 &
           vit_3states %in% 1:2)

## WHY WOULD THIS MATTER?
# # Take only deploy_ID that have data for both winter and summer
all_ids <- unique(c(summ_res$deploy_ID, wint_res$deploy_ID))
# not_in_wint <- setdiff(unique(summ_res$deploy_ID), unique(wint_res$deploy_ID))
# not_in_summ <- setdiff(unique(wint_res$deploy_ID), unique(summ_res$deploy_ID))
# focal_ids <- all_ids[!all_ids %in% c(not_in_wint, not_in_summ)]

n_by_unit_focal <- pop_data %>%
  as.data.frame() %>%
  #filter(deploy_ID %in% focal_ids) %>%
  filter(deploy_ID %in% all_ids & !is.na(csu)) %>%
  dplyr::select(deploy_ID, csu) %>%
  distinct() %>%
  group_by(csu) %>%
  tally() %>%
  rename(n_both_seasons = n)

focal_known_cap <- pop_data %>%
  as.data.frame() %>%
  # replaced focal with all
  filter(deploy_ID %in% all_ids & !is.na(csu)) %>%
  pull(deploy_ID) %>%
  unique()

summ_res <- summ_res %>%
  filter(deploy_ID %in% focal_known_cap)
wint_res <- wint_res %>%
  filter(deploy_ID %in% focal_known_cap)

# ggplot(summ_res, aes(x = x, y = y, color = csu)) +
#   geom_point()
# ggplot(wint_res, aes(x = x, y = y, color = csu)) +
#   geom_point()

# summ_res %>%
#   mutate(season = "Summer") %>%
#   bind_rows(wint_res) %>%
#   mutate(season = case_when(
#     is.na(season) ~ "Winter",
#     TRUE ~ season
#   )) %>%
#   ggplot(aes(x = x, y = y, color = season)) +
#   geom_point(alpha = 0.1)

# I need to have residency data for at least one animal in each CSU that's in
# the iSSA data. I need one summer for the spring iSSA CSUs and one winter for
# the fall iSSA CSUs. Check:
need_summer <- unique(unlist(lapply(1:nrow(mig_spring_amt), FUN = function(x) {
  unique(mig_spring_amt$rsteps[[x]]$csu)})))
need_winter <- unique(unlist(lapply(1:nrow(mig_fall_amt), FUN = function(x) {
  unique(mig_fall_amt$rsteps[[x]]$csu)})))
needed_csu <- data.frame(csu = sort(unique(c(need_summer, need_winter)))) %>%
  mutate(needed_summer = csu %in% need_summer,
         needed_winter = csu %in% need_winter) %>%
  # Do I have what I need from the 12h data?
  mutate(got_summer = csu %in% summ_res$csu,
         got_winter = csu %in% wint_res$csu)

# For the CSUs that I don't have info for in the 12h data, I can subset a few
# individuals from the 2h data.
miss_summer <- needed_csu %>%
  filter(needed_summer == TRUE & got_summer == FALSE) %>%
  pull(csu)
miss_winter <- needed_csu %>%
  filter(needed_winter == TRUE & got_winter == FALSE) %>%
  pull(csu)

# Get residency data from 2h data
mig_data_2h <- spring_mig_data %>%
  bind_rows(fall_mig_data) %>%
  mutate(mig = TRUE)

res_data_2h <- mule2 %>%
  left_join(mig_data_2h) %>%
  filter(is.na(mig)) %>%
  dplyr::select(-mig)

# How many individuals do I have for each of those in the 2h data?
n_summ_2h <- res_data_2h %>%
  filter(csu %in% miss_summer,
         month(t_) %in% 7:9) %>%
  select(csu, deploy_ID) %>%
  distinct() %>%
  group_by(csu) %>%
  tally()

n_wint_2h <- res_data_2h %>%
  filter(csu %in% miss_winter,
         month(t_) %in% 1:3) %>%
  select(csu, deploy_ID) %>%
  distinct() %>%
  group_by(csu) %>%
  tally()

# Randomly sample 2 of them per CSU

summ_options <- res_data_2h %>%
  filter(csu %in% miss_summer,
         month(t_) %in% 7:9) %>%
  select(csu, deploy_ID) %>%
  distinct()

summ_samp <- c()
for (u in 1:length(unique(summ_options$csu))) {
  choices <- summ_options[summ_options$csu == unique(summ_options$csu)[u], ]
  if (nrow(choices) > 1) {
    summ_samp <- c(summ_samp,
                   sample(choices$deploy_ID, size = 2, replace = FALSE))
  } else {
    summ_samp[u] <- choices$deploy_ID
  }
}

wint_options <- res_data_2h %>%
  filter(csu %in% miss_winter,
         month(t_) %in% 7:9) %>%
  select(csu, deploy_ID) %>%
  distinct()

wint_samp <- c()
for (u in 1:length(unique(wint_options$csu))) {
  choices <- wint_options[wint_options$csu == unique(wint_options$csu)[u], ]
  if (nrow(choices) > 1) {
    wint_samp <- c(wint_samp,
                   sample(choices$deploy_ID, size = 2, replace = FALSE))
  } else {
    wint_samp[u] <- choices$deploy_ID
  }
}

# Get them
summ_2h <- res_data_2h %>%
  filter(deploy_ID %in% summ_samp,
         month(t_) %in% 7:9)

wint_2h <- res_data_2h %>%
  filter(deploy_ID %in% wint_samp,
         month(t_) %in% 1:3)

# Merge all residency data together now
summ_res <- summ_2h %>%
  rename(x = x_, y = y_) %>%
  select(deploy_ID, x, y, csu) %>%
  bind_rows(summ_res[, c("deploy_ID", "x", "y", "csu")])

wint_res <- wint_2h %>%
  rename(x = x_, y = y_) %>%
  select(deploy_ID, x, y, csu) %>%
  bind_rows(wint_res[, c("deploy_ID", "x", "y", "csu")])

# Make table of how many individuals per CSU
summ_res %>%
  select(deploy_ID, csu) %>%
  distinct() %>%
  group_by(csu) %>%
  tally() %>%
  as.data.frame()

wint_res %>%
  select(deploy_ID, csu) %>%
  distinct() %>%
  group_by(csu) %>%
  tally() %>%
  as.data.frame()

saveRDS(summ_res, "output/residency_data_summer.rds")
saveRDS(wint_res, "output/residency_data_winter.rds")

# Save new version of iSSA data with the sampled individuals removed ####

mig_spring_amt <- mig_spring_amt %>%
  mutate(csu = lapply(rsteps, FUN = function(x) {x$csu[1]}))
mig_spring_amt <- mig_spring_amt %>%
  filter(!deploy_ID %in% summ_2h$deploy_ID &
           !is.na(csu))

mig_fall_amt <- mig_fall_amt %>%
  mutate(csu = lapply(rsteps, FUN = function(x) {x$csu[1]}))
mig_fall_amt <- mig_fall_amt %>%
  filter(!deploy_ID %in% wint_2h$deploy_ID &
           !is.na(csu))

saveRDS(mig_spring_amt, "output/mig_spring_100rsteps_with-dyn-covs_pranges-sample_2024-06-19.rds")
saveRDS(mig_fall_amt, "output/mig_fall_100rsteps_with-dyn-covs_pranges-sample_2024-06-19.rds")
