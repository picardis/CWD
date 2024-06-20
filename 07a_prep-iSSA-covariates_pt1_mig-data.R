# Load packages ####

library(tidyverse)
library(amt)

# Load data ####

mule2 <- readRDS("output/mule-deer_regularized-2h_with-NAs_2023-07-22.rds")
mule2 <- mule2[!is.na(mule2$deploy_ID), ] # Temporary fix. This needs to be addressed upstream

# Load HMM output ####

mule12 <- readRDS("output/mule-deer_regularized-12h_CTMM_with-viterbi.rds")

# Load new information on capture sites ####

cap_lasal <- read.csv("input/lasals_animals_with_cap_areas_krmmanual.csv") %>%
  mutate(deploy_ID = paste0(uniqueID, "_", collarID)) %>%
  select(deploy_ID, realCaptur)
cap_sjuan <- read.csv("input/sanjuan_animals_with_cap_areas_krmmanual.csv") %>%
  mutate(deploy_ID = paste0(uniqueID, "_", collarID)) %>%
  select(deploy_ID, realCaptur)

cap_all <- bind_rows(cap_lasal, cap_sjuan)

mule2 <- mule2 %>%
  left_join(cap_all)

mule12 <- mule12 %>%
  left_join(cap_all)

# Now create new column that uses realCaptur when available and the old
# capture subunit when realCaptur is not available
mule2 <- mule2 %>%
  mutate(csu = case_when(
    !is.na(realCaptur) ~ realCaptur,
    is.na(realCaptur) ~ captureSubUnit
  )) %>%
  mutate(csu = case_when(
    csu == "Abajo Mtns" ~ "Abajo",
    TRUE ~ csu
  ))

mule12 <- mule12 %>%
  mutate(csu = case_when(
    !is.na(realCaptur) ~ realCaptur,
    is.na(realCaptur) ~ captureSubUnit
  )) %>%
  mutate(csu = case_when(
    csu == "Abajo Mtns" ~ "Abajo",
    TRUE ~ csu
  ))

saveRDS(mule2, "output/mule-deer_2h_with-csu.rds")
saveRDS(mule12, "output/mule-deer_12h_with-csu.rds")

# Get migration intervals ####

mule_mig_bouts <- mule12 %>%
  select(ID, t_, vit_3states) %>%
  arrange(ID, t_) %>%
  mutate(start_and_end = case_when(
    # Single-location bouts
    # Within-individual
    (lead(vit_3states) != 3 &
       lag(vit_3states) != 3 &
       vit_3states == 3) |
      # At the beginning
      (lag(ID) != ID &
         lead(vit_3states) != 3 &
         vit_3states == 3) |
      # At the end
      (lead(ID) != ID &
         lag(vit_3states) != 3 &
         vit_3states == 3) ~ "single",
    # First location of one ID
    lag(ID) != ID & vit_3states == 3 ~ "start",
    # First location of a new migration bout
    lag(ID) == ID &
      lag(vit_3states) != 3 & vit_3states == 3 ~ "start",
    # Last location of a migration bout
    lag(ID) == ID & lead(ID) == ID &
      lead(vit_3states) != 3 & vit_3states == 3 ~ "end",
    # Last location of one ID
    lead(ID) != ID & vit_3states == 3 ~ "end",
    # Else, NA
    TRUE ~ NA_character_
  ))

mule_mig_starts <- mule_mig_bouts %>%
  filter(start_and_end == "start") %>%
  select(ID, start_ = t_)

mule_mig_ends <- mule_mig_bouts %>%
  filter(start_and_end == "end") %>%
  select(ID, end_ = t_)

mig_intervals <- mule_mig_starts %>%
  select(-ID) %>%
  bind_cols(mule_mig_ends) %>%
  select(ID, start_, end_)

table(month(mig_intervals$start_))
# April-May and December are peak migrations
# Could restrict analysis to just these
# Or could do it for all/most
# Split:
# Spring migration = April-June
# Fall migration = October-December

spring_mig <- mig_intervals %>%
  filter(month(start_) %in% c(4:6)) %>%
  mutate(mig_int = start_ %--% end_)

fall_mig <- mig_intervals %>%
  filter(month(start_) %in% c(10:12)) %>%
  mutate(mig_int = start_ %--% end_)

# This is one spot where I lose a bunch of individuals. I go from 596 to
# 387 in the spring and 511 in the fall. Others have migration steps in
# different months.

# Filter data within migration intervals ####

# Spring

# Test
# i <- spring_mig$ID[352]

spring_mig_data <- data.frame()

for (i in unique(spring_mig$ID)) {

  ints <- spring_mig %>%
    filter(ID == i)

  depID <- word(i, 1, 2, "_")

  all_res <- data.frame()

  for (r in 1:nrow(ints)) {

    res <- mule2 %>%
      filter(deploy_ID == depID &
               t_ %within% ints$mig_int[r])

    all_res <- bind_rows(all_res, res)

  }

  spring_mig_data <- bind_rows(spring_mig_data, all_res)

}

# Fall

fall_mig_data <- data.frame()

for (i in unique(fall_mig$ID)) {

  ints <- fall_mig %>%
    filter(ID == i)

  depID <- word(i, 1, 2, "_")

  all_res <- data.frame()

  for (r in 1:nrow(ints)) {

    res <- mule2 %>%
      filter(deploy_ID == depID &
               t_ %within% ints$mig_int[r])

    all_res <- bind_rows(all_res, res)

  }

  fall_mig_data <- bind_rows(fall_mig_data, all_res)

}

# Save ####

saveRDS(spring_mig_data, "output/mig-data_spring_2h-res.rds")
saveRDS(fall_mig_data, "output/mig-data_fall_2h-res.rds")

# Format for use in amt ####

mig_spring_amt <- spring_mig_data %>%
  filter(!is.na(x_) & !is.na(y_)) %>%
  # Make the unit the deployment-year rather than the deployment ID
  mutate(year = year(t_),
         deploy_ID_year = paste0(deploy_ID, "_", year)) %>%
  nest(cols = -deploy_ID_year) %>%
  mutate(deploy_ID = stringr::word(deploy_ID_year, 1, 2, "_")) %>%
  mutate(trk = lapply(cols, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      make_track(.x = x_, .y = y_, .t = t_,
                 all_cols = TRUE) %>%
      track_resample(rate = hours(2))
  })) %>%
  mutate(steps = lapply(trk, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      steps_by_burst(keep_cols = "end")
  })) %>%
  select(deploy_ID_year, deploy_ID, cols, trk, steps)

# Check that all steps are at 2h resolution (FALSE means they all are)
any(unlist(lapply(mig_spring_amt$steps, function(x) {unique(x$dt_)})) != 2)

mig_fall_amt <- fall_mig_data %>%
  filter(!is.na(x_) & !is.na(y_)) %>%
  # Make the unit the deployment-year rather than the deployment ID
  mutate(year = year(t_),
         deploy_ID_year = paste0(deploy_ID, "_", year)) %>%
  nest(cols = -deploy_ID_year) %>%
  mutate(deploy_ID = stringr::word(deploy_ID_year, 1, 2, "_")) %>%
  mutate(trk = lapply(cols, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      make_track(.x = x_, .y = y_, .t = t_,
                 all_cols = TRUE) %>%
      track_resample(rate = hours(2))
  })) %>%
  mutate(steps = lapply(trk, FUN = function (x) {
    x %>%
      arrange(t_) %>%
      steps_by_burst(keep_cols = "end")
  })) %>%
  select(deploy_ID_year, deploy_ID, cols, trk, steps)

# Check that all steps are at 2h resolution (FALSE means they all are)
any(unlist(lapply(mig_fall_amt$steps, function(x) {unique(x$dt_)})) != 2)

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

# Save ####

saveRDS(mig_spring_amt, "output/mig-trk_spring_2h-res.rds")
saveRDS(mig_fall_amt, "output/mig-trk_fall_2h-res.rds")
