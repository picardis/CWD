# The HMM is unable to pick up migration with data at 2-hour resolution.
# Here, I resample data at a coarser time resolution.

# Load packages ####

library(tidyverse)
library(sf)
library(lubridate)
library(amt)
library(momentuHMM)

# Save? ####

save <- FALSE

# Path to file folder ####

path <- "../MuleDeerUtah/CleanedFiles/"

list.files(path)

# Load and format data ####

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

# Ask Kezia:
# 1. what are lat_new and long_new?
# 2. am I correct to assume Mountain time zone or is it UTC?

# Merge ####

mule <- bind_rows(la_sal, north_slope, san_juan, south_manti, south_slope) %>%
  rename(timestamp = dateYearAndJulian) %>%
  dplyr::select(-projectName, -uniqueID_1, -uniqueID_Range_Year, -X, -da) %>%
  # The same individual was sometimes re-collared with a new collar programmed
  # with a different fix rate (e.g., MD15F0016).
  # Treat these as separate deployments
  mutate(deploy_ID = paste0(uniqueID, "_", collarID))

summary(mule)

rm(list = ls()[!ls() %in% c("mule", "path")])

# Transform to spatial ####

mule_sf <- mule %>%
  dplyr::select(uniqueID, deploy_ID, latitude, longitude, timestamp, study_area) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Add UTM coordinates ####

mule_utm <- st_transform(mule_sf, crs = 32612)

mule <- mule %>%
  cbind(st_coordinates(mule_utm)) %>%
  rename(utm_x = X, utm_y = Y)

rm(mule_sf)
rm(mule_utm)

# Round timestamps ####

# Get rid of any fixes taken > 15 min away from the expected timestamp and
# round timestamp to the closest hour
mule_round_timestamp <- mule %>%
  mutate(exp_timestamp = round_date(timestamp, unit = "hour"),
         time_error = abs(as.numeric(
           difftime(timestamp, exp_timestamp, unit = "secs")
         ))) %>%
  filter(time_error <= 900) %>%
  dplyr::select(-timestamp, time_error) %>%
  rename(timestamp = exp_timestamp)

# Keep only deployments that have 12-hour data or finer
keepers <- mule_round_timestamp %>%
  group_by(deploy_ID) %>%
  mutate(lag = as.numeric(difftime(timestamp, lag(timestamp), unit = "hours"))) %>%
  group_by(deploy_ID, lag) %>%
  tally() %>%
  arrange(deploy_ID, desc(n)) %>%
  slice(1) %>%
  filter(lag %in% c(2, 12)) %>%
  pull(deploy_ID)

mule_12h <- mule_round_timestamp %>%
  filter(deploy_ID %in% keepers)

# CTMM ####

mule_12h <- mule_12h %>%
  mutate(ID = deploy_ID) %>%
  mutate(timestamp_utc = with_tz(timestamp, tzone = "UTC"))

set.seed(1)

rm(list = setdiff(ls(), "mule_12h"))

for (i in 1:length(unique(mule_12h$ID))) {

  cat(i, "of", length(unique(mule_12h$ID)), "\n")

  try({sub <- mule_12h[mule_12h$ID == unique(mule_12h$ID)[i], ]

  mule_ctmm <- crawlWrap(obsData = sub,
                          timeStep = "12 hours",
                          coord = c("utm_x", "utm_y"),
                          Time.name = "timestamp_utc")

  ctmm_12h <- mule_ctmm$crwPredict %>%
    as_tibble() %>%
    filter(locType == "p") %>%
    mutate(uncertainty = case_when(
      se.mu.x > 0 | se.mu.y > 0 ~ "Uncertainty",
      TRUE ~ "None"
    ))

  saveRDS(ctmm_12h, paste0("output/ctmm/", unique(mule_12h$ID)[i], "_12h.rds"))

  })

}

# Bind back together
files <- list.files("output/ctmm", full.names = TRUE)

ctmm <- do.call("rbind", lapply(files, function(x) {
  readRDS(x)
  }))

locs <- ctmm %>%
  select(ID:uncertainty)

unc <- data.frame(n = rle(locs$uncertainty)[[1]],
                  value = rle(locs$uncertainty)[[2]]) %>%
  filter(value == "Uncertainty")

# Add flag for data gaps ####

# Start and end date of each deployment
start_end <- locs %>%
  group_by(ID) %>%
  summarize(start = min(as_date(timestamp_utc)),
            end = max(as_date(timestamp_utc)))

# Initialize flag column
locs_flag <- locs %>%
  mutate(data_gap = NA)

# Loop over bursts
for (i in 1:nrow(start_end)) {

  cat(i, "of", nrow(start_end), "\n")
  # Which ID are we on
  who <- start_end$ID[i]
  # Make sequence of days between start and end date
  dseq <- seq(from = start_end$start[i], to = start_end$end[i], by = 1)
  # Grab original data
  orig_dat <- mule_12h %>%
    filter(ID == who)
  # On which days of the sequence did we NOT have data?
  days_without_data <- dseq[which(!dseq %in% as_date(orig_dat$timestamp))]
  # On which days of the sequence did we have data?
  days_with_data <- dseq[which(dseq %in% as_date(orig_dat$timestamp))]
  # Flag days without data with data_gap = TRUE, else FALSE
  locs_flag <- locs_flag %>%
    mutate(data_gap = case_when(
      ID == who & as_date(timestamp_utc) %in% days_without_data ~ TRUE,
      ID == who & as_date(timestamp_utc) %in% days_with_data ~ FALSE,
      TRUE ~ data_gap
    ))
}

# Split bursts at data gaps ####

locs_no_gaps <- locs_flag %>%
  filter(!data_gap)

locs_resample <- locs_no_gaps %>%
  nest(trk = -ID) %>%
  mutate(trk = lapply(trk, function(x) {
    x %>%
      make_track(.x = mu.x, .y = mu.y,
                 .t = timestamp_utc, all_cols = TRUE) %>%
      track_resample(rate = hours(12))
  })) %>%
  unnest(cols = trk) %>%
  rename(deploy_ID = ID) %>%
  mutate(ID = paste0(deploy_ID, "_", burst_))

keepers <- locs_resample %>%
  group_by(ID) %>%
  tally() %>%
  arrange(n) %>%
  filter(n > 9) %>%
  pull(ID)

locs_ctmm_no_gaps <- locs_resample %>%
  filter(ID %in% keepers)

# Summaries ####

# Number of bursts
length(unique(locs_ctmm_no_gaps$ID)) # 6164

# Years
sort(unique(year(locs_ctmm_no_gaps$t_))) # from 2014 to 2021

# Check that every lag between locations is 12 hours
# (excluding the first location of each burst which is necessarily NA)
locs_ctmm_no_gaps %>%
  group_by(ID) %>%
  arrange(t_) %>%
  mutate(lag_ = t_ - lag(t_)) %>%
  group_by(ID) %>%
  slice(-1) %>%
  ungroup() %>%
  group_by(lag_) %>%
  tally()

# Bring timestamps back to Mountain time ####

locs_ctmm_no_gaps <- locs_ctmm_no_gaps %>%
  mutate(t_ = with_tz(t_, tzone = "America/Denver"))

# Save ####

if (save) {saveRDS(locs_ctmm_no_gaps, "output/mule-deer_regularized-12h_CTMM_no-gaps.rds")}

# Insert NA coordinates in data gaps ####

locs_ctmm_nas <- locs_flag %>%
  mutate(x_ = case_when(
    data_gap == TRUE ~ NA_real_,
    data_gap == FALSE ~ mu.x
  ),
  y_ = case_when(
    data_gap == TRUE ~ NA_real_,
    data_gap == FALSE ~ mu.y
  ))

# Bring timestamps back to Mountain time ####

locs_ctmm_nas <- locs_ctmm_nas %>%
  mutate(t_ = with_tz(timestamp_utc, tzone = "America/Denver"))

# Discard short bursts ####

# Count only points where coordinates are not NA
keepers <- locs_ctmm_nas %>%
  filter(!is.na(x_) & !is.na(y_)) %>%
  group_by(ID) %>%
  tally() %>%
  arrange(n) %>%
  filter(n >= 10) %>%
  pull(ID)

locs_ctmm_nas <- locs_ctmm_nas %>%
  filter(ID %in% keepers)

# Save ####

if (save) {saveRDS(locs_ctmm_nas, "output/mule-deer_regularized-12h_CTMM_NA-gaps.rds")}
