# Load packages ####

library(tidyverse)
library(sf)
library(lubridate)
library(amt)
library(raster)

# Save? ####

save <- FALSE

# Path to file folder ####

path <- "../../../MuleDeerUtah/CleanedFiles/"

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

st_write(mule_sf, "output/mule-deer-raw.shp")

# Add UTM coordinates ####

mule_utm <- st_transform(mule_sf, crs = 32612)

mule <- mule %>%
  cbind(st_coordinates(mule_utm)) %>%
  rename(utm_x = X, utm_y = Y)

rm(mule_sf)
rm(mule_utm)

# Check temporal resolution ####

# Check how far off points are from the expected time
mule %>%
  mutate(exp_timestamp = round_date(timestamp, unit = "hour"),
         time_error = abs(as.numeric(
           difftime(timestamp, exp_timestamp, unit = "secs")
         ))) %>%
  pull(time_error) %>%
  table() # the vast majority are within 5 minutes (900 secs)
# but there are some that are up to 30 minutes away

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

# Temporal resolution
mule_round_timestamp %>%
  group_by(deploy_ID) %>%
  arrange(timestamp) %>%
  mutate(lag = as.numeric(difftime(timestamp, lag(timestamp), unit = "hours"))) %>%
  pull(lag) %>%
  table() %>%
  sort() # data are mostly at 2-hour intervals

# Were all collars programmed with a 2-hour lag?
mule_round_timestamp %>%
  group_by(deploy_ID) %>%
  arrange(timestamp) %>%
  mutate(lag = as.numeric(difftime(timestamp, lag(timestamp), unit = "hours"))) %>%
  group_by(deploy_ID, lag) %>%
  tally() %>%
  arrange(deploy_ID, desc(n)) %>%
  slice(1) %>%
  pull(lag) %>%
  table()
# It looks like 475 collars were programmed to collect data every 2 hours,
# all others at every 12 or more. Is this right? Check with Kezia.

# Assuming that is right, keep only deployments that have 2-hour data.
keepers <- mule_round_timestamp %>%
  group_by(deploy_ID) %>%
  mutate(lag = as.numeric(difftime(timestamp, lag(timestamp), unit = "hours"))) %>%
  group_by(deploy_ID, lag) %>%
  tally() %>%
  arrange(deploy_ID, desc(n)) %>%
  slice(1) %>%
  filter(lag == 2) %>%
  pull(deploy_ID)

mule_2h <- mule_round_timestamp %>%
  filter(deploy_ID %in% keepers)

# Remove duplicates ####

sum(duplicated(mule_2h[, c("deploy_ID", "timestamp")]))

mule_2h <- mule_2h[!duplicated(mule_2h[, c("deploy_ID", "timestamp")]), ]

# Speed filter ####

mule_sdr <- mule_2h %>%
  nest(trk = -deploy_ID) %>%
  mutate(trk = lapply(trk, function(x) {
    x %>%
      make_track(.x = utm_x, .y = utm_y, .t = timestamp) %>%
      mutate(sl_ = step_lengths(.))
  }))

sdr_values <- do.call("c", lapply(mule_sdr$trk, sdr))
hist(sdr_values, breaks = 5000, xlim = c(0, 2000))

mule_2h <- mule_2h %>%
  nest(trk = -deploy_ID) %>%
  mutate(trk = lapply(trk, function(x) {
    x %>%
      make_track(.x = utm_x, .y = utm_y, .t = timestamp, all_cols = TRUE) %>%
      mutate(sl_ = step_lengths(.)) %>%
      flag_fast_steps(delta = 1000)
  })) %>%
 unnest(cols = trk) %>%
 filter(!fast_step_) %>%
 mutate(uniqueID = stringr::word(deploy_ID, 1, 1, "_"))

# Summaries ####

# Number of individuals
length(unique(mule_2h$uniqueID)) # 468

# Number of deployments
length(unique(mule_2h$deploy_ID)) # 475

# Years
sort(unique(year(mule_2h$t_))) # from 2017 to 2021

# Data points per deployment
pts_per_dep <- mule_2h %>%
  group_by(deploy_ID) %>%
  tally() %>%
  arrange(n)

hist(pts_per_dep$n)

# Drop deployments with few data points ####

# Could be more aggressive here; remove anyone with < 50 points
keepers <- pts_per_dep %>%
  filter(n >= 50) %>%
  pull(deploy_ID)

mule_2h <- mule_2h %>%
  filter(deploy_ID %in% keepers)

# Number of individuals
length(unique(mule_2h$uniqueID)) # 462

# Number of deployments
length(unique(mule_2h$deploy_ID)) # 469

# Introduce NAs at expected timestamps when fix was not taken ####

# I am using track_resample because it will assign a new burst number
# incrementally each time there is a gap in the data longer than 2 hours.
# However, I have to deal with the fact that the schedule sometimes shifts by an
# odd number of hours. I want to introduce NAs for as long as the individual
# stays on either odd or even hours, and treat it as a different burst when it
# switches from odd to even or viceversa.

# Identify places where there is a gap > 2 hours. Each becomes a burst
mule_resample <- mule_2h %>%
  nest(trk = -deploy_ID) %>%
  mutate(trk = lapply(trk, function(x) {
    x %>%
      make_track(.x = x_, .y = y_, .t = t_) %>%
      track_resample(rate = hours(2))
  })) %>%
  unnest(cols = trk)

# Start and end timestamps of each burst
ranges <- mule_resample %>%
  mutate(burst_ID = paste0(deploy_ID, "_", burst_)) %>%
  group_by(burst_ID, deploy_ID) %>%
  summarize(start = min(t_),
            end = max(t_)) %>%
  arrange(deploy_ID, start)

# Identify places where schedule switches by an odd number of hours between
# consecutive bursts of the same individual.
ranges <- ranges %>%
  group_by(deploy_ID) %>%
  mutate(gap = as.numeric(difftime(start, lag(end), unit = "hours")),
         odd = gap %% 2)

# Now assign a streak ID that only changes if the schedule switches by an odd
# number of hours for a single deployment

# Each deployment starts with streak 1
first_streak_start <- ranges %>%
  group_by(deploy_ID) %>%
  slice(1) %>%
  mutate(streak = 1)

ranges <- ranges %>%
  left_join(first_streak_start)

# Incrementally increase streak number when schedule switches by odd number of
# hours
for (i in 1:nrow(ranges)) {
  print(i)
  if (is.na(ranges$streak[i])) {
    ranges$streak[i] <- ranges$streak[i - 1] + ranges$odd[i]
}}

# Create streak ID
ranges <- ranges %>%
  mutate(streak_ID = paste0(deploy_ID, "_", streak))

# Get rid of streaks composed of only one fix
ranges <- ranges %>%
  filter(start != end)

# Now I'm ready to introduce NAs within each streak (not between streaks)

# Create empty data frame to store results
mule_nas <- rep(NA, ncol(mule_2h)) %>%
  as.matrix() %>%
  t() %>%
  as.data.frame()

names(mule_nas) <- names(mule_2h)

mule_nas <- cbind.data.frame(streak_ID = NA,
                             mule_nas)

mule_nas$t_ <- ymd_hms(mule_nas$t_, tz = "America/Denver")

# Loop over streak ID ranges
for (i in 1:length(unique(ranges$streak_ID))) {

  print(i)

  who <- unique(ranges$streak_ID)[i]

  rng <- ranges %>%
    filter(streak_ID == who)

  ts <- ymd_hms(NA, tz = "America/Denver")

  if(nrow(rng) > 1) {

    for (j in 1:(nrow(rng) - 1)) {

      ts <- c(ts,
              seq(rng[j, ]$end, rng[j + 1, ]$start, by = "2 hours"))

    }

    for (k in 1:nrow(rng)) {

      ts <- c(ts,
              seq(rng[k, ]$start, rng[k, ]$end, by = "2 hours"))

    }


  } else {

    ts <- c(ts,
            seq(rng$start, rng$end, by = "2 hours"))

  }

  ts <- ts[-1]

  all_ts <- sort(unique(ts))

  coord <- mule_2h %>%
    filter(deploy_ID == stringr::word(who, 1, 2, "_") &
             t_ >= min(rng$start) &
             t_ <= max(rng$end)) %>%
    mutate(streak_ID = who)

  res <- data.frame(streak_ID = who,
                    t_ = all_ts) %>%
    left_join(coord, by = c("streak_ID", "t_"))

  mule_nas <- rbind(mule_nas, res)

}

mule_nas <- mule_nas[-1, ]

# Save ####

if (save) {saveRDS(mule_nas, "output/mule-deer_regularized-2h.rds")}

# Checks ####

# Check that every lag between locations is 2 hours
# (excluding the first location of each burst which is necessarily NA)
mule_nas %>%
  group_by(streak_ID) %>%
  arrange(t_) %>%
  mutate(lag_ = t_ - lag(t_)) %>%
  group_by(streak_ID) %>%
  slice(-1) %>%
  ungroup() %>%
  group_by(lag_) %>%
  tally()

# Remove streaks of missing fixes > 1 day ####

# For each streak ID and date, do we have at least 1 fix on that day?
test <- mule_nas %>%
  mutate(date = as_date(t_),
         fix = ifelse(is.na(x_), FALSE, TRUE)) %>%
  group_by(streak_ID, date, fix) %>%
  tally() %>%
  ungroup() %>%
  mutate(have_data = case_when(
    !fix & n == 12 ~ FALSE,
    TRUE ~ TRUE
  )) %>%
  distinct(streak_ID, date, have_data)

# Flag first consecutive run for each streak ID
first_run <- test %>%
  group_by(streak_ID) %>%
  arrange(date) %>%
  slice(1) %>%
  mutate(run = 1)

test <- test %>%
  left_join(first_run)

for (r in 2:nrow(test)) {

  print(r)

  test[r, ]$run <- case_when(
    test[r, ]$run == 1 ~ test[r, ]$run,
    test[r, ]$have_data & test[r - 1, ]$have_data ~ test[r - 1, ]$run,
    !test[r, ]$have_data ~ test[r - 1, ]$run,
    test[r, ]$have_data & !test[r - 1, ]$have_data ~ test[r - 1, ]$run + 1
  )

  }

mule_nas <- mule_nas %>%
  mutate(date = as_date(t_)) %>%
  left_join(test, by = c("streak_ID", "date"))

# Assign burst by combining streak ID and run
mule_runs <- mule_nas %>%
  mutate(burst = paste0(streak_ID, "_", run)) %>%
  # Get rid of days without data
  filter(have_data)

# Make sure the individual info are associated to the missing locations too
ind_info <- mule %>%
  dplyr::select(uniqueID, species, sex, captureUnit, uniqueID_Range) %>%
  distinct() %>%
  filter(uniqueID %in% unique(mule_runs$uniqueID)) %>%
  # this individual is duplicated because one time the sex is F, one time it is
  # unknown. Get rid of the unknown
  slice(-308)

mule_final <- mule_runs %>%
  dplyr::select(uniqueID, burst, deploy_ID, t_, x_, y_,
         latitude, longitude,
         mortality, currentAge, currentCohort) %>%
  left_join(ind_info)

# Check
nrow(mule_runs) == nrow(mule_final)

# Get rid of short streaks ####

# Count only points where coordinates are not NA
keepers <- mule_final %>%
  filter(!is.na(x_) & !is.na(y_)) %>%
  group_by(burst) %>%
  tally() %>%
  arrange(n) %>%
  filter(n >= 50) %>%
  pull(burst)

mule_final <- mule_final %>%
  filter(burst %in% keepers)

# Save ####

if (save) {saveRDS(mule_final, "output/mule-deer_regularized-2h_with-NAs.rds")}

# Checks ####

# Check that every lag between locations is 2 hours
# (excluding the first location of each burst which is necessarily NA)
mule_final %>%
  group_by(burst) %>%
  arrange(t_) %>%
  mutate(lag_ = t_ - lag(t_)) %>%
  group_by(burst) %>%
  slice(-1) %>%
  ungroup() %>%
  group_by(lag_) %>%
  tally()

# Map KDE + CWD occurrences ####

# Load data on CWD occurrences (2021-2022)
cwd <- readxl::read_xlsx("input/CWD_June2021_Present.xlsx")

# Load UT boundary
utah <- sf::read_sf("../../../Other People/Ben Crabb/From Ben/clean_data/boundaries/Utah_NAD83_12N.shp")

# Create template raster
template <- raster(ext = extent(utah), res = 1000, crs = crs(utah))

# Subsample to one centroid per individual and format for amt
mule_amt <- mule %>%
  group_by(deploy_ID) %>%
  summarize(centroid_x = mean(utm_x), centroid_y = mean(utm_y)) %>%
  as.data.frame() %>%
  make_track(centroid_x, centroid_y)

# Calculate KDE

hr_all <- hr_kde(mule_amt, trast = template, levels = c(0.90),
                 h = c(20000, 20000))

# Crop to Utah extent

utah_ud <- crop(hr_all$ud, utah)

# Create dummy raster to say whether a point is inside or outside UT
dummy <- rasterize(utah, utah_ud)
values(utah_ud) <- ifelse(is.na(values(dummy)), NA, values(utah_ud))

# Take CWD positive cases
cwd_pos <- cwd %>%
  filter(Results %in% c("Positive") & Species == "Deer") %>%
  st_as_sf(coords = c("GeoLong", "GeoLat"), crs = 4326) %>%
  st_transform(crs = 32612)

# Plot
ggplot() +
  geom_raster(data = as.data.frame(utah_ud, xy = TRUE),
              mapping = aes(x = x, y = y, fill = layer)) +
  geom_sf(data = utah, fill = NA) +
  geom_sf(data = cwd_pos, fill = "chartreuse", shape = 21) +
  scale_fill_viridis_c(option = "magma", na.value = NA,
                       breaks = c(0, 3e-11, 6e-11), labels = c("Low", "Med", "High")) +
  theme_void() +
  labs(fill = "Collared\nmule deer\ndensity", title = "Confirmed CWD cases")

if (save) {ggsave("output/confirmed-CWD_2023-01-16_v2.tiff", compression = "lzw",
       width = 3.5, height = 6, dpi = 400)}
