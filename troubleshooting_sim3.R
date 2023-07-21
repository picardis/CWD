mig_spring_amt <- readRDS("output/mig_spring_100rsteps_with-covs.rds")
ndvi <- rast("output/processed_layers/NDVI.tiff")

x <- mig_spring_amt$rsteps[[1]]
x <- x %>%
  filter(as_date(t2_) == ymd("2019-04-27"))

res1 <- x %>%
  extract_covariates_var_time(ndvi,
                              where = "both",
                              when = "before",
                              max_time = days(16)) %>%
  rename(ndvi_start = time_var_covar_start,
         ndvi_end = time_var_covar_end)

source("FUN_ndvi.R")

res2 <- x %>%
  attach_ndvi(ndvi)

find_ndvi(ymd("2019-04-27"), time(ndvi))

time(ndvi$MOD13Q1_NDVI_2019_113)

res3 <- x
res3$ndvi_end <- extract(ndvi$MOD13Q1_NDVI_2019_113, x[, c("x2_", "y2_")])[, 2]

res1_focal <- res1 %>%
  select(x2_, y2_, t2_, land_cover_end, ndvi_end)
res2_focal <- res2 %>%
  select(x2_, y2_, t2_, land_cover_end, ndvi_end)
res3_focal <- res3 %>%
  select(x2_, y2_, t2_, land_cover_end, ndvi_end)

comparison <- res1_focal %>%
  left_join(res2_focal, by = c("x2_", "y2_", "t2_", "land_cover_end")) %>%
  left_join(res3_focal, by = c("x2_", "y2_", "t2_", "land_cover_end"))

View(comparison)
