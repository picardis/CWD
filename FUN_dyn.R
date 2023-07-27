# testing_data <- mig_spring_amt[mig_spring_amt$deploy_ID == "MD17F0171_39994", ]$rsteps[[1]] %>%
#   select(-ndvi_start, -ndvi_end)
#
# # testing_data <- readRDS("output/mig_spring_100rsteps_with-dyn-covs.rds") %>%
# #   unnest(rsteps) %>%
# #   filter(is.na(ndvi_end)) %>%
# #   slice(2:10) %>%
# #   nest(rsteps = -deploy_ID)
#
# test_ndvi <- ndvi$MOD13Q1_NDVI_2019_081
#
# extract(test_ndvi, testing_data$rsteps[[1]][, c("x2_", "y2_")])

# steps <- testing_data$rsteps[[1]]
# class(steps) <- class(old_data$rsteps[[1]])
#
#   steps %>%
#     select(-ndvi_start, -ndvi_end) %>%
#     extract_covariates_var_time(ndvi,
#                               where = "both",
#                               when = "any",
#                               max_time = days(16)) %>%
#   rename(ndvi_start = time_var_covar_start,
#          ndvi_end = time_var_covar_end) %>%
#     View()
#
# attach_dyn(steps, ndvi) %>%
#   View()
#
# steps <- steps %>%
#   select(-ndvi_start, -ndvi_end)

# This function finds the right layer to use from a dynamic raster
# depending on the date of the step.
# It uses the closest date preceding the current timestamp.
# date is the timestamp of the step
# ndvi_times is the time column of a SpatRast
find_layer <- function(date, layer_times) {
  # Only handle one step at a time
  if(length(date) > 1) {
    stop("'date' must be length 1.")
  }
  # Calculate difference between step date and all available dates
  d <- as.numeric(difftime(date, layer_times, units = "days"))
  # Find smallest non-negative difference
  pos <- d[which(d >= 0)]
  min_d <- min(pos)
  # Return the closest preceding date
  return(layer_times[which(d == min_d)])
}

# This function extracts values from the chosen layer
# steps is the rsteps object
# dyn is a SpatRast with a z column (time)
# col_name is a character string indicating the desired column name in output
attach_dyn <- function(stps, dyn, col_name) {
  # Start date
  stps$date1_ <- as.Date(stps$t1_)
  # End date
  stps$date2_ <- as.Date(stps$t2_)

  # Attach values for start dates
  stps2 <- lapply(split(stps, stps$date1_), function(dd) {
    # Get appropriate date
    layer_date <- find_layer(dd$date1_[1], time(dyn))
    # Get that layer
    N <- dyn[[which(time(dyn) == layer_date)]]
    names(N) <- paste0(col_name, "_start")
    # Attach to steps
    res <- extract_covariates(dd, N, where = "start")
    # Return
    return(res)
  }) %>%
    bind_rows()

  # Attach NDVI for end dates
  stps3 <- lapply(split(stps2, stps2$date2_), function(dd) {
    # Get NDVI date
    layer_date <- find_layer(dd$date2_[1], time(dyn))
    # Get that layer
    N <- dyn[[which(time(dyn) == layer_date)]]
    names(N) <- paste0(col_name, "_end")
    # Attach to steps
    res <- extract_covariates(dd, N, where = "end")
    # Return
    return(res)
  }) %>%
    bind_rows() %>%
    select(-c(date1_, date2_))

  # Return
  return(stps3)

}
