# testing_data <- readRDS("output/mig_spring_100rsteps_with-dyn-covs.rds") %>%
#   unnest(rsteps) %>%
#   filter(is.na(ndvi_end)) %>%
#   slice(2:10) %>%
#   nest(rsteps = -deploy_ID)
#
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
# attach_ndvi(steps, ndvi) %>%
#   View()
#
# steps <- steps %>%
#   select(-ndvi_start, -ndvi_end)

# This function finds the right NDVI layer to use
# depending on the date of the step.
# It uses the closest date preceding the current timestamp.
# date is the timestamp of the step
# ndvi_times is the time column of a SpatRast
find_ndvi <- function(date, ndvi_times) {
  # Only handle one step at a time
  if(length(date) > 1) {
    stop("'date' must be length 1.")
  }
  # Calculate difference between step date and all available NDVI dates
  d <- as.numeric(difftime(date, ndvi_times, units = "days"))
  # Find smallest non-negative difference
  pos <- d[which(d >= 0)]
  min_d <- min(pos)
  # Return the closest preceding NDVI date
  return(ndvi_times[which(d == min_d)])
}

# This function extracts NDVI values from the chosen layer
# steps is the rsteps object
# ndvi is a SpatRast
attach_ndvi <- function(stps, ndvi) {
  # Start date
  stps$date1_ <- as.Date(stps$t1_)
  # End date
  stps$date2_ <- as.Date(stps$t2_)

  # Attach NDVI for start dates
  stps2 <- lapply(split(stps, stps$date1_), function(dd) {
    # Get NDVI date
    ndvi_date <- find_ndvi(dd$date1_[1], time(ndvi))
    # Get that layer
    N <- ndvi[[which(time(ndvi) == ndvi_date)]]
    names(N) <- "ndvi_start"
    # Attach to steps
    res <- extract_covariates(dd, N, where = "start")
    # Return
    return(res)
  }) %>%
    bind_rows()

  # Attach NDVI for end dates
  stps3 <- lapply(split(stps2, stps2$date2_), function(dd) {
    # Get NDVI date
    ndvi_date <- find_ndvi(dd$date2_[1], time(ndvi))
    # Get that layer
    N <- ndvi[[which(time(ndvi) == ndvi_date)]]
    names(N) <- "ndvi_end"
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
