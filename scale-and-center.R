# Function to scale and center covariates
scale_data <- function(dat, means, sds){
  res <- dat
  for (i in 1:ncol(means)) {
    nm <- names(means)[i]
    nm_sc <- paste0(nm, "_sc")
    res[[nm_sc]] <- (dat[[nm]] - means[[i]])/sds[[i]]
  }
  return(res)
}
