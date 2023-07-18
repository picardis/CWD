x = sim_mod_spring;
map = rasts;
start = start;
fun = function (xy, map) {
  xy %>%
    extract_covariates(map, where = "both") %>%
    mutate(log_sl_ = log(sl_),
           cos_ta_ = cos(ta_),
           dist_to_roads_end_log = log(dist_to_roads_end + 0.001)) %>%
    scale_data(means = means, sds = sds)
}

covars = NULL; max.dist = amt:::get_max_dist(x); n.control = 1e+06;
n.sample = 1; stochastic = FALSE; compensate.movement = !stochastic;
normalize = TRUE; interpolate = FALSE; as.rast = FALSE;
tolerance.outside = 0
{
  arguments <- as.list(environment())
  checkmate::assert_class(start, "sim_start")
  if (stochastic) {
    xy <- random_steps_simple(start, sl_model = x$sl_, ta_model = x$ta_,
                              n.control = n.control)
  }
  else {
    xy <- amt:::kernel_setup(map, max.dist, start, covars)
  }
  bb.map <- as.vector(terra::ext(map))
  fraction.outside <- mean(xy$x2_ < bb.map["xmin"] | xy$x2_ >
                             bb.map["xmax"] | xy$y2_ < bb.map["ymin"] | xy$y2_ >
                             bb.map["ymax"])
  if (fraction.outside > tolerance.outside) {
    warning(paste0(round(fraction.outside * 100, 3), "% of steps are ending outside the study area but only ",
                   round(tolerance.outside * 100, 3), "% is allowed. ",
                   "Terminating simulations here."))
    return(NULL)
  }
  xy$t1_ <- start$t_
  xy$t2_ <- start$t_ + start$dt
  xy <- fun(xy, map)
  w <- amt:::ssf_weights(xy, x, compensate.movement = compensate.movement)
  r <- if (!as.rast) {
    dplyr::select(xy[sample.int(nrow(xy), size = n.sample,
                                prob = w), ], x_ = x2_, y_ = y2_, t2_)
  }
  else {
    if (stochastic) {
      stop("`as.rast` not implemented for `stochastic = TRUE`")
    }
    else {
      terra::rast(data.frame(xy[, c("x2_", "y2_")], w))
    }
  }
  if (as.rast & normalize) {
    r <- normalize(r)
  }
  res <- list(args = arguments, redistribution.kernel = r)
  class(res) <- c("redistribution_kernel", "list")
  res
}
