# Need this function with "- 1" added to formula for simulation to work properly.
# Need 'redistribution_kernel()' and 'ssf_weights()' here so that the R search
# properly uses this updated 'ssf_formula()' function.
# Also need 'simulate_path()' here so that it calls proper
# 'redistribution_kernel()"
# (But they are unchanged except for 'amt:::' where needed).
ssf_formula <- function (formula)
{
  rhs <- strsplit(as.character(formula)[3], "\\+")[[1]]
  rhs <- rhs[-grep("strata", rhs)]
  stats::as.formula(paste("~", paste(rhs, collapse = "+"), "- 1"))
}

redistribution_kernel <- function (x = make_issf_model(), start = make_start(), map,
                                   fun = function(xy, map) {
                                     extract_covariates(xy, map, where = "both")
                                   }, covars = NULL, max.dist = get_max_dist(x), n.control = 1e+06,
                                   n.sample = 1, stochastic = FALSE, compensate.movement = !stochastic,
                                   normalize = TRUE, interpolate = FALSE, as.rast = FALSE,
                                   tolerance.outside = 0)
{
  arguments <- as.list(environment())
  checkmate::assert_class(start, "sim_start")
  if (stochastic) {
    xy <- amt:::random_steps_simple(start, sl_model = x$sl_, ta_model = x$ta_,
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
  w <- ssf_weights(xy, x, compensate.movement = compensate.movement)
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
    r <- amt:::normalize(r)
  }
  res <- list(args = arguments, redistribution.kernel = r)
  class(res) <- c("redistribution_kernel", "list")
  res
}

ssf_weights <- function (xy, object, compensate.movement = FALSE)
{
  checkmate::assert_class(xy, "data.frame")
  checkmate::assert_class(object, "fit_clogit")
  checkmate::assert_logical(compensate.movement)
  coefs <- coef(object)
  ff <- ssf_formula(object$model$formula)
  newdata <- xy
  attr(newdata, "na.action") <- "na.pass"
  xyz <- stats::model.matrix.default(ff, data = newdata, na.action = stats::na.pass)
  w <- as.matrix(xyz[, names(coefs)]) %*% coefs
  if (compensate.movement) {
    phi <- amt:::movement_kernel1(xy, object$sl_, object$ta_)
    w <- w + phi - log(xy$sl_)
  }
  w <- exp(w - mean(w[is.finite(w)], na.rm = TRUE))
  w[!is.finite(w)] <- 0
  w
}

simulate_path.redistribution_kernel <-
  function (x, n.steps = 100, start = x$args$start, verbose = FALSE, ...)
{
  mod <- x$args
  xy <- tibble(x_ = rep(NA, n.steps + 1), y_ = NA_real_, t_ = start$t_ +
                 start$dt * (0:n.steps), dt = start$dt)
  xy$x_[1] <- start$x_
  xy$y_[1] <- start$y_
  for (i in 1:n.steps) {
    rk <- redistribution_kernel(x = mod$x, start = start,
                                map = mod$map, fun = mod$fun, max.dist = mod$max.dist,
                                n.control = mod$n.control, n.sample = 1, stochastic = mod$stochastic,
                                normalize = TRUE, interpolate = FALSE, as.rast = FALSE,
                                tolerance.outside = mod$tolerance.outside)
    if (is.null(rk)) {
      warning(paste0("Simulation stopped after ", i -
                       1, " time steps, because the animal stepped out of the landscape."))
      return(xy)
    }
    rk <- rk$redistribution.kernel
    new.ta <- atan2(rk$y_[1] - start$y_[1], rk$x_[1] - start$x_[1])
    xy$x_[i + 1] <- rk$x_[1]
    xy$y_[i + 1] <- rk$y_[1]
    start <- make_start(as.numeric(xy[i + 1, c("x_", "y_")]),
                        new.ta, crs = attr(x$args$start, "crs"))
  }
  return(xy)
}

