xx <- mig_spring_amt$rsteps[[which(mig_spring_amt$deploy_ID == ind)]]

xy <- mig_spring_amt$rsteps[[which(mig_spring_amt$deploy_ID == ind)]][1, ] %>%
  select(x1_:ta_, t1_:dt_) %>%
  extract_covariates(rasts_spring, where = "both") %>%
  attach_dyn(dyn = ndvi, col_name = "ndvi") %>%
  mutate(log_sl_ = log(sl_),
         cos_ta_ = cos(ta_),
         dist_to_roads_end_log = log(dist_to_roads_end + 0.001),
         dist_to_range_end_log = log(dist_to_range_end + 0.001),
         land_cover_end = factor(land_cover_end, levels = c("Forest",
                                                            "Agricultural",
                                                            "Developed",
                                                            "Grassland",
                                                            "Shrubland",
                                                            "Unsuitable")),
         ndvi_start = case_when(
           land_cover_start == "Open Water" ~ 0,
           TRUE ~ ndvi_start
         ),
         ndvi_end = case_when(
           land_cover_end == "Open Water" ~ 0,
           TRUE ~ ndvi_end
         )) %>%
  scale_data(means = means, sds = sds)

object <- issa_spring$issa[[which(issa_spring$deploy_ID == ind)]]

function (xy, object, compensate.movement = FALSE)
{
  checkmate::assert_class(xy, "data.frame")
  checkmate::assert_class(object, "fit_clogit")
  checkmate::assert_logical(compensate.movement)
  coefs <- coef(object)
  ff <- ssf_formula(object$model$formula)
  newdata <- xy
  attr(newdata, "na.action") <- "na.pass"
  xyz <- stats::model.matrix.default(ff, data = newdata,
                                     na.action = stats::na.pass,
                                     contrasts.arg = object$model$contrasts)
  xyz <- xyz[, 2:ncol(xyz), drop = FALSE]
  w <- as.matrix(xyz) %*% coefs
  if (compensate.movement) {
    phi <- movement_kernel1(xy, object$sl_, object$ta_)
    w <- w + phi - log(xy$sl_)
  }
  w <- exp(w - mean(w[is.finite(w)], na.rm = TRUE))
  w[!is.finite(w)] <- 0
  w
}
