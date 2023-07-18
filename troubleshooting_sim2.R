object = x; compensate.movement = FALSE
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
    phi <- movement_kernel1(xy, object$sl_, object$ta_)
    w <- w + phi - log(xy$sl_)
  }
  w <- exp(w - mean(w[is.finite(w)], na.rm = TRUE))
  w[!is.finite(w)] <- 0
  w
}
