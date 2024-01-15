ssf_formula <- function (formula) {
  rhs <- strsplit(as.character(formula)[3], "\\+")[[1]]
  rhs <- rhs[-grep("strata", rhs)]
  stats::as.formula(paste("~ 0 +", paste(rhs, collapse = "+")))
}
