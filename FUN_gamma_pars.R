# Convert gamma parameters
gamma_pars <- function(mean = NULL, sd = NULL,
                       rate = NULL, shape = NULL, scale = NULL) {

  if(is.numeric(mean) && is.numeric(sd)) {

    rate <- mean/sd^2
    shape <- (mean/sd)^2

    pars <- c(rate = rate,
              shape = shape)

  } else if(is.numeric(rate) && is.numeric(shape)) {

    mean <- shape/rate
    sd <- sqrt(shape)/rate

    pars <- c(mean = mean,
              sd = sd)

  } else if(is.numeric(scale) && is.numeric(shape)) {

    rate <- 1/scale

    mean <- shape/rate
    sd <- sqrt(shape)/rate

    pars <- c(mean = mean,
              sd = sd)

  }

  return(pars)

}
