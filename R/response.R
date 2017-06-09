#'Generator of a response function based on the EWOC PH model
#
#'@param rho a numerical value indicating the true value of the parameter rho.
#'@param mtd a numerical value indicating the true value of the parameter mtd.
#'@param theta a numerical value defining the proportion of expected patients
#'to experience a medically unacceptable, dose-limiting toxicity (DLT) if
#'administered the MTD.
#'@param min_dose a numerical value defining the lower bound of the support of
#'the MTD.
#'@param max_dose a numerical value defining the upper bound of the support of
#'the MTD.
#'@param tau a numerical value defining the period of time for a possible
#'toxicity be observed.
#'@param distribution a character establishing the distribution for the time of
#'events.
#'@param shape a numerical value indicating the true value of the parameter shape.
#'It is only necessary if 'distribution' = "weibul".
#'@export
response_ph <- function(rho, mtd, theta, min_dose, max_dose,
                        tau, distribution, shape = NULL) {

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose,
                         max_dose = max_dose)

  response_sim <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose,
                          max_dose = max_dose)

    if (distribution == "Weibull") {
      if (is.null(shape))
        stop("Weibull distribution requires a shape parameter.")
    } else {
      shape <- 1
    }

    u <- runif(length(dose))
    design <- cbind(1, dose)

    beta <- rep(NA, 2)
    beta[1] <- log(-log(1 - rho)/(tau^shape))
    beta[2] <- log(log(1 - theta)/log(1 - rho))/gamma
    out <- as.numeric((-log(u)*exp(-design%*%beta))^(1/shape))

    return(out)
  }

  out <- response_sim
  return(out)
}


