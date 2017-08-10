#'Generating a binary response function based on the EWOC basic model
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
#'@export
response_d1basic <- function(rho, mtd, theta, min_dose, max_dose) {

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose,
                         max_dose = max_dose)

  response_sim <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose,
                          max_dose = max_dose)

    beta <- rep(NA, 2)
    beta[1] <- logit(rho)
    beta[2] <- (logit(theta) - logit(rho))/gamma

    design_matrix <- cbind(1, dose)
    lp <- design_matrix %*% beta
    p <- as.numeric(plogis(lp))
    out <- rbinom(n = length(dose), size = 1, prob = p)
    return(out)
  }

  out <- response_sim
  return(out)
}


#'Generating a binary response function based on the EWOC extended model
#
#'@param rho a numerical vector indicating the true value of the parameters
#'rho_0 and rho_1.
#'@param theta a numerical value defining the proportion of expected patients
#'to experience a medically unacceptable, dose-limiting toxicity (DLT) if
#'administered the MTD.
#'@param min_dose a numerical value defining the lower bound of the support of
#'the MTD.
#'@param max_dose a numerical value defining the upper bound of the support of
#'the MTD.
#'@export
response_d1extended <- function(rho, theta, min_dose, max_dose) {

  response_sim <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose,
                          max_dose = max_dose)

    beta <- rep(NA, 2)
    beta[1] <- logit(rho[1])
    beta[2] <- logit(rho[2]) - logit(rho[1])

    design_matrix <- cbind(1, dose)
    lp <- design_matrix %*% beta
    p <- as.numeric(plogis(lp))
    out <- rbinom(n = length(dose), size = 1, prob = p)
    return(out)
  }

  out <- response_sim
  return(out)
}

#'Generating a response function based on the EWOC PH model
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
#'It is only necessary if 'distribution' = "weibull".
#'@export
response_d1ph <- function(rho, mtd, theta, min_dose, max_dose,
                          tau, distribution, shape = NULL) {

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose,
                         max_dose = max_dose)

  response_sim <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose,
                          max_dose = max_dose)

    if (distribution == "weibull") {
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

#'Generating a binary response function based on the EWOC model with a
#'multinomial covariable
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
#'@export
response_d1multinomial <- function(rho, mtd, theta, min_dose, max_dose,
                                   levels_cov) {

  response_sim <- function(dose, cov){

    gamma <- standard_dose(dose = mtd,
                           min_dose = min_dose(cov),
                           max_dose = max_dose(cov))

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose(cov),
                          max_dose = max_dose(cov))

    beta <- rep(NA, length(mtd))
    beta[1] <- logit(rho[1])
    beta[2] <- (logit(theta) - logit(rho[1]))/gamma[1]
    for (i in 3:(np+1)) {
      beta[i] <- logit(theta) - logit(rho[1]) -
        gamma[i-1]*(logit(theta) - beta[1])/gamma[i-2]
    }

    cov <- factor(cov, levels = levels_cov)

    design_matrix <- cbind(1, dose, model.matrix(~ cov))
    lp <- design_matrix %*% beta
    p <- as.numeric(plogis(lp))
    out <- rbinom(n = length(dose), size = 1, prob = p)
    return(out)
  }

  out <- response_sim
  return(out)
}

