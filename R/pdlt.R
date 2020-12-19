#'Generating a probability of DLT function based on the EWOC classical model
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
#'
#'@return A function with dose as an input and a probability based on the
#'logistic regression and parameters as an output.
#'
#'@examples
#'pdlt <- pdlt_d1classical(rho = 0.05, mtd = 60, theta = 0.33,
#'                       min_dose = 20, max_dose = 100)
#'
#'pdlt(20)
#'
#'
#'@export
pdlt_d1classical <- function(rho, mtd, theta, min_dose, max_dose){

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose,
                         max_dose = max_dose)

  pdlt <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose,
                          max_dose = max_dose)

    beta <- rep(NA, 2)
    beta[1] <- logit(rho)
    beta[2] <- (logit(theta) - logit(rho))/gamma

    design_matrix <- cbind(1, dose)
    lp <- design_matrix %*% beta
    out <- as.numeric(plogis(lp))
    return(out)
  }

  return(pdlt)
}

#'Generating a probability of DLT function based on the EWOC extended model
#
#'@param rho a numerical vector indicating the true value of the parameters
#'rho_0 and rho_1.
#'@param min_dose a numerical value defining the lower bound of the support of
#'the MTD.
#'@param max_dose a numerical value defining the upper bound of the support of
#'the MTD.
#'
#'@return A function with dose as an input and a probability based on the
#'logistic regression and parameters as an output.
#'
#'@examples
#'pdlt <- pdlt_d1extended(rho = c(0.05, 0.5),
#'                        min_dose = 10, max_dose = 50)
#'pdlt(20)
#'
#'@export
pdlt_d1extended <- function(rho, min_dose, max_dose){

  pdlt <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose,
                          max_dose = max_dose)

    beta <- rep(NA, 2)
    beta[1] <- logit(rho[1])
    beta[2] <- logit(rho[2]) - logit(rho[1])

    design_matrix <- cbind(1, dose)
    lp <- design_matrix%*%beta
    out <- as.numeric(plogis(lp))
    return(out)
  }

  return(pdlt)
}

#'Generating a probability of DLT function based on the EWOC Proportional Hazards model
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
#'events. It can be defined as 'exponential' or 'weibull'.
#'@param shape a numerical value indicating the true value of the parameter shape.
#'It is only necessary if 'distribution' = "weibull".
#'
#'@return A function with dose as an input and a probability based on the
#'logistic regression and parameters as an output.
#'
#'@examples
#'pdlt <- pdlt_d1ph(rho = 0.05, mtd = 40, theta = 0.33,
#'                  min_dose = 30, max_dose = 50,
#'                  tau = 10, distribution = "exponential")
#'pdlt(40)
#'
#'@export
pdlt_d1ph <- function(rho, mtd, shape = NULL, theta, min_dose, max_dose,
                      tau, distribution) {

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose,
                         max_dose = max_dose)

  pdlt <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose,
                          max_dose = max_dose)

    if (distribution != "weibull")
      shape <- 1

    beta <- rep(NA, 2)
    beta[1] <- log(-log(1 - rho[1])/(tau^shape))
    beta[2] <- log(log(1 - theta)/log(1 - rho[1]))/gamma

    design_matrix <- cbind(1, dose)

    out <- 1 - exp(-(exp(design_matrix%*%beta)*tau)^shape)

    return(out)
  }

  return(pdlt)
}

#'Generating a probability of DLT function based on the EWOC with dichotomous covariate model
#
#'@param rho a numerical value indicating the true value of the parameter rho.
#'@param mtd a numerical vector indicating the true values of the parameters mtd associated with \code{level_cov}.
#'@param theta a numerical value defining the proportion of expected patients
#'to experience a medically unacceptable, dose-limiting toxicity (DLT) if
#'administered the MTD.
#'@param min_dose either a numerical value or a function of the covariate
#'defining the lower bound of the support of the MTD used to standardize the doses.
#'@param max_dose either a numerical value or a function of the covariate
#'defining the upper bound of the support of the MTD used to standardize the doses.
#'@param levels_cov a character vector of the possible values for the covariate.
#'
#'@return A function with dose as an input and a probability based on the
#'logistic regression and parameters as an output.
#'
#'@examples
#'pdlt <- pdlt_d1dicov(rho = 0.05, mtd = c(40,60), theta = 0.33,
#'                     min_dose = 20, max_dose = 100,
#'                     levels_cov = c("A", "B"))
#'
#'pdlt(20, "A")
#'
#'@export
pdlt_d1dicov <- function(rho, mtd, theta, min_dose, max_dose, levels_cov) {

  if (length(mtd) == 1)
    stop("mtd has length 1.")

  if (!is.function(min_dose)){
    min_dose_value <- min_dose

    min_dose <- function(covariate = NULL){
      out <- min_dose_value
      return(out)
    }

    min_dose <- Vectorize(min_dose)
  }

  if (!is.function(max_dose)) {
    max_dose_value <- max_dose

    max_dose <- function(covariate = NULL){
        out <- max_dose_value
        return(out)
    }
    max_dose <- Vectorize(max_dose)
  }

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose(levels_cov),
                         max_dose = max_dose(levels_cov))

  pdlt <- function(dose, cov){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose(cov),
                          max_dose = max_dose(cov))

  }
  return(pdlt)
}

#'@export
pdlt_d1excontinuous <- function(rho, direction,
                                min_dose, max_dose, min_cov, max_cov) {

  pdlt <- function(dose, cov){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose(cov),
                          max_dose = max_dose(cov))
    beta <- rep(NA, 3)

    if (direction == "positive"){
      beta[1] <- logit(rho[1])
      beta[2] <- logit(rho[2]) - logit(rho[1])
      beta[3] <- logit(rho[3]) - logit(rho[1])
    } else {
      beta[1] <- logit(rho[1])
      beta[2] <- logit(rho[2]) - logit(rho[1])
      beta[3] <- -(logit(rho[3]) - logit(rho[1]))
    }

    cov <- (cov - min_cov)/(max_cov - min_cov)

    design_matrix <- cbind(1, dose, cov)
    lp <- design_matrix%*%beta
    out <- as.numeric(plogis(lp))
    return(out)
  }

  return(pdlt)
}

