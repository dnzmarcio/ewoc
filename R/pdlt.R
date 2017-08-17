#'@export
pdlt_d1basic <- function(rho, mtd, theta, min_dose, max_dose){

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

