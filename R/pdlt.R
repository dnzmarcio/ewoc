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

#'@export
pdlt_d1multinomial <- function(rho, mtd, theta, min_dose, max_dose, levels_cov) {

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose(levels_cov),
                         max_dose = max_dose(levels_cov))

  pdlt <- function(dose, cov){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose(cov),
                          max_dose = max_dose(cov))

    beta <- rep(NA, (length(mtd) + 1))
    beta[1] <- logit(rho[1])
    beta[2] <- (logit(theta) - logit(rho[1]))/mtd[1]

    for (i in 3:(length(mtd) + 1)) {
      beta[i] <- logit(theta) - logit(rho[1]) -
        mtd[i-1]*(logit(theta) - beta[1])/mtd[i-2]
    }

    cov <- factor(cov, levels = levels_cov)
    cov <- matrix(model.matrix( ~ cov)[-1], nrow = 1)

    design_matrix <- cbind(1, dose, cov)
    lp <- design_matrix%*%beta

    out <- as.numeric(plogis(lp))
    return(out)
  }

  return(pdlt)
}


#'@export
pdlt_d1continuous <- function(mtd, rho, theta, direction,
                              min_dose, max_dose,
                              min_cov, max_cov) {

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose(cov),
                         max_dose = max_dose(cov))

  pdlt <- function(dose, cov){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose(cov),
                          max_dose = max_dose(cov))

    beta <- rep(NA, 3)
    beta[1] <- logit(rho[1])
    beta[2] <- (logit(theta) - logit(rho[1]))/gamma
    beta[3] <- (logit(rho[2]) - logit(rho[1]))

    cov <- (cov - min_cov)/(max_cov - min_cov)

    design_matrix <- cbind(1, dose, cov)
    lp <- design_matrix%*%beta
    out <- as.numeric(plogis(lp))
    return(out)

  }
  return(pdlt)
}

#'@export
pdlt_d1excontinuous <- function(rho, direction,
                                min_dose, max_dose, min_cov, max_cov) {

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose(max_cov),
                         max_dose = max_dose(max_cov))

  pdlt <- function(dose, cov){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose(cov),
                          max_dose = max_dose(cov))

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

