#'@export
pdlt_d1basic <- function(rho, mtd, theta, min_dose, max_dose){

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose,
                         max_dose = max_dose)

  parm <- cbind(rho, gamma)
  parm <- as.matrix(parm)
  parm_names <- c("rho", "mtd")

  pdlt <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose,
                          max_dose = max_dose)

    pdlt_aux <- function(parm, theta, dose, parm_names) {
      rho <- parm[which(parm_names == "rho")]
      mtd <- parm[which(parm_names == "mtd")]

      beta <- rep(NA, 2)
      beta[1] <- logit(rho)
      beta[2] <- (logit(theta) - logit(rho))/mtd

      design_matrix <- cbind(1, dose)
      lp <- design_matrix %*% beta
      out <- as.numeric(plogis(lp))
      return(out)
    }

    out <- apply(parm, 1, pdlt_aux, theta = theta,
                 dose = dose, parm_names = parm_names)
    return(out)
  }

  return(pdlt)
}

#'@export
pdlt_d1extended <- function(rho, min_dose, max_dose){

  parm <- rho
  parm <- as.matrix(parm)
  parm_names <- "rho"

  pdlt <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose,
                          max_dose = max_dose)

    pdlt_aux <- function(parm, theta, dose, parm_names){

      rho <- parm
      beta <- rep(NA, 2)
      beta[1] <- logit(rho[1])
      beta[2] <- logit(rho[2]) - logit(rho[1])

      design_matrix <- cbind(1, dose)
      lp <- design_matrix%*%beta
      out <- as.numeric(plogis(lp))
      return(out)
    }

    out <- apply(parm, 1, pdlt_aux, theta = theta, dose = dose,
                 parm_names = parm_names)
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

  if (distribution == "exponential") {
    parm <- cbind(rho, gamma)
    parm_names <- c("rho", "mtd")
  } else {
    parm <- cbind(rho, gamma, shape)
    parm_names <- c("rho", "mtd", "shape")
  }
  parm <- as.matrix(parm)

  pdlt <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose,
                          max_dose = max_dose)

    pdlt_aux <- function(parm, dose, theta, tau, distribution, parm_names) {

      rho <- parm[which(parm_names == "rho")]
      gamma <- parm[which(parm_names == "mtd")]

      if (distribution == "weibull") {
        shape <- parm[which(parm_names == "shape")]
      } else {
        shape <- 1
      }

      beta <- rep(NA, 2)
      beta[1] <- log(-log(1 - rho[1])/(tau^shape))
      beta[2] <- log(log(1 - theta)/log(1 - rho[1]))/gamma

      design_matrix <- cbind(1, dose)

      out <- 1 - exp(-exp(design_matrix%*%beta)*(tau^shape))

      return(out)
    }

    out <- apply(parm, 1, pdlt_aux, dose = dose, theta = theta, tau = tau,
                 distribution = distribution, parm_names = parm_names)
    return(out)
  }

  return(pdlt)
}

#'@export
pdlt_d1multinomial <- function(rho, mtd, theta, min_dose, max_dose, levels_cov) {

  if (!is.matrix(mtd))
    mtd <- matrix(mtd, nrow = 1)

  gamma <- matrix(NA, ncol = ncol(mtd), nrow = nrow(mtd))

  for (i in 1:length(levels_cov)){
    gamma[, i] <- standard_dose(dose = mtd[, i],
                                min_dose = min_dose(levels_cov[i]),
                                max_dose = max_dose(levels_cov[i]))
  }

  parm <- cbind(gamma, rho)
  parm <- as.matrix(parm)
  parm.names <- c(rep("mtd", ncol(gamma)), "rho")

  pdlt <- function(dose, cov){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose(cov),
                          max_dose = max_dose(cov))

    aux_pdlt <- function(parm, dose, theta, cov, parm.names) {

      rho <- parm[which(parm.names == "rho")]
      mtd <- parm[which(parm.names == "mtd")]

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

    out <- apply(parm, 1, aux_pdlt, dose = dose, theta = theta,
                 cov = cov, parm.names = parm.names)

    return(out)
  }

  return(pdlt)
}


#'@export
pdlt_d1ordinal <- function(rho, mtd, theta, min_dose, max_dose, levels_cov) {

  if (!is.matrix(mtd))
    mtd <- matrix(mtd, nrow = 1)

  gamma <- matrix(NA, ncol = ncol(mtd), nrow = nrow(mtd))

  for (i in 1:length(levels_cov)){
    gamma[, i] <- standard_dose(dose = mtd[, i],
                                min_dose = min_dose(levels_cov[i]),
                                max_dose = max_dose(levels_cov[i]))
  }

  parm <- cbind(gamma, rho)
  parm <- as.matrix(parm)
  parm.names <- c(rep("mtd", ncol(gamma)), "rho")

  pdlt <- function(dose, cov){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose(cov),
                          max_dose = max_dose(cov))

    aux_pdlt <- function(parm, dose, theta, cov, parm.names) {

      rho <- parm[which(parm.names == "rho")]
      mtd <- parm[which(parm.names == "mtd")]

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

    out <- apply(parm, 1, aux_pdlt, dose = dose, theta = theta,
                 cov = cov, parm.names = parm.names)
    return(out)
  }

  return(pdlt)
}

#'@export
pdlt_d1continuous <- function(mtd, rho, theta, min_dose, max_dose,
                              min_cov, max_cov, cov) {

  gamma <- standard_dose(dose = mtd,
                         min_dose = min_dose(cov),
                         max_dose = max_dose(cov))

  parm <- cbind(gamma, rho)
  parm <- as.matrix(parm)
  parm.names <- c("mtd", "rho", "rho")

  dose <- standard_dose(dose = dose,
                        min_dose = min_dose(cov),
                        max_dose = max_dose(cov))

  pdlt <- function(dose){

    aux_pdlt <- function(parm, dose, mtd, theta, cov, min_cov, max_cov,
                         parm.names) {

      rho <- parm[which(parm.names == "rho")]
      mtd <- parm[which(parm.names == "mtd")]

      beta <- rep(NA, 3)
      beta[1] <- logit(rho[1]) - min_cov*(logit(rho[2]) - logit(rho[1]))/
        (max_cov - min_cov)
      beta[2] <- (logit(theta) - logit(rho[1]))/mtd
      beta[3] <- (logit(rho[2]) - logit(rho[1]))/(max_cov - min_cov)

      design_matrix <- cbind(1, dose, cov)
      lp <- design_matrix%*%beta
      out <- as.numeric(plogis(lp))
      return(out)
    }

    out <- apply(parm, 1, aux_pdlt, mtd = mtd, dose = dose, theta = theta,
                 cov = cov, min_cov = min_cov, max_cov = max_cov,
                 parm.names = parm.names)
    return(out)
  }
  return(pdlt)
}

#'@export
pdlt_d1excontinuous <- function(rho, theta, min_dose, max_dose,
                                min_cov, max_cov, cov) {

  parm <- cbind(rho)
  parm <- matrix(rho, ncol = 3)
  parm.names <- c("rho", "rho", "rho")

  pdlt <- function(dose){

    dose <- standard_dose(dose = dose,
                          min_dose = min_dose(cov),
                          max_dose = max_dose(cov))

    aux_pdlt <- function(parm, dose, mtd, theta, cov, min_cov, max_cov,
                         parm.names) {

      rho <- parm[which(parm.names == "rho")]

      beta <- rep(NA, 3)
      beta[1] <- logit(rho[1]) - min_cov*(logit(rho[2]) - logit(rho[1]))/(max_cov - min_cov)
      beta[2] <- logit(rho[3]) - logit(rho[1])
      beta[3] <- (logit(rho[2]) - logit(rho[1]))/(max_cov - min_cov)

      design_matrix <- cbind(1, dose, cov)
      lp <- design_matrix%*%beta
      out <- as.numeric(plogis(lp))
      return(out)
    }

    out <- apply(parm, 1, aux_pdlt, mtd = mtd, dose = dose, theta = theta,
                 cov = cov, min_cov = min_cov, max_cov = max_cov,
                 parm.names = parm.names)
    return(out)
  }
  return(pdlt)
}

