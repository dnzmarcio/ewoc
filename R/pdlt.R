#'@export
pdlt_d1basic <- function(dose, rho, gamma, theta, min_dose, max_dose){

  parm <- cbind(rho, gamma)
  parm_names <- c("rho", "mtd")

  dose <- standard_dose(dose = dose,
                        min_dose = min_dose,
                        max_dose = max_dose)

  if (!is.matrix(parm))
    parm <- as.matrix(parm)

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

#'@export
pdlt_d1extended <- function(dose, rho, min_dose, max_dose){

  parm <- rho
  parm_names <- "rho"

  dose <- standard_dose(dose = dose,
                        min_dose = min_dose,
                        max_dose = max_dose)

  if (!is.matrix(parm))
    parm <- as.matrix(parm)

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

#'@export
pdlt_d1ph <- function(dose, rho, gamma, shape = NULL, theta, min_dose, max_dose,
                      tau, distribution) {

  if (distribution == "exponential") {
    parm <- cbind(rho, gamma)
    parm_names <- c("rho", "mtd")
  } else {
    parm <- cbind(rho, gamma, shape)
    parm_names <- c("rho", "mtd", "shape")
  }

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

#'@export
pdlt_d1multinomial <- function(dose, rho, gamma, theta, min_dose, max_dose, cov) {

  dose <- standard_dose(dose = dose,
                        min_dose = min_dose,
                        max_dose = max_dose)

  parm <- cbind(gamma, rho)
  parm.names <- c(rep("mtd", ncol(gamma)), "rho")

  aux_pdlt <- function(parm, dose, theta, cov, parm.names) {

    rho <- parm[which(parm.names == "rho")]
    mtd <- parm[which(parm.names == "mtd")]

    beta <- rep(NA, (length(mtd) + 1))
    beta[1] <- logit(rho[1])
    beta[2] <- (logit(theta) - logit(rho[1]))/mtd[1]

    for (i in 3:(length(mtd) + 1)) {
      beta[i] <- logit(theta) - logit(rho[1]) - mtd[i-1]*(logit(theta) - beta[1])/mtd[i-2]
    }

    design_matrix <- cbind(1, dose, matrix(cov[-1], nrow = 1))
    lp <- design_matrix%*%beta
    out <- as.numeric(plogis(lp))
    return(out)
  }

  out <- apply(parm, 1, aux_pdlt, dose = dose, theta = theta,
               cov = cov, parm.names = parm.names)
  return(out)
}

#'@export
pdlt_d1ordinal <- function(dose, gamma, rho, theta,  min_dose, max_dose, cov) {

  dose <- standard_dose(dose = dose,
                        min_dose = min_dose,
                        max_dose = max_dose)

  parm <- cbind(gamma, rho)
  parm.names <- c(rep("mtd", ncol(gamma)), "rho")

  aux_pdlt <- function(parm, dose, mtd, theta, cov, parm.names) {

    rho <- parm[which(parm.names == "rho")]
    mtd <- parm[which(parm.names == "mtd")]

    beta <- rep(NA, (length(mtd) + 1))
    beta[1] <- logit(rho[1])
    beta[2] <- (logit(theta) - logit(rho[1]))/mtd[1]

    for (i in 3:(length(mtd) + 1)) {
      beta[i] <- logit(theta) - logit(rho[1]) - mtd[i-1]*(logit(theta) - beta[1])/mtd[i-2]
    }

    design_matrix <- cbind(1, dose, matrix(cov[-1], nrow = 1))
    lp <- design_matrix%*%beta
    out <- as.numeric(plogis(lp))
    return(out)
  }

  out <- apply(parm, 1, aux_pdlt, mtd = mtd, dose = dose, theta = theta,
               cov = cov, parm.names = parm.names)
  return(out)
}

