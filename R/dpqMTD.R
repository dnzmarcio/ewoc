#'@export
dmtd <- function(mtd, data) {

  #Posterior distribution of MTD without the constant of integration
  aux_dmtd <- function (mtd, data) {

    #Bivariate posterior distribution of rho and MTD without the constant
    #of integration
    drho_mtd <- function (rho, mtd, data) {

      #Parameters
      theta <- data$theta
      min_dose <- 0
      max_dose <- 1
      rho_prior <- data$rho_prior
      mtd_prior <- data$mtd_prior
      design_matrix <- data$design_matrix
      DLT <- data$response[, 1]
      npatients <- data$response[, 2]

      if (rho <= 0 | rho >= theta |
          mtd <= 0 | mtd >= 1) {
        out <- 0
      } else {

        #Log-prior
        rho_prior <- theta*dbeta(rho, rho_prior[1], rho_prior[2], log = TRUE)
        mtd_prior <- dbeta(mtd, mtd_prior[2], mtd_prior[2], log = TRUE)

        #Log-likelihood
        beta_0 <- qlogis(rho)
        beta_1 <- (qlogis(theta) - qlogis(rho))/mtd
        beta <- c(beta_0, beta_1)
        eta <- tcrossprod(beta, design_matrix)
        p <- plogis(eta)
        LL <- sum(dbinom(DLT, npatients, p, log = TRUE))

        #Log-Posterior
        LP <- LL + rho_prior + mtd_prior

        out <- exp(LP)
      }

      return(out)
    }

    theta <- data$theta
    auxVec_drho_mtd <- Vectorize(drho_mtd, vectorize.args = "rho")

    out <- distrExIntegrate(auxVec_drho_mtd, lower = 0, upper = theta,
                     mtd = mtd, data = data)

    return(out)
  }

  auxVec_dmtd <- Vectorize(aux_dmtd, vectorize.args = "mtd")

  #Constant of integration of the bivariate distribution
  C <- distrExIntegrate(auxVec_dmtd,
                 lower = 0, upper = 1,
                 data = data)

  out <- auxVec_dmtd(mtd, data)/C

  return(out)
}

#'@export
pmtd <- function(mtd, data) {

  out <- ifelse(mtd > 0 & mtd < 1,
                distrExIntegrate(dmtd, lower = 0, upper = mtd,
                          data = data),
                ifelse(mtd <= 0, 0, 1))

  return(out)
}

#'@export
qmtd <- function(data) {

  aux_qmtd <- function(mtd, alpha, data) {

    out <- pmtd(mtd, data) - alpha
    return(out)
  }

  gamma <- uniroot(aux_qmtd, lower = 0,
                       upper = 1,
                       alpha = data$alpha, data = data)$root

  next_dose <-
    inv_standard_dose(dose = gamma,
                      min_dose = data$limits$min_dose(),
                      max_dose = data$limits$max_dose())
  out <-list(next_dose = next_dose)

  return(out)
}

#'@export
qmtd_jags <- function(data, n_adapt, burn_in, n_mcmc, n_thin, n_chains) {

  data$mcmc <- ewoc_jags(data, n_adapt, burn_in, n_mcmc, n_thin, n_chains)
  out <- next_dose(data)

  if (data$type == "discrete")
    out$next_dose <- rounding_system(dose = out$next_dose,
                                     grid = data$dose_set,
                                     rounding = data$rounding)

  return(out)
}


