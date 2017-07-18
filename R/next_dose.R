#'@export
next_dose.d1basic <- function(data){

  rho <- data$mcmc$rho
  gamma <- data$mcmc$gamma
  mtd <-
    inv_standard_dose(dose = gamma,
                      min_dose = data$limits$min_dose,
                      max_dose = data$limits$max_dose)

  next_dose <- quantile(mtd, probs = data$alpha)

    out <- list(mtd = mtd, next_dose = next_dose,
                rho = rho, gamma = gamma,
                sample = data$mcmc$sample)
  return(out)
}


#'@export
next_dose.d1extended <- function(data){

  rho <- data$mcmc$rho

  scale <- logit(rho[, 2]) - logit(rho[, 1])
  gamma <- (logit(data$theta) - logit(rho[, 1]))/scale
  mtd <-
    inv_standard_dose(dose = gamma,
                      min_dose = data$limits$min_dose,
                      max_dose = data$limits$max_dose)

  next_dose <- quantile(mtd, probs = data$alpha)

  out <- list(mtd = mtd, next_dose = next_dose,
              rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}

#'@export
next_dose.d1ph <- function(data){

  gamma <- data$mcmc$gamma
  shape <- data$mcmc$shape
  rho <- data$mcmc$rho
  mtd <-
    inv_standard_dose(dose = gamma,
                      min_dose = data$limits$min_dose,
                      max_dose = data$limits$max_dose)

  next_dose <- quantile(mtd, probs = data$alpha)

  out <- list(mtd = mtd, next_dose = next_dose,
              rho = rho, shape = shape, gamma = gamma,
              sample = data$mcmc$sample)
  return(out)
}

#'@export
next_dose.d1multinomial <- function(data){

  rho <- data$mcmc$rho
  gamma <- data$mcmc$gamma

  mtd <- matrix(NA, nrow = nrow(gamma), ncol = ncol(gamma))
  next_dose <- rep(NA, ncol(gamma))

  for (i in 1:length(data$levels_cov)){
    mtd[, i] <-
      inv_standard_dose(dose = gamma[, i],
                        min_dose = data$limits$min_dose(data$levels_cov[i]),
                        max_dose = data$limits$max_dose(data$levels_cov[i]))
    next_dose[i] <- quantile(mtd[, i], probs = data$alpha)
  }

  out <- list(mtd = mtd, next_dose = next_dose,
              rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}

#'@export
next_dose.d1ordinal <- function(data){

  rho <- data$mcmc$rho
  gamma <- data$mcmc$gamma

  mtd <- matrix(NA, nrow = nrow(gamma), ncol = ncol(gamma))
  next_dose <- rep(NA, ncol(gamma))

  for (i in 1:length(data$levels_cov)){
    mtd[, i] <-
      inv_standard_dose(dose = gamma[, i],
                        min_dose = data$limits$min_dose(data$levels_cov[i]),
                        max_dose = data$limits$max_dose(data$levels_cov[i]))
    next_dose[i] <- quantile(mtd[, i], probs = data$alpha)
  }

  out <- list(mtd = mtd, next_dose = next_dose,
              rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}

#'@export
next_dose.d1continuous <- function(data){

  rho <- data$mcmc$rho
  gamma <- data$mcmc$gamma

  mtd <-
    inv_standard_dose(dose = gamma,
                      min_dose = data$limits$min_dose(data$next_patient_cov),
                      max_dose = data$limits$max_dose(data$next_patient_cov))

  next_dose <- quantile(mtd, probs = data$alpha)

  out <- list(mtd = mtd, next_dose = next_dose,
              rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}


#'@export
next_dose.d1excontinuous <- function(data){

  rho <- data$mcmc$rho

  scale_p0 <- (data$next_patient_cov - data$min_cov)/(data$max_cov - data$min_cov)
  scale_p1 <- logit(rho[, 2]) - logit(rho[, 1])
  scale_p2 <- logit(data$theta) - logit(rho[, 1])
  scale_p3 <- logit(rho[, 3]) - logit(rho[, 1])
  gamma <- (scale_p2 - scale_p0*scale_p1)/scale_p3
  mtd <-
    inv_standard_dose(dose = gamma,
                      min_dose = data$limits$min_dose(data$next_patient_cov),
                      max_dose = data$limits$max_dose(data$next_patient_cov))

  next_dose <- quantile(mtd, probs = data$alpha)

  out <- list(mtd = mtd, next_dose = next_dose,
              rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}
