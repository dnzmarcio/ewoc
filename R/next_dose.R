#'@export
next_dose.d1basic <- function(data){

  rho <- data$mcmc$rho
  gamma <- data$mcmc$gamma
  beta <- data$mcmc$beta

  mtd <- inv_standard_dose(dose = gamma,
                           min_dose = data$limits$min_dose,
                           max_dose = data$limits$max_dose)

  next_dose <- quantile(mtd, probs = data$alpha)

  next_dose <- ifelse(next_dose > data$limits$last_dose,
                      data$limits$last_dose,
                      ifelse(next_dose < data$limits$first_dose,
                             data$limits$first_dose, next_dose))

  next_gamma <- standard_dose(dose = next_dose,
                              min_dose = data$limits$min_dose,
                              max_dose = data$limits$max_dose)

  pdlt <- as.numeric(plogis(cbind(1, next_gamma)%*%t(beta)))

    out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
                rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}


#'@export
next_dose.d1extended <- function(data){

  rho <- data$mcmc$rho
  beta <- data$mcmc$beta

  scale <- logit(rho[, 2]) - logit(rho[, 1])
  gamma <- (logit(data$theta) - logit(rho[, 1]))/scale
  mtd <- inv_standard_dose(dose = gamma,
                           min_dose = data$limits$min_dose,
                           max_dose = data$limits$max_dose)

  next_dose <- quantile(mtd, probs = data$alpha)

  next_dose <- ifelse(next_dose > data$limits$last_dose,
                      data$limits$last_dose,
                      ifelse(next_dose < data$limits$first_dose,
                             data$limits$first_dose, next_dose))

  next_gamma <- standard_dose(dose = next_dose,
                              min_dose = data$limits$min_dose,
                              max_dose = data$limits$max_dose)

  pdlt <- as.numeric(plogis(cbind(1, next_gamma)%*%t(beta)))

  out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
              rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}

#'@export
next_dose.d1ph <- function(data){

  gamma <- data$mcmc$gamma - 10^(-2)
  shape <- data$mcmc$shape
  rho <- data$mcmc$rho
  beta <- data$mcmc$beta

  mtd <- inv_standard_dose(dose = gamma,
                           min_dose = data$limits$min_dose,
                           max_dose = data$limits$max_dose)

  next_dose <- quantile(mtd, probs = data$alpha)

  next_dose <- ifelse(next_dose > data$limits$last_dose,
                      data$limits$last_dose,
                      ifelse(next_dose < data$limits$first_dose,
                             data$limits$first_dose, next_dose))

  next_gamma <- standard_dose(dose = next_dose,
                              min_dose = data$limits$min_dose,
                              max_dose = data$limits$max_dose)

  if (data$distribution != "weibull")
    shape <- 1
  pdlt <- as.numeric(1 -
                       exp(-exp(cbind(1, next_gamma)%*%t(beta))*
                             (data$tau^shape)))

  out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
              rho = rho, shape = shape, gamma = gamma,
              sample = data$mcmc$sample)
  return(out)
}

#'@export
next_dose.d1multinomial <- function(data){

  rho <- data$mcmc$rho
  gamma <- data$mcmc$gamma
  beta <- data$mcmc$beta

  mtd <- matrix(NA, nrow = nrow(gamma), ncol = length(data$next_patient_cov))
  next_dose <- rep(NA, length(data$next_patient_cov))
  next_gamma <- rep(NA, length(data$next_patient_cov))

  for (i in 1:length(data$next_patient_cov)){
    index <- which(data$next_patient_cov[i] == data$levels_cov)

    mtd[, i] <- inv_standard_dose(dose = gamma[, i],
                        min_dose =
                          data$limits$min_dose(data$next_patient_cov[i]),
                        max_dose =
                          data$limits$max_dose(data$next_patient_cov[i]))

    next_dose[i] <- quantile(mtd[, i], probs = data$alpha)

    next_dose[i] <- ifelse(next_dose[i] >
                           data$limits$last_dose(data$levels_cov[index]),
                           data$limits$last_dose(data$levels_cov[index]),
                           ifelse(next_dose[i] <
                                  data$limits$first_dose(data$levels_cov[index]),
                                  data$limits$first_dose(data$levels_cov[index]),
                                  next_dose[i]))

    next_gamma[i] <- standard_dose(dose = next_dose[i],
                        min_dose =
                          data$limits$min_dose(data$next_patient_cov[i]),
                        max_dose =
                          data$limits$max_dose(data$next_patient_cov[i]))
  }

  cov <- factor(data$next_patient_cov, levels = data$levels_cov)
  cov <- matrix(model.matrix(~ cov)[-1], nrow = 1)

  temp <- cbind(1, next_gamma, cov)
  pdlt <- as.numeric(plogis(temp%*%t(beta)))


  out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
              rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}

#'@export
next_dose.d1continuous <- function(data){

  rho <- data$mcmc$rho
  gamma <- data$mcmc$gamma

  mtd <- inv_standard_dose(dose = gamma,
                           min_dose =
                             data$limits$min_dose(data$next_patient_cov),
                           max_dose =
                             data$limits$max_dose(data$next_patient_cov))

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
  mtd <- inv_standard_dose(dose = gamma,
                           min_dose =
                             data$limits$min_dose(data$next_patient_cov),
                           max_dose =
                             data$limits$max_dose(data$next_patient_cov))

  next_dose <- quantile(mtd, probs = data$alpha)

  out <- list(mtd = mtd, next_dose = next_dose,
              rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}
