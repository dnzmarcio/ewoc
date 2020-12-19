next_dose.ewoc_d1classical <- function(data){

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

  if (data$type == "continuous")
    if ((next_dose - data$current_dose) > data$max_increment)
      next_dose <- data$current_dose + data$max_increment

  if (data$type == "discrete"){
    next_dose <- rounding_system(dose = next_dose,
                                 grid = data$dose_set,
                                 rounding = data$rounding)

    if (data$no_skip_dose)
      if (which(data$dose_set == next_dose) -
          which(data$dose_set == data$current_dose) > 1)
        next_dose <- data$dose_set[(which(data$dose_set == data$current_dose) + 1)]
  }

  next_gamma <- standard_dose(dose = next_dose,
                              min_dose = data$limits$min_dose,
                              max_dose = data$limits$max_dose)

  pdlt <- as.numeric(plogis(cbind(1, next_gamma)%*%t(beta)))

    out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
                rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}


next_dose.ewoc_d1extended <- function(data){

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

  if (data$type == "continuous")
    if ((next_dose - data$current_dose) > data$max_increment)
      next_dose <- data$current_dose + data$max_increment

  if (data$type == "discrete"){
    next_dose <- rounding_system(dose = next_dose,
                                 grid = data$dose_set,
                                 rounding = data$rounding)

    if (data$no_skip_dose)
      if (which(data$dose_set == next_dose) -
          which(data$dose_set == data$current_dose) > 1)
        next_dose <- data$dose_set[(which(data$dose_set == data$current_dose) + 1)]
  }

  next_gamma <- standard_dose(dose = next_dose,
                              min_dose = data$limits$min_dose,
                              max_dose = data$limits$max_dose)

  pdlt <- as.numeric(plogis(cbind(1, next_gamma)%*%t(beta)))

  out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
              rho = rho, gamma = gamma, sample = data$mcmc$sample)
  return(out)
}

next_dose.ewoc_d1ph <- function(data){

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

  if (data$type == "continuous")
    if ((next_dose - data$current_dose) > data$max_increment)
      next_dose <- data$current_dose + data$max_increment

  if (data$type == "discrete"){
    next_dose <- rounding_system(dose = next_dose,
                                 grid = data$dose_set,
                                 rounding = data$rounding)

    if (data$no_skip_dose)
      if (which(data$dose_set == next_dose) -
          which(data$dose_set == data$current_dose) > 1)
        next_dose <- data$dose_set[(which(data$dose_set == data$current_dose) + 1)]
  }

  next_gamma <- standard_dose(dose = next_dose,
                              min_dose = data$limits$min_dose,
                              max_dose = data$limits$max_dose)

  if (data$distribution != "weibull")
    shape <- 1

  pdlt <- as.numeric(1 - exp(-exp(cbind(1, next_gamma)%*%t(beta))*
                             (data$tau^shape)))

  out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
              rho = rho, shape = shape, gamma = gamma,
              sample = data$mcmc$sample)
  return(out)
}

#'@export
next_dose.d1dicov <- function(data){

  rho <- data$mcmc$rho
  gamma <- data$mcmc$gamma
  beta <- data$mcmc$beta

  mtd <- matrix(NA, nrow = nrow(gamma), ncol = length(data$levels_cov))
  next_dose <- rep(NA, length(data$levels_cov))
  aux_gamma <- matrix(NA, nrow = nrow(gamma), ncol = length(data$levels_cov))
  next_gamma <- rep(NA, length(data$levels_cov))


  for (i in 1:length(data$levels_cov)){
    index <- which(data$levels_cov[i] == data$levels_cov)

    next_gamma[i] <- quantile(gamma[, index], probs = data$alpha)

    mtd[, i] <- inv_standard_dose(dose = gamma[, index],
                        min_dose =
                          data$limits$min_dose(data$levels_cov[i]),
                        max_dose =
                          data$limits$max_dose(data$levels_cov[i]))

    aux_gamma[, i] <- gamma[, index]

    next_dose[i] <- inv_standard_dose(dose = next_gamma[i],
                                      min_dose =
                                        data$limits$min_dose(data$levels_cov[i]),
                                      max_dose =
                                        data$limits$max_dose(data$levels_cov[i]))
    next_dose[i] <- ifelse(next_dose[i] >
                           data$limits$last_dose(data$levels_cov[index]),
                           data$limits$last_dose(data$levels_cov[index]),
                           ifelse(next_dose[i] <
                                  data$limits$first_dose(data$levels_cov[index]),
                                  data$limits$first_dose(data$levels_cov[index]),
                                  next_dose[i]))
  }

  cov <- factor(data$levels_cov, levels = data$levels_cov)
  cov <- matrix(model.matrix(~ cov)[, -1])

  temp <- cbind(1, next_gamma, cov)
  pdlt <- as.matrix(plogis(temp%*%t(beta)))


  out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
              rho = rho, gamma = aux_gamma, sample = data$mcmc$sample)
  return(out)
}

