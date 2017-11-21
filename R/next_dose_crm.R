next_dose.crm_d1classic <- function(data){

  beta <- data$mcmc$beta
  beta_est <- apply(beta, 2, median)

  mtd <- inv_standard_dose(dose = data$mcmc$gamma,
                           min_dose = data$limits$min_dose,
                           max_dose = data$limits$max_dose)

  crm_criterion <- function(dose){
    design_matrix <- cbind(1, dose)
    lp <- design_matrix %*% beta_est
    out <- abs(as.numeric(plogis(lp)) - data$theta)
    return(out)
  }

  next_gamma <- optimize(f = crm_criterion,
                        lower = 0,
                        upper = 1)$minimum
  next_dose <- inv_standard_dose(dose = next_gamma,
                                 min_dose = data$limits$min_dose,
                                 max_dose = data$limits$max_dose)
  next_dose <- ifelse(next_dose > data$limits$last_dose,
                      data$limits$last_dose,
                      ifelse(next_dose < data$limits$first_dose,
                             data$limits$first_dose, next_dose))

  if (abs(next_dose - data$current_dose) > data$max_increment)
    next_dose <- data$current_dose + data$max_increment

  if (data$type == "discrete")
    next_dose <- rounding_system(dose = next_dose,
                                 grid = data$dose_set,
                                 rounding = data$rounding)

  next_gamma <- standard_dose(dose = next_dose,
                              min_dose = data$limits$min_dose,
                              max_dose = data$limits$max_dose)

  pdlt <- as.numeric(plogis(cbind(1, next_gamma)%*%t(beta)))

    out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
                rho = data$mcmc$rho, gamma = data$mcmc$gamma,
                sample = data$mcmc$sample)
  return(out)
}


next_dose.crm_d1extended <- function(data){

  beta <- data$mcmc$beta
  beta_est <- apply(beta, 2, median)

  scale <- logit(data$mcmc$rho[, 2]) - logit(data$mcmc$rho[, 1])
  gamma <- (logit(data$theta) - logit(data$mcmc$rho[, 1]))/scale

  mtd <- inv_standard_dose(dose = gamma,
                           min_dose = data$limits$min_dose,
                           max_dose = data$limits$max_dose)

  crm_criterion <- function(dose){
    design_matrix <- cbind(1, dose)
    lp <- design_matrix %*% beta_est
    out <- abs(as.numeric(plogis(lp)) - data$theta)
    return(out)
  }

  next_gamma <- optimize(f = crm_criterion,
                         lower = 0,
                         upper = 1)$minimum
  next_dose <- inv_standard_dose(dose = next_gamma,
                                 min_dose = data$limits$min_dose,
                                 max_dose = data$limits$max_dose)

  next_dose <- ifelse(next_dose > data$limits$last_dose,
                      data$limits$last_dose,
                      ifelse(next_dose < data$limits$first_dose,
                             data$limits$first_dose, next_dose))

  if (abs(next_dose - data$current_dose) > data$max_increment)
    next_dose <- data$current_dose + data$max_increment

  if (data$type == "discrete")
    next_dose <- rounding_system(dose = next_dose,
                                 grid = data$dose_set,
                                 rounding = data$rounding)

  next_gamma <- standard_dose(dose = next_dose,
                              min_dose = data$limits$min_dose,
                              max_dose = data$limits$max_dose)

  pdlt <- as.numeric(plogis(cbind(1, next_gamma)%*%t(beta)))

  out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
              rho = data$mcmc$rho, gamma = data$mcmc$gamma,
              sample = data$mcmc$sample)
  return(out)
}

next_dose.crm_d1ph <- function(data){

  shape <- data$mcmc$shape
  beta <- data$mcmc$beta
  beta_est <- apply(beta, 2, median)

  mtd <- inv_standard_dose(dose = data$mcmc$gamma,
                           min_dose = data$limits$min_dose,
                           max_dose = data$limits$max_dose)

  crm_criterion <- function(dose){
    design_matrix <- cbind(1, dose)
    lp <- design_matrix %*% beta_est

    if (data$distribution != "weibull")
      shape <- 1

    out <- abs(as.numeric(1 - exp(-exp(lp)*(data$tau^shape))) - data$theta)
    return(out)
  }

  next_gamma <- optimize(f = crm_criterion,
                         lower = 0,
                         upper = 1)$minimum
  next_dose <- inv_standard_dose(dose = next_gamma,
                                 min_dose = data$limits$min_dose,
                                 max_dose = data$limits$max_dose)

  next_dose <- ifelse(next_dose > data$limits$last_dose,
                      data$limits$last_dose,
                      ifelse(next_dose < data$limits$first_dose,
                             data$limits$first_dose, next_dose))

  if (abs(next_dose - data$current_dose) > data$max_increment)
    next_dose <- data$current_dose + data$max_increment

  if (data$type == "discrete")
    next_dose <- rounding_system(dose = next_dose,
                                 grid = data$dose_set,
                                 rounding = data$rounding)

  next_gamma <- standard_dose(dose = next_dose,
                              min_dose = data$limits$min_dose,
                              max_dose = data$limits$max_dose)

  if (data$distribution != "weibull")
    shape <- 1
  pdlt <- as.numeric(1 - exp(-exp(cbind(1, next_gamma)%*%t(beta))*
                               (data$tau^shape)))

  out <- list(mtd = mtd, pdlt = pdlt, next_dose = next_dose,
              rho = data$mcmc$rho, shape = data$mcmc$shape, gamma = data$mcmc$gamma,
              sample = data$mcmc$sample)
  return(out)
}


