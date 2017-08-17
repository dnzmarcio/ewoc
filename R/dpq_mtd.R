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


