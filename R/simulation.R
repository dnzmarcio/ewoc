#'Simulation of a trial using Escalation Over with Dose Control
#'
#'Performing a simulation of several phase I clinical trial based on the
#'Escalation Over Dose Control (EWOC).
#'
#'@param step_zero an object from the classes 'ewoc_d1basic', 'ewoc_d1extended',
#''ewoc_d1ph', 'ewoc_d1multinomial', 'ewoc_d1ordinal', and 'ewoc_d1continuous'
#'created using dummy data.
#'@param n_sim a number indicating the number of phase I clinical trials
#'to be simulated.
#'@param sample_size a number indicating the number of patients enrolled for
#'each clinical trial.
#'@param alpha_strategy a character indicating the strategy to apply for the
#'feasibility value. Default is "constant". Options are "increasing" and
#'"conditional".
#'@param alpha_rate a numerical value indicating the rate of the
#'feasibility strategy. Only necessary if alpha_strategy is either
#''increasing' or 'conditional'.
#'@param response_sim a function which is self-contained and will be used
#'as a generator function of the response variables in the simulation.
#'Its only imput is 'dose' and output is the DLT time.
#'
#'@examples
#'\dontrun{
### Basic EWOC
#'DLT <- 0
#'dose <- 30
#'step_zero <- ewoc_d1basic(DLT ~ dose, type = 'discrete',
#'                          theta = 0.33, alpha = 0.25,
#'                          min_dose = 0, max_dose = 100,
#'                          dose_set = seq(0, 100, 20),
#'                          rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                          mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                          rounding = "nearest")
#'response_sim <- response_d1basic(rho = 0.05, mtd = 20, theta = 0.33,
#'                                 min_dose = 10, max_dose = 50)
#'sim <- trial_simulation(step_zero = step_zero,
#'                        n_sim = 1, sample_size = 30,
#'                        alpha_strategy = "constant",
#'                        response_sim = response_sim)
#'
### Extended EWOC
#'DLT <- 0
#'dose <- 30
#'step_zero <- ewoc_d1extended(DLT ~ dose, type = 'discrete',
#'                            theta = 0.33, alpha = 0.25,
#'                            min_dose = 0, max_dose = 100,
#'                            dose_set = seq(0, 100, 20),
#'                            rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                            mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                            rounding = "nearest")
#'response_sim <- response_d1extended(rho = c(0.05, 0.5), theta = 0.33,
#'                                    min_dose = 10, max_dose = 50)
#'sim <- trial_simulation(step_zero = step_zero,
#'                        n_sim = 1, sample_size = 30,
#'                        alpha_strategy = "constant",
#'                        response_sim = response_sim)
#'
#'### PH EWOC
#'time <- 0
#'status <- 0
#'dose <- 30
#'
#'step_zero <- ewoc_d1ph(cbind(time, status) ~ dose, type = 'discrete',
#'                       theta = 0.33, alpha = 0.25, tau = 10,
#'                       min_dose = 30, max_dose = 50,
#'                       dose_set = seq(30, 50, 5),
#'                       rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                       mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                       distribution = 'exponential',
#'                       rounding = 'nearest')
#'response_sim <- response_d1ph(rho = 0.05, mtd = 40, theta = 0.33,
#'                                 min_dose = 30, max_dose = 50,
#'                                 tau = 10, distribution = "exponential")
#'sim <- trial_simulation(step_zero = step_zero,
#'                        n_sim = 1, sample_size = 30,
#'                        alpha_strategy = "increasing",
#'                        response_sim = response_sim)
#'
#'DLT <- rep(0, 1)
#'group <- "B"
#'dose <- rep(30, 1)
#'step_zero <- ewoc_d1multinomial(DLT ~ dose | group,
#'                                type = 'continuous',
#'                                theta = 0.33, alpha = 0.25,
#'                                min_dose = 30, max_dose = 50,
#'                                levels_cov = c("A", "B", "C"),
#'                                next_patient_cov = "A",
#'                                mtd_prior = matrix(1, nrow = 3, ncol = 2),
#'                                rho_prior = matrix(1, nrow = 1, ncol = 2))
#'
#'}
#'
#'
#'@export
trial_simulation <- function(step_zero, n_sim, sample_size,
                             alpha_strategy = "fixed",
                             alpha_rate = NULL, response_sim,
                             covariable_sim = NULL,
                             stop_rule_sim = NULL){
  UseMethod("trial_simulation")
}

#'@export
trial_simulation.d1basic <- function(step_zero, n_sim, sample_size,
                                     alpha_strategy =
                                       c("fixed", "increasing", "conditional"),
                                     alpha_rate = 0.05,
                                     response_sim = NULL,
                                     stop_rule_sim = NULL){

  if (is.null(response_sim))
    stop("'response_sim' function should be defined.")

  n_dose <- sample_size + 1
  dose_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  dlt_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  mtd_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  rho_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  alpha_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)


  for (i in 1:n_sim){

    dlt <- as.numeric(step_zero$trial$response)
    dose <- as.numeric(step_zero$trial$design_matrix[, 2])
    alpha_sim[, 1:length(dose)] <- as.numeric(step_zero$trial$alpha)

    for (j in (length(dose)+1):n_dose) {

      formula <- dlt[1:(j-1)] ~ dose[1:(j-1)]
      resolution <- ifelse(!is.na(dlt), 1, 0)

      if (j <= sample_size){
        alpha_sim[i, j] <- feasibility(alpha = alpha_sim[i, 1:(j-1)],
                                       strategy = alpha_strategy,
                                       dlt = dlt[1:(j-1)],
                                       resolution = resolution[1:(j-1)],
                                       rate = alpha_rate)

        update <- ewoc_d1basic(formula,
                               type = step_zero$trial$type,
                               theta = step_zero$trial$theta,
                               alpha = alpha_sim[i, j],
                               min_dose = step_zero$trial$min_dose,
                               max_dose = step_zero$trial$max_dose,
                               first_dose = step_zero$trial$first_dose,
                               last_dose = step_zero$trial$last_dose,
                               dose_set = step_zero$trial$dose_set,
                               rho_prior = step_zero$trial$rho_prior,
                               mtd_prior = step_zero$trial$mtd_prior,
                               rounding = step_zero$trial$rounding)

        if (!is.null(stop_rule_sim))
          if (stop_rule_sim(update)){
            dose[j:sample_size] <- NA
            dlt[j:sample_size] <- NA
            mtd_estimate <- NA
            rho_estimate <- NA
            break
          }

        dose[j] <- update$next_dose
        dlt[j] <- response_sim(dose = dose[j])
        mtd_estimate <- update$next_dose
        rho_estimate <- median(update$rho)
      }
    }

    dose_sim[i, ] <- dose
    dlt_sim[i, ] <- dlt
    mtd_sim[i] <- mtd_estimate
    rho_sim[i] <- rho_estimate
  }

  out <- list(dose_sim = dose_sim, dlt_sim = dlt_sim,
              mtd_sim = mtd_sim, rho_sim = rho_sim, alpha_sim = alpha_sim)
  return(out)
}

#'@export
trial_simulation.d1extended <- function(step_zero, n_sim, sample_size,
                                        alpha_strategy =
                                          c("fixed", "increasing", "conditional"),
                                        alpha_rate = 0.05,
                                        response_sim = NULL,
                                        stop_rule_sim = NULL){

  if (is.null(response_sim))
    stop("'response_sim' function should be defined.")

  n_dose <- sample_size + 1
  dose_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  dlt_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  mtd_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  rho_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  alpha_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)


  for (i in 1:n_sim){

    dlt <- as.numeric(step_zero$trial$response)
    dose <- as.numeric(step_zero$trial$design_matrix[, 2])
    alpha_sim[, 1:length(dose)] <- as.numeric(step_zero$trial$alpha)

    for (j in (length(dose)+1):n_dose) {

      formula <- dlt[1:(j-1)] ~ dose[1:(j-1)]
      resolution <- ifelse(!is.na(dlt), 1, 0)

      if (j <= sample_size){
        alpha_sim[i, j] <- feasibility(alpha = alpha_sim[i, 1:(j-1)],
                                       strategy = alpha_strategy,
                                       dlt = dlt[1:(j-1)],
                                       resolution = resolution[1:(j-1)],
                                       rate = alpha_rate)

        update <- ewoc_d1extended(formula,
                                  type = step_zero$trial$type,
                                  theta = step_zero$trial$theta,
                                  alpha = alpha_sim[i, j],
                                  min_dose = step_zero$trial$min_dose,
                                  max_dose = step_zero$trial$max_dose,
                                  first_dose = step_zero$trial$first_dose,
                                  last_dose = step_zero$trial$last_dose,
                                  dose_set = step_zero$trial$dose_set,
                                  rho_prior = step_zero$trial$rho_prior,
                                  mtd_prior = step_zero$trial$mtd_prior,
                                  rounding = step_zero$trial$rounding)

        if (!is.null(stop_rule_sim))
          if (stop_rule_sim(update)){
            dose[j:sample_size] <- NA
            dlt[j:sample_size] <- NA
            mtd_estimate <- NA
            rho_estimate <- NA
            break
          }

        dose[j] <- update$next_dose
        dlt[j] <- response_sim(dose = dose[j])
        mtd_estimate <- update$next_dose
        rho_estimate <- median(update$rho)
      }
    }

    dose_sim[i, ] <- dose
    dlt_sim[i, ] <- dlt
    mtd_sim[i] <- mtd_estimate
    rho_sim[i] <- rho_estimate
  }

  out <- list(dose_sim = dose_sim, dlt_sim = dlt_sim,
              mtd_sim = mtd_sim, rho_sim = rho_sim, alpha_sim = alpha_sim)
  return(out)
}


#'@export
trial_simulation.d1ph <- function(step_zero, n_sim, sample_size,
                                  alpha_strategy =
                                    c("fixed", "increasing", "conditional"),
                                  alpha_rate = 0.05,
                                  response_sim = NULL,
                                  stop_rule_sim = NULL){

  if (is.null(response_sim))
    stop("'response_sim' function should be defined.")
  if (length(alpha_strategy) != 1)
    stop("'alpha_strategy' should be defined.")

  n_dose <- sample_size + 1
  dose_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  dlt_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  time_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  mtd_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  rho_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  alpha_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)

  for (i in 1:n_sim){

    dlt <- as.numeric(step_zero$trial$response[, 2])
    dose <- as.numeric(step_zero$trial$design_matrix[, 2])
    alpha_sim[, 1:length(dose)] <- as.numeric(step_zero$trial$alpha)

    event_time <- ifelse(dlt == 1, as.numeric(step_zero$trial$response[, 1]),
                            (response_sim(dose = dose) +
                              max(as.numeric(step_zero$trial$response[, 1]))))
    event_time <- c(event_time, rep(NA, (sample_size - length(event_time))))
    current_time <- max(step_zero$trial$response[, 1])
    initial_time <- rep(0, sample_size)
    j <- length(dose)

    while ((current_time - initial_time[sample_size]) <= step_zero$trial$tau) {

      current_time <- current_time + rexp(1, 1)

      j <- j + 1

      if (j <= sample_size)
        initial_time[j:sample_size] <- current_time

      time_cens <- ifelse(event_time > (current_time - initial_time),
                          ifelse((current_time - initial_time) >
                                   step_zero$trial$tau, step_zero$trial$tau,
                                 (current_time - initial_time)),
                          ifelse(event_time > step_zero$trial$tau,
                                 step_zero$trial$tau, event_time))

      dlt <- ifelse(event_time > (current_time - initial_time), 0,
                    ifelse(event_time > step_zero$trial$tau, 0, 1))

      resolution <- ifelse(dlt == 1, 1,
                           ifelse((current_time - initial_time) >
                                    step_zero$trial$tau, 1, 0))

      if (j <= sample_size){
        alpha_sim[i, j] <- feasibility(alpha = alpha_sim[i, 1:(j-1)],
                                       strategy = alpha_strategy,
                                       dlt = dlt[1:(j-1)],
                                       resolution = resolution[1:(j-1)],
                                       rate = alpha_rate)

        formula <- cbind(time_cens[1:(j-1)], dlt[1:(j-1)]) ~ dose[1:(j-1)]
        update <- ewoc_d1ph(formula,
                            type = step_zero$trial$type,
                            theta = step_zero$trial$theta,
                            alpha = alpha_sim[i, j],
                            tau = step_zero$trial$tau,
                            min_dose = step_zero$trial$min_dose,
                            max_dose = step_zero$trial$max_dose,
                            first_dose = step_zero$trial$first_dose,
                            last_dose = step_zero$trial$last_dose,
                            dose_set = step_zero$trial$dose_set,
                            rho_prior = step_zero$trial$rho_prior,
                            mtd_prior = step_zero$trial$mtd_prior,
                            shape_prior = step_zero$trial$shape_prior,
                            distribution = step_zero$trial$distribution,
                            rounding = step_zero$trial$rounding)

        if (!is.null(stop_rule_sim))
          if (stop_rule_sim(update)){
            dose[j:sample_size] <- NA
            dlt[j:sample_size] <- NA
            event_time[j:sample_size] <- NA
            mtd_estimate <- NA
            rho_estimate <- NA
            break
          }

        dose[j] <- update$next_dose
        event_time[j] <- response_sim(dose = dose[j])
        mtd_estimate <- update$next_dose
        rho_estimate <- median(update$rho)
      }
    }

    dose_sim[i, ] <- dose
    dlt_sim[i, ] <- dlt
    time_sim[i, ] <- event_time
    mtd_sim[i] <- mtd_estimate
    rho_sim[i, ] <- rho_estimate
  }

  out <- list(time_sim = time_sim, dose_sim = dose_sim, dlt_sim = dlt_sim,
              mtd_sim = mtd_sim, rho_sim = rho_sim, alpha_sim = alpha_sim)
}


#'@export
trial_simulation.d1multinomial <-
  function(step_zero, n_sim, sample_size,
           alpha_strategy = c("fixed", "increasing", "conditional"),
           alpha_rate = 0.05,
           response_sim = NULL,
           covariable_sim = NULL,
           stop_rule_sim = NULL){

  if (is.null(response_sim))
    stop("'response_sim' function should be defined.")

  if (is.null(covariable_sim))
    stop("'covariable_sim' function should be defined.")

  n_dose <- sample_size + 1
  dose_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  covariable_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  dlt_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  mtd_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  rho_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  alpha_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)


  for (i in 1:n_sim){

    dlt <- as.numeric(step_zero$trial$response[, 1])
    npatients <- as.numeric(step_zero$trial$response[, 2])
    dose <- as.numeric(step_zero$trial$design_matrix[, 2])
    covariable <- as.numeric(step_zero$trial$covariable)
    alpha_sim[, 1:length(dose)] <- as.numeric(step_zero$trial$alpha)

    for (j in (length(dose)+1):n_dose) {

      formula <-
        cbind(dlt[1:(j-1)], npatients) ~ dose[1:(j-1)] | covariable[1:(j-1)]
      resolution <- ifelse(!is.na(dlt), 1, 0)

      if (j <= sample_size){
        alpha_sim[i, j] <- feasibility(alpha = alpha_sim[i, 1:(j-1)],
                                       strategy = alpha_strategy,
                                       dlt = dlt[1:(j-1)],
                                       resolution = resolution[1:(j-1)],
                                       rate = alpha_rate)
        covariable[j:(j + npatients)] <- covariable_sim(n = npatients)

        update <- ewoc_d1multinomial(formula,
                                     type = step_zero$trial$type,
                                     theta = step_zero$trial$theta,
                                     alpha = alpha_sim[i, j],
                                     min_dose = step_zero$trial$min_dose,
                                     max_dose = step_zero$trial$max_dose,
                                     first_dose = step_zero$trial$first_dose,
                                     last_dose = step_zero$trial$last_dose,
                                     dose_set = step_zero$trial$dose_set,
                                     levels_cov = step_zero$trial$levels_cov,
                                     next_patient_cov = covariable[j],
                                     rho_prior = step_zero$trial$rho_prior,
                                     mtd_prior = step_zero$trial$mtd_prior,
                                     rounding = step_zero$trial$rounding)

        if (!is.null(stop_rule_sim))
          if (stop_rule_sim(update)){
            dose[j:sample_size] <- NA
            dlt[j:sample_size] <- NA
            covariable[j:sample_size] <- NA
            mtd_estimate <- NA
            rho_estimate <- NA
            break
          }

        dose[j] <- update$next_dose
        dlt[j] <- response_sim(dose = dose[j], cov = covariable[j])
        npatients[j] <- npatients[j-1]
        mtd_estimate <- update$next_dose
        rho_estimate <- median(update$rho)
      }
    }

    dose_sim[i, ] <- dose
    dlt_sim[i, ] <- dlt
    covariable_sim[i, ] <- covariable
    mtd_sim[i] <- mtd_estimate
    rho_sim[i] <- rho_estimate
  }

  out <- list(dose_sim = dose_sim, covariable_sim = covariable_sim,
              dlt_sim = dlt_sim, mtd_sim = mtd_sim, rho_sim = rho_sim,
              alpha_sim = alpha_sim)
  return(out)
}
