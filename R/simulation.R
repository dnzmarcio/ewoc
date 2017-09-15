#'Generic EWOC simulation
#'
#'Generic function for simulating EWOC trials.
#'
#'@param step_zero an object from the classes 'ewoc_d1classic', 'ewoc_d1extended',
#''ewoc_d1ph' created using dummy data.
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
#'Its only imput is 'dose' and output is the indicator of DLT for classical and
#'extended EWOC and the time until DLT for PH EWOC.
#'@param stop_rule_sim a function having as an imput an object containing all
#'the information related to the trial as the returned object trial from either
#'\code{ewoc_d1classic}, \code{ewoc_d1extended}, \code{ewoc_d1ph} and as
#'output a logical valuel indicating the trial should be stopped.
#'
#'@examples
#'\dontrun{
#'### Classic EWOC
#'DLT <- 0
#'dose <- 30
#'step_zero <- ewoc_d1classic(DLT ~ dose, type = 'discrete',
#'                            theta = 0.33, alpha = 0.25,
#'                            min_dose = 0, max_dose = 100,
#'                            dose_set = seq(0, 100, 20),
#'                            rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                            mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                            rounding = "nearest")
#'response_sim <- response_d1classic(rho = 0.05, mtd = 20, theta = 0.33,
#'                                   min_dose = 10, max_dose = 50)
#'sim <- trial_simulation(step_zero = step_zero,
#'                        n_sim = 1, sample_size = 30,
#'                        alpha_strategy = "constant",
#'                        response_sim = response_sim)
#'
#'### Extended EWOC
#'DLT <- 0
#'dose <- 30
#'step_zero <- ewoc_d1extended(DLT ~ dose, type = 'discrete',
#'                            theta = 0.33, alpha = 0.25,
#'                            min_dose = 0, max_dose = 100,
#'                            dose_set = seq(0, 100, 20),
#'                            rho_prior = matrix(1, ncol = 2, nrow = 2),
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
#'                              min_dose = 30, max_dose = 50,
#'                              tau = 10, distribution = "exponential")
#'sim <- trial_simulation(step_zero = step_zero,
#'                        n_sim = 1, sample_size = 30,
#'                        alpha_strategy = "increasing",
#'                        response_sim = response_sim)
#'}
#'
#'@export
trial_simulation <- function(step_zero, n_sim, sample_size,
                             alpha_strategy = "fixed",
                             alpha_rate = NULL, response_sim,
                             stop_rule_sim = NULL){
  UseMethod("trial_simulation")
}


#'Simulation of trials using classical EWOC
#'
#'Performing a simulation of several phase I clinical trials based on classical EWOC.
#'
#'@param step_zero an object from the classes 'ewoc_d1classic'
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
#'Its only imput is 'dose' and output is the ndicator of DLT.
#'@param stop_rule_sim a function having as an imput an object containing all
#'the information related to the trial as the returned object trial from
#'\code{ewoc_d1classic} and as output a logical valuel indicating the trial
#'should be stopped.
#'
#'@export
trial_simulation.d1classic <- function(step_zero, n_sim, sample_size,
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

        update <- ewoc_d1classic(formula,
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

#'Simulation of trials using extended EWOC
#'
#'Performing a simulation of several phase I clinical trials based on the extended EWOC.
#'
#'@param step_zero an object from the class 'ewoc_d1extended'
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
#'Its only imput is 'dose' and output is the indicator of DLT.
#'@param stop_rule_sim a function having as an imput an object containing all
#'the information related to the trial as the returned object trial from
#'\code{ewoc_d1extended} and as output a logical valuel indicating the trial
#'should be stopped.
#'
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

#'Simulation of trials using proportional hazards EWOC
#'
#'Performing a simulation of several phase I clinical trials based on the proportional hazards EWOC.
#'
#'@param step_zero an object from the class 'ewoc_d1ph'
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
#'Its only imput is 'dose' and output is the time of DLT.
#'@param stop_rule_sim a function having as an imput an object containing all
#'the information related to the trial as the returned object trial from
#'\code{ewoc_d1ph} and as output a logical valuel indicating the trial
#'should be stopped.
#'
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


