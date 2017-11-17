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
#'Its only input is 'dose' and output is the indicator of DLT for classical and
#'extended EWOC and the time until DLT for PH EWOC.
#'@param stop_rule_sim a function having as an input an object containing all
#'the information related to the trial as the returned object trial from either
#'\code{ewoc_d1classic}, \code{ewoc_d1extended}, \code{ewoc_d1ph} and as
#'output a logical value indicating the trial should be stopped.
#'@param ncores a numeric value indicating the number of cores to be used in the
#'simulation performed in parallel.
#'
#'@return \code{alpha_sim} a matrix \code{n_sim} x \code{sample_size} containing
#'the values of feasibility used for each step in the trial and each trial in
#'the simulation.
#'@return \code{dlt_sim} a matrix \code{n_sim} x \code{sample_size} containing
#'ones and zeros indicating the occurrence of DLT (1) and the absence of DLT (0)
#'for each step in the trial and each trial in the simulation.
#'@return \code{dose_sim} a matrix \code{n_sim} x \code{sample_size} containing
#'the doses assigned for each step in the trial and each trial in the simulation.
#'@return \code{mtd_sim} a numeric vector \code{n_sim} x 1 containing
#'the recommended MTD for each trial in the simulation.
#'@return \code{rho_sim} a numeric vector \code{n_sim} x k containing
#'the estimated rho parameter(s) for each trial in the simulation, where k = 1
#'for ewoc_d1classic, ewoc_d1ph, and k = 2 for ewoc_d1extended.
#'
#'@examples
#'\dontshow{
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
#'sim <- ewoc_simulation(step_zero = step_zero,
#'                       n_sim = 1, sample_size = 2,
#'                       alpha_strategy = "increasing",
#'                       response_sim = response_sim,
#'                       ncores = 2)
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
#'response_sim <- response_d1extended(rho = c(0.05, 0.5),
#'                                    min_dose = 10, max_dose = 50)
#'sim <- ewoc_simulation(step_zero = step_zero,
#'                       n_sim = 1, sample_size = 2,
#'                       alpha_strategy = "increasing",
#'                       response_sim = response_sim,
#'                       ncores = 2)
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
#'sim <- ewoc_simulation(step_zero = step_zero,
#'                       n_sim = 1, sample_size = 2,
#'                       alpha_strategy = "increasing",
#'                       response_sim = response_sim,
#'                       ncores = 2)
#'}
#'
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
#'sim <- ewoc_simulation(step_zero = step_zero,
#'                        n_sim = 2, sample_size = 30,
#'                        alpha_strategy = "increasing",
#'                        response_sim = response_sim,
#'                        ncores = 2)
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
#'response_sim <- response_d1extended(rho = c(0.05, 0.5),
#'                                    min_dose = 10, max_dose = 50)
#'sim <- ewoc_simulation(step_zero = step_zero,
#'                        n_sim = 2, sample_size = 30,
#'                        alpha_strategy = "increasing",
#'                        response_sim = response_sim,
#'                        ncores = 2)
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
#'sim <- ewoc_simulation(step_zero = step_zero,
#'                        n_sim = 2, sample_size = 30,
#'                        alpha_strategy = "increasing",
#'                        response_sim = response_sim,
#'                        ncores = 2)
#'}
#'
#'@importFrom foreach foreach %dopar%
#'@importFrom doParallel registerDoParallel
#'
#'@export
ewoc_simulation <- function(step_zero, n_sim, sample_size,
                            alpha_strategy = "fixed",
                            alpha_rate = NULL, response_sim,
                            stop_rule_sim = NULL,
                            ncores = 1){
  UseMethod("ewoc_simulation")
}


#'@export
ewoc_simulation.ewoc_d1classic <- function(step_zero, n_sim, sample_size,
                                      alpha_strategy =
                                        c("fixed", "increasing", "conditional"),
                                      alpha_rate = 0.05,
                                      response_sim = NULL,
                                      stop_rule_sim = NULL,
                                      ncores = 1){

  if (is.null(response_sim))
    stop("'response_sim' function should be defined.")

  n_dose <- sample_size + 1
  dose_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  dlt_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  mtd_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  rho_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  alpha_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)

  registerDoParallel(ncores)
  result <-
    foreach(i = 1:n_sim,
            .combine='comb',
            .multicombine=TRUE,
            .init=list(list(), list(), list(), list(), list())) %dopar% {

              dlt <- as.numeric(step_zero$trial$response)
              dose <- as.numeric(step_zero$trial$design_matrix[, 2])
              alpha <- as.numeric(step_zero$trial$alpha)

              for (j in (length(dose)+1):n_dose) {

                formula <- dlt[1:(j-1)] ~ dose[1:(j-1)]
                resolution <- ifelse(!is.na(dlt), 1, 0)

                if (j <= sample_size){
                  alpha[j] <- feasibility(alpha = alpha[(j-1)],
                                          strategy = alpha_strategy,
                                          dlt = dlt[1:(j-1)],
                                          resolution = resolution[1:(j-1)],
                                          rate = alpha_rate)

                  update <- ewoc_d1classic(formula,
                                           type = step_zero$trial$type,
                                           theta = step_zero$trial$theta,
                                           alpha = alpha[j],
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

              list(dose, dlt, mtd_estimate, rho_estimate, alpha)
            }


  dose_sim <- result[[1]]
  dlt_sim <- result[[2]]
  mtd_sim <- result[[3]]
  rho_sim <- result[[4]]
  alpha_sim <- result[[5]]

  out <- list(dose_sim = dose_sim, dlt_sim = dlt_sim,
              mtd_sim = mtd_sim, rho_sim = rho_sim, alpha_sim = alpha_sim)
  return(out)
}

#'@importFrom foreach foreach %dopar%
#'@importFrom doParallel registerDoParallel
#'@export
ewoc_simulation.ewoc_d1extended <- function(step_zero, n_sim, sample_size,
                                       alpha_strategy =
                                         c("fixed", "increasing", "conditional"),
                                       alpha_rate = 0.05,
                                       response_sim = NULL,
                                       stop_rule_sim = NULL,
                                       ncores = 1){

  if (is.null(response_sim))
    stop("'response_sim' function should be defined.")

  n_dose <- sample_size + 1
  dose_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  dlt_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)
  mtd_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  rho_sim <- matrix(NA, ncol = 1, nrow = n_sim)
  alpha_sim <- matrix(NA, ncol = sample_size, nrow = n_sim)


  registerDoParallel(ncores)
  result <-
    foreach(i = 1:n_sim,
            .combine='comb',
            .multicombine=TRUE,
            .init=list(list(), list(), list(), list(), list())) %dopar% {

              dlt <- as.numeric(step_zero$trial$response)
              dose <- as.numeric(step_zero$trial$design_matrix[, 2])
              alpha <- as.numeric(step_zero$trial$alpha)

              for (j in (length(dose)+1):n_dose) {

                formula <- dlt[1:(j-1)] ~ dose[1:(j-1)]
                resolution <- ifelse(!is.na(dlt), 1, 0)

                if (j <= sample_size){
                  alpha[j] <- feasibility(alpha = alpha[(j-1)],
                                              strategy = alpha_strategy,
                                              dlt = dlt[1:(j-1)],
                                              resolution = resolution[1:(j-1)],
                                              rate = alpha_rate)

                  update <- ewoc_d1extended(formula,
                                            type = step_zero$trial$type,
                                            theta = step_zero$trial$theta,
                                            alpha = alpha[j],
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

              list(dose, dlt, mtd_estimate, rho_estimate, alpha)
            }

  dose_sim <- result[[1]]
  dlt_sim <- result[[2]]
  mtd_sim <- result[[3]]
  rho_sim <- result[[4]]
  alpha_sim <- result[[5]]

  out <- list(dose_sim = dose_sim, dlt_sim = dlt_sim,
              mtd_sim = mtd_sim, rho_sim = rho_sim, alpha_sim = alpha_sim)
  return(out)
}

#'@importFrom foreach foreach %dopar%
#'@importFrom doParallel registerDoParallel
#'@export
ewoc_simulation.ewoc_d1ph <- function(step_zero, n_sim, sample_size,
                                 alpha_strategy =
                                   c("fixed", "increasing", "conditional"),
                                 alpha_rate = 0.05,
                                 response_sim = NULL,
                                 stop_rule_sim = NULL,
                                 ncores = 1){

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

  registerDoParallel(ncores)
  result <-
    foreach(i = 1:n_sim,
            .combine='comb',
            .multicombine=TRUE,
            .init=list(list(), list(), list(), list(), list(), list())) %dopar% {

              dlt <- as.numeric(step_zero$trial$response[, 2])
              dose <- as.numeric(step_zero$trial$design_matrix[, 2])
              alpha <- as.numeric(step_zero$trial$alpha)

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
                  alpha[j] <- feasibility(alpha = alpha[(j-1)],
                                                 strategy = alpha_strategy,
                                                 dlt = dlt[1:(j-1)],
                                                 resolution = resolution[1:(j-1)],
                                                 rate = alpha_rate)

                  formula <- cbind(time_cens[1:(j-1)], dlt[1:(j-1)]) ~ dose[1:(j-1)]
                  update <- ewoc_d1ph(formula,
                                      type = step_zero$trial$type,
                                      theta = step_zero$trial$theta,
                                      alpha = alpha[j],
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

              list(event_time, dose, dlt, mtd_estimate, rho_estimate, alpha)
            }

  time_sim <- result[[1]]
  dose_sim <- result[[2]]
  dlt_sim <- result[[3]]
  mtd_sim <- result[[4]]
  rho_sim <- result[[5]]
  alpha_sim <- result[[6]]

  out <- list(time_sim = time_sim, dose_sim = dose_sim, dlt_sim = dlt_sim,
              mtd_sim = mtd_sim, rho_sim = rho_sim, alpha_sim = alpha_sim)
}

