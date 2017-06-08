#'@export
trial_simulation.d1ph <- function(step_zero, n_sim, sample_size,
                                  alpha_strategy =
                                    c("fixed", "increasing", "conditional"),
                                  alpha_rate = NULL,
                                  response_sim = NULL){

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
  alpha_sim <- matrix(NA, ncol = n_dose, nrow = n_sim)

  alpha_sim[, 1] <- step_zero$trial$alpha

  for (i in 1:n_sim){

    dlt <- rep(NA, sample_size)
    dlt[1] <- step_zero$trial$response[1, 1]
    dose <- rep(NA, sample_size)
    dose[1] <- step_zero$trial$design_matrix[1, 2]

    event_time <- rep(NA, sample_size)
    event_time[1] <- response_sim(dose = dose[1])

    current_time <- 0
    initial_time <- rep(0, sample_size)
    j <- 1

    while ((current_time - initial_time[sample_size]) <=
           step_zero$trial$tau) {

      current_time <- current_time + rexp(1, 1)

      if (j <= sample_size){
        j <- j + 1
        initial_time[j:sample_size] <- current_time
      }


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

      alpha_sim[i, j] <- feasibility(current_alpha = alpha_sim[i, (j - 1)],
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
                          distribution = step_zero$trial$distribution,
                          rounding = step_zero$trial$rounding)

      if (j <= sample_size){
        dose[j] <- update$next_dose
        event_time[j] <- response_sim(dose = dose[j])
      }
    }

    dose_sim[i, ] <- dose
    dlt_sim[i, ] <- dlt
    time_sim[i, ] <- time_cens
    mtd_sim[i] <- update$next_dose
    rho_sim[i, ] <- median(update$rho)
  }

  out <- list(time_sim = time_sim, dose_sim = dose_sim, dlt_sim = dlt_sim,
              mtd_sim = mtd_sim, rho_sim = rho_sim, alpha_sim = alpha_sim)
}
