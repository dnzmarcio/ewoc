#'Escalation With Overdose Control
#'
#'Finding the next dose for a phase I clinical trial based on the
#'Escalation with Overdose Control (EWOC) design considering the
#'extended parametrization for binary response and single agent.
#'
#'@param formula an object of class \code{\link[Formula]{Formula}}: a symbolic
#'description of the model to be fitted with only one regressor term
#'corresponding to the dose for the right side and a numeric vector as a response
#'containing number of DLT for the left side.
#'@param theta a numerical value defining the proportion of expected patients
#'to experience a medically unacceptable, dose-limiting toxicity (DLT) if
#'administered the MTD.
#'@param alpha a numerical value defining the probability that the dose selected
#'by EWOC is higher than the MTD.
#'@param rho_prior a matrix 3 x 2 of hyperparameters for the Beta prior
#'distribution associated with each parameter rho. Each row corresponds to a parameter.
#'@param type a character describing the type of the Maximum Tolerable Dose
#'(MTD) variable.
#'@param min_dose a numerical value defining the lower bound of the support of
#'the MTD.
#'@param max_dose a numerical value defining the upper bound of the support of
#'the MTD.
#'@param first_dose a numerical value for the first allowable dose in the trial.
#'@param last_dose a numerical value for the last allowable dose in the trial.
#'@param dose_set a numerical vector of allowable doses in the trial. It is only
#'necessary if type = "discrete".
#'@param max_increment a numerical value indicating the maximum increment from the current dose to the next dose.
#'@param rounding a character indicating how to round a continuous dose to the
#'one of elements of the dose set.
#'It is only necessary if type = "discrete".
#'@param n_adapt the number of iterations for adaptation.
#'See \code{\link[rjags]{adapt}} for details.
#'@param burn_in the number of iterations before to start monitoring.
#'@param n_mcmc the number of iterations to monitor.
#'@param n_thin thinning interval for monitors.
#'@param n_chains the number of parallel chains for the model.
#'
#'@return \code{next_dose} the next recommend dose.
#'@return \code{mtd} a numerical vector for the posterior MTD distribution considering the next patient covariable.
#'@return \code{rho} a matrix for the posterior rho_0 and rho_1 distributions.
#'@return \code{sample} a list of the MCMC chains distribution.
#'@return \code{trial} a list of the trial conditions.
#'
#'@examples
#'DLT <- 0
#'dose <- 30
#'
#'test <- ewoc_d1extended(DLT ~ dose, type = 'discrete',
#'                        theta = 0.33, alpha = 0.25,
#'                        dose_set = c(30, 40, 50),
#'                        min_dose = 20, max_dose = 100,
#'                        rho_prior = matrix(1, ncol = 2, nrow = 2),
#'                        rounding = "nearest")
#'summary(test)
#'plot(test)
#'
#'@import stats
#'
#'@export
ewoc_d1extended <- function(formula, theta, alpha,
                            rho_prior,
                            min_dose, max_dose,
                            type = c('continuous', 'discrete'),
                            first_dose = NULL, last_dose = NULL,
                            dose_set = NULL, max_increment = NULL,
                            rounding = c("down", "nearest"),
                            n_adapt = 5000, burn_in = 1000,
                            n_mcmc = 1000, n_thin = 1, n_chains = 1) {

  formula <- Formula::Formula(formula)
  if (class(formula)[2] != "formula")
    stop("Invalid formula! \n")

  data_base <- model.frame(formula, na.action = na.exclude,
                                  drop.unused.levels = FALSE)

  dose_matrix <- model.matrix(formula, data_base, rhs = 1)

  if (length(formula)[2] == 1){
    design_matrix <- dose_matrix
    colnames(design_matrix) <- c("intercept", "dose")
  } else {
    stop("This design cannot accommodate a covariable.")
  }
  colnames(design_matrix) <- c("intercept", "dose")

  response <- model.response(data_base)

  if (length(type) > 1 | !(type == "continuous" | type == "discrete"))
    stop("'type' should be either 'continuous' or 'discrete'.")

  if (type == "discrete") {
    if (is.null(dose_set))
      stop("'dose_set' should be informed for type = 'discrete'.")

    if (length(rounding) > 1 | !(rounding == "down" | rounding == "nearest"))
      stop("'rounding' should be either 'down' or 'nearest'.")

    if (is.null(max_increment))
      max_increment <- max(diff(dose_set))
  }

  if (!(alpha > 0 & alpha < 1))
    stop("'alpha' should be in the interval (0, 1).")

  if (!(theta > 0 & theta < 1))
    stop("'theta' should be in the interval (0, 1).")

  if (nrow(rho_prior) != 2 | ncol(rho_prior) != 2)
    stop(paste0("'rho_prior' should be a matrix with 2 columns and 2 rows."))

  limits <- limits_d1nocov(first_dose = first_dose, last_dose = last_dose,
                           min_dose = min_dose, max_dose = max_dose,
                           type = type, rounding = rounding,
                           dose_set = dose_set)

  if (is.null(max_increment))
    max_increment <- limits$last_dose - limits$first_dose

  current_dose <- design_matrix[nrow(design_matrix), 2]

  design_matrix[, 2] <-
    standard_dose(dose = design_matrix[, 2],
                  min_dose = limits$min_dose,
                  max_dose = limits$max_dose)

  my_data <- list(response = response, design_matrix = design_matrix,
                  theta = theta, alpha = alpha,
                  limits = limits,
                  dose_set = dose_set, current_dose = current_dose,
                  rho_prior = rho_prior,
                  type = type, rounding = rounding)
  class(my_data) <- c("ewoc_d1extended", "d1extended")

  my_data$mcmc <- jags(my_data, n_adapt, burn_in, n_mcmc, n_thin, n_chains)
  out <- next_dose(my_data)

  trial <- list(response = response, design_matrix = design_matrix,
                theta = theta, alpha = alpha,
                first_dose = limits$first_dose, last_dose = limits$last_dose,
                min_dose = limits$min_dose, max_dose = limits$max_dose,
                dose_set = dose_set,
                rho_prior = rho_prior,
                type = type, rounding = rounding,
                n_adapt = n_adapt, burn_in = burn_in, n_mcmc = n_mcmc,
                n_thin = n_thin, n_chains = n_chains)
  out$trial <- trial

  class(out) <- c("ewoc_d1extended", "d1extended")

  return(out)
}




