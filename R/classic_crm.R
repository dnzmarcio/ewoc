#'Continual Reassessment Method
#'
#'Finding the next dose for a phase I clinical trial based on the
#'Continual Reassessment Method (CRM) design considering the
#'classic parametrization of EWOC for binary responses and single agent.
#'
#'@param formula an object of class \code{\link[Formula]{Formula}}: a symbolic
#'description of the model to be fitted with only one regressor term
#'corresponding to the dose for the right side and a numeric vector a response
#'containing number of DLT for the left side.
#'@param type a character describing the type of the Maximum Tolerable Dose
#'(MTD) variable.
#'@param theta a numerical value defining the proportion of expected patients
#'to experience a medically unacceptable, dose-limiting toxicity (DLT) if
#'administered the MTD.
#'@param mtd_prior a matrix 1 x 2 of hyperparameters for the Beta prior
#'distribution associated with the parameter MTD.
#'@param rho_prior a matrix 1 x 2 of hyperparameters for the Beta prior distribution
#'associated with the parameter rho.
#'@param min_dose a numerical value defining the lower bound of the support of
#'the MTD.
#'@param max_dose a numerical value defining the upper bound of the support of
#'the MTD.
#'@param first_dose a numerical value for the first allowable dose in the trial.
#'It is only necessary if type = 'continuous'.
#'@param last_dose a numerical value for the last allowable dose in the trial.
#'It is only necessary if type = 'continuous'.
#'@param dose_set a numerical vector of allowable doses in the trial.
#'It is only necessary if type = 'discrete'.
#'@param max_increment a numerical value indicating the maximum increment from the current dose to the next dose.
#'@param rounding a character indicating how to round a continuous dose to the
#'one of elements of the dose set. It is only necessary if type = 'discrete'.
#'@param n_adapt the number of iterations for adaptation.
#'See \code{\link[rjags]{adapt}} for details.
#'@param burn_in numerical value indicating the number of iterations before to start monitoring.
#'@param n_mcmc numerical value indicating the number of iterations to monitor.
#'@param n_thin numerical value corresponding to the thinning interval for monitors.
#'@param n_chains numerical value indicating the number of parallel chains for the model.
#'
#'@return \code{next_dose} the next recommend dose.
#'@return \code{mtd} the posterior MTD distribution.
#'@return \code{rho} the posterior rho_0 distribution.
#'@return \code{sample} a list of the MCMC chains distribution.
#'@return \code{trial} a list of the trial conditions.
#'
#'@references Babb, J., Rogatko, A. and Zacks, S., 1998.
#'Cancer phase I clinical trials: efficient dose escalation with overdose
#'control. Statistics in medicine, 17(10), pp.1103-1120.
#'
#'@examples
#'DLT <- 0
#'dose <- 20
#'test <- crm_d1classic(DLT ~ dose, type = 'discrete',
#'                       theta = 0.33,
#'                       min_dose = 0, max_dose = 100,
#'                       dose_set = seq(0, 100, 20),
#'                       rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                       mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                       rounding = "nearest")
#'summary(test)
#'plot(test)
#'
#'@import stats
#'
#'@export
crm_d1classic <- function(formula, theta,
                           mtd_prior, rho_prior,
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

  if (!(theta > 0 & theta < 1))
    stop("'theta' should be in the interval (0, 1).")

  if (nrow(mtd_prior) != 1 | ncol(mtd_prior) != 2)
    stop(paste0("'mtd_prior' should be a matrix with 2 columns and 1 row."))

  if (nrow(rho_prior) != 1 | ncol(rho_prior) != 2)
    stop(paste0("'rho_prior' should be a matrix with 2 columns and 1 row."))

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
                  theta = theta, limits = limits,
                  dose_set = dose_set,
                  max_increment = max_increment, current_dose = current_dose,
                  rho_prior = rho_prior, mtd_prior = mtd_prior,
                  type = type[1], rounding = rounding)
  class(my_data) <- c("crm_d1classic", "d1classic")

  my_data$mcmc <- jags(my_data, n_adapt, burn_in, n_mcmc, n_thin, n_chains)
  out <- next_dose(my_data)

  design_matrix[, 2] <-
    inv_standard_dose(dose = design_matrix[, 2],
                      min_dose = limits$min_dose,
                      max_dose = limits$max_dose)

  trial <- list(response = response, design_matrix = design_matrix,
                theta = theta,
                first_dose = limits$first_dose, last_dose = limits$last_dose,
                min_dose = limits$min_dose, max_dose = limits$max_dose,
                dose_set = dose_set,  max_increment = max_increment,
                rho_prior = rho_prior, mtd_prior = mtd_prior,
                type = type, rounding = rounding,
                n_adapt = n_adapt, burn_in = burn_in, n_mcmc = n_mcmc,
                n_thin = n_thin, n_chains = n_chains)
  out$trial <- trial

  class(out) <- c("ewoc_d1classic", "d1classic")
  return(out)
}
