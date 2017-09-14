#'Escalation Over with Dose Control
#'
#'Finding the next dose for a phase I clinical trial based on Escalation Over
#'Dose Control (EWOC) design considering extended parametrization
#'for binary response with continuous covariable.
#'
#'@param formula an object of class \code{\link[Formula]{Formula}}: a symbolic
#'description of the model to be fitted with two regressor parts separated by `|`
#'corresponding to the dose and covariable, respectively, for the right side  and
#'a numeric vector as a response containing number of DLT for the left side.
#'@param theta a numerical value defining the proportion of expected patients
#'to experience a medically unacceptable, dose-limiting toxicity (DLT) if
#'administered the MTD.
#'@param alpha a numerical value defining the probability that the dose selected
#'by EWOC is higher than the MTD.
#'@param rho_prior a matrix of hyperparameters for the Beta prior distribution
#'associated with each rho. Each row corresponds to a paramater.
#'@param min_cov a numerical value defining the minimum value for the covariate.
#'@param max_cov a numerical value defining the maximum value for the covariate.
#'@param next_patient_cov a character value indicating the covariable associated to
#'next patient.
#'@param direction a character defining the relationship between the covariate and the
#'probability of DLT.
#'@param min_dose either a numerical value or a function of the covariable
#'defining the lower bound of the support of the MTD used to standardize the doses.
#'@param max_dose either a numerical value or a function of the covariable
#'defining the upper bound of the support of the MTD used to standardize the doses.
#'@param type a character describing the type of the Maximum Tolerable Dose
#'(MTD) variable.
#'@param first_dose a numerical value for the first allowable dose in the trial.
#'It is only necessary if type = 'continuous'.
#'@param last_dose a numerical value for the last allowable dose in the trial.
#'It is only necessary if type = 'continuous'.
#'@param dose_set a numerical vector of allowable doses in the trial. It is only
#'necessary if type = 'discrete'.
#'@param rounding a character indicating how to round a continuous dose to the
#'one of elements of the dose set. It is only necessary if type = 'discrete'.
#'@param n_adapt the number of iterations for adaptation.
#'See \code{\link[rjags]{adapt}} for details.
#'@param burn_in the number of iterations before to start monitoring.
#'@param n_mcmc the number of iterations to monitor.
#'@param n_thin thinning interval for monitors.
#'@param n_chains the number of parallel chains for the model.
#'
#'@return \code{next_dose} a numerical value corresponding to the next recommend dose.
#'@return \code{mtd} a numerical vector for the posterior MTD distribution considering the next patient covariable.
#'@return \code{rho} a matrix for the posterior rho_00 and rho_01 distributions.
#'@return \code{gamma} a numerical vector for the posterior standardized MTD distribution considering the next patient covariable.
#'@return \code{sample} a list of the MCMC chains distribution.
#'@return \code{trial} a list of the trial conditions.
#'
#'@references Babb JS, Rogatko A. Patient specific dosing in a cancer phase I clinical trial. Statistics in medicine. 2001 Jul 30;20(14):2079-90.
#'
#'@export
ewoc_d1excontinuous <- function(formula, theta, alpha,
                                rho_prior,
                                min_dose, max_dose,
                                min_cov, max_cov, next_patient_cov,
                                direction = c('positive', 'negative'),
                                type = c('continuous', 'discrete'),
                                first_dose = NULL, last_dose = NULL,
                                dose_set = NULL,
                                rounding = c("down", "nearest"),
                                n_adapt = 5000, burn_in = 1000,
                                n_mcmc = 1000, n_thin = 1, n_chains = 1) {

  formula <- Formula::Formula(formula)
  if (class(formula)[2] != "formula")
    stop("Invalid formula! \n")

  data_base <- model.frame(formula, na.action = na.exclude,
                           drop.unused.levels = FALSE)

  dose_matrix <- model.matrix(formula, data_base, rhs = 1)

  if (length(formula)[2] == 1) {
    stop("This design requires a continuous covariable.")
  } else {
    covariable_matrix <- model.matrix(formula, data_base, rhs = 2)
    covariable <- covariable_matrix[, 2]
  }

  design_matrix <- cbind(dose_matrix, covariable)
  colnames(design_matrix) <- c("intercept", "dose", "covariable")
  response <- model.response(data_base)

  if (length(type) > 1 | !(type == "continuous" | type == "discrete"))
    stop("'type' should be either 'continuous' or 'discrete'.")

  if (type == "discrete") {
    if (is.null(dose_set))
      stop("'dose_set' should be informed for type = 'discrete'.")

    if (length(rounding) > 1 | !(rounding == "down" | rounding == "nearest"))
      stop("'rounding' should be either 'down' or 'nearest'.")
  }

  if (!(alpha > 0 & alpha < 1))
    stop("'alpha' should be in the interval (0, 1).")

  if (!(theta > 0 & theta < 1))
    stop("'theta' should be in the interval (0, 1).")

  if (any(covariable < min_cov | covariable > max_cov))
    stop("'covariable' in formula has to be between 'min_cov' and 'max_cov'.")

  if (next_patient_cov < min_cov | next_patient_cov > max_cov)
    stop("'next_patient_cov' has to be between 'min_cov' and 'max_cov'.")

  limits <- limits_d1cov(first_dose = first_dose, last_dose = last_dose,
                         min_dose = min_dose, max_dose = max_dose,
                         type = type, rounding = rounding,
                         dose_set = dose_set,
                         covariable = covariable)

  lb <- - limits$min_dose(next_patient_cov)/(limits$max_dose(next_patient_cov) -
                                         limits$min_dose(next_patient_cov))

  design_matrix[, 2] <- standard_dose(dose = design_matrix[, 2],
                                      min_dose = limits$min_dose(covariable),
                                      max_dose = limits$max_dose(covariable))

  design_matrix[, 3] <- standard_dose(dose = design_matrix[, 3],
                                      min_dose = min_cov,
                                      max_dose = max_cov)

  my_data <- list(response = response, design_matrix = design_matrix,
                  next_patient_cov = next_patient_cov,
                  direction = direction,
                  theta = theta, alpha = alpha, limits = limits,
                  dose_set = dose_set,
                  rho_prior = rho_prior,
                  max_cov = max_cov, min_cov = min_cov, lb = lb,
                  type = type, rounding = rounding)
  class(my_data) <- "d1excontinuous"

  out <- qmtd_jags(my_data, n_adapt, burn_in, n_mcmc, n_thin, n_chains)

  trial <- list(response = response, design_matrix = design_matrix,
                theta = theta, alpha = alpha,
                first_dose = limits$first_dose, last_dose = limits$last_dose,
                min_dose = limits$min_dose, max_dose = limits$max_dose,
                dose_set = dose_set,
                rho_prior = rho_prior,
                covariable = covariable, next_patient_cov = next_patient_cov,
                direction = direction,
                min_cov = min_cov, max_cov = max_cov,
                type = type, rounding = rounding,
                n_adapt = n_adapt, burn_in = burn_in, n_mcmc = n_mcmc,
                n_thin = n_thin, n_chains = n_chains)
  out$trial <- trial

  class(out) <- c("ewoc_d1excontinuous", "d1excontinuous")
  return(out)
}

#'@export
ewoc_jags.d1excontinuous <- function(data, n_adapt, burn_in,
                                     n_mcmc, n_thin, n_chains) {

  # JAGS model function
  if(data$direction == "positive"){
    jfun <- "model {

      for(i in 1:nobs) {
        dlt[i] ~ dbin(p[i], 1)
        p[i] <- 1/(1 + exp(-lp[i]))
        lp[i] <- inprod(design_matrix[i, ], beta)
      }

      beta[1] <- logit(rho[1])
      beta[2] <- logit(rho[3]) - logit(rho[1])
      beta[3] <- logit(rho[2]) - logit(rho[1])

      rho[1] <- min*r[1]
      min <- min(exp.limit, rho[2], rho[3])
      r[1] ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])
      rho[2] <- theta*r[2]
      r[2] ~ dbeta(rho_prior[2, 1], rho_prior[2, 2])
      rho[3] <- r[3]
      r[3] ~ dbeta(rho_prior[3, 1], rho_prior[3, 2])

      exp.limit <- plogis(limit, 0, 1)
      limit <- (logit(theta) -
                ((next_patient_cov - min_cov)/(max_cov - min_cov))*logit(rho[2]) -
                lb*logit(rho[3]))/
                (1 - lb - (next_patient_cov - min_cov)/(max_cov - min_cov))

    }"
  } else {

    jfun <- "model {

      for(i in 1:nobs) {
        dlt[i] ~ dbin(p[i], 1)
        p[i] <- 1/(1 + exp(-lp[i]))
        lp[i] <- inprod(design_matrix[i, ], beta)
      }

      beta[1] <- logit(rho[1])
      beta[2] <- logit(rho[3]) - logit(rho[1])
      beta[3] <- logit(rho[2]) - logit(rho[1])

      rho[1] <- min01*r[1]
      min01 <- min(exp.limit, rho[3])
      r[1] ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])
      rho[2] <- min02*r[2]
      min02 <- min(rho[2], theta)
      r[2] ~ dbeta(rho_prior[2, 1], rho_prior[2, 2])
      rho[3] <- r[3]
      r[3] ~ dbeta(rho_prior[3, 1], rho_prior[3, 2])

      exp.limit <- plogis(limit, 0, 1)
      limit <- (logit(theta) -
                ((next_patient_cov - min_cov)/(max_cov - min_cov))*logit(rho[2]) -
                lb*logit(rho[3]))/
                (1 - lb - (next_patient_cov - min_cov)/(max_cov - min_cov))
    }"

  }

  data_base <- list('dlt' = data$response,
                    'design_matrix' = data$design_matrix,
                    'theta' = data$theta, 'nobs' = length(data$response),
                    'rho_prior' = data$rho_prior,
                    'lb' = data$lb, 'next_patient_cov' = data$next_patient_cov,
                    'min_cov' = data$min_cov, 'max_cov' = data$max_cov)

  inits <- function() {
    out <- list(r = rbeta(nrow(data$rho_prior),
                          data$rho_prior[, 1], data$rho_prior[, 2]))
    return(out)
  }

  # Calling JAGS
  j <- rjags::jags.model(textConnection(jfun),
                         data = data_base,
                         inits = inits(),
                         n.chains = n_chains,
                         n.adapt = n_adapt)
  update(j, burn_in)
  sample <- rjags::coda.samples(j, variable.names = c("beta", "rho"),
                                n.iter = n_mcmc, thin = n_thin,
                                n.chains = n_chains)
  beta <- sample[[1]][, 1:3]
  rho <- sample[[1]][, 4:6]

  out <- list(rho = rho, beta = beta, sample = sample)

  return(out)
}
