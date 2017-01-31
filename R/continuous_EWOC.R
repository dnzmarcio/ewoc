#'Escalation Over with Dose Control
#'
#'Finding the next dose for a phase I clinical trial based on Escalation Over
#'Dose Control (EWOC) design considering extended parametrization
#'for binary response with continuous covariable.
#'
#'@param formula an object of class \code{\link[Formula]{Formula}}: a symbolic
#'description of the model to be fitted with two regressor parts separated by `|`
#'corresponding to the dose and covariable, respectively, for the right side and
#'a matrix as a response containing number of DLT and number of patients for
#'the left side.
#'@param theta a numerical value defining the proportion of expected patients
#'to experience a medically unacceptable, dose-limiting toxicity (DLT) if
#'administered the MTD.
#'@param alpha a numerical value defining the probability that the dose selected
#'by EWOC is higher than the MTD.
#'@param rho_prior a matrix of hyperparameters for the Beta prior distribution
#'associated with each rho. Each row corresponds to a paramater.
#'@param mtd_prior a matrix of hyperparameters for the Beta prior distribution
#'associated with the MTD. Each row corresponds to a paramater.
#'@param next_patient_cov a character value indicating the covariable associated to
#'next patient.
#'@param levels_cov a character vector of the possible values for the ordinal covariable.
#'@param first_dose a numerical value for the first allowable dose in the trial.
#'It is only necessary if type = 'continuous'.
#'@param last_dose a numerical value for the last allowable dose in the trial.
#'It is only necessary if type = 'continuous'.
#'@param dose_set a numerical vector of allowable doses in the trial. It is only
#'necessary if type = 'discrete'.
#'@param min_dose a numerical value defining the lower bound of the support of
#'the MTD.
#'@param max_dose a numerical value defining the upper bound of the support of
#'the MTD.
#'@param type a character describing the type of the Maximum Tolerable Dose
#'(MTD) variable.
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
#'@return \code{hpd_dose} a numerical vector containing  the 95\% HPD for the next dose.
#'@return \code{pdlt} a numerical value corresponding to the probability of DLT for the next dose.
#'@return \code{hpd_pdlt} a numerical vector containing  the 95\% HPD for the probability of DLT for the next dose.
#'@return \code{mtd} a numerical vector for the posterior MTD distribution considering the next patient covariable.
#'@return \code{rho} a matrix for the posterior rho_00 and rho_01 distributions.
#'@return \code{gamma} a numerical vector for the posterior standardized MTD distribution considering the maximum value of the covariable.
#'@return \code{sample} a list of the MCMC chains distribution.
#'@return \code{trial} a list of trial conditions.
#'
#'@examples
#'DLT <- 0
#'npatients <- 1
#'dose <- 30
#'tumor_size <- 2
#'test <- ewoc1d_continuous(cbind(DLT, npatients) ~ dose | tumor_size,
#'                          type = 'continuous',
#'                          theta = 0.33, alpha = 0.25, direction = 'positive',
#'                          dose_set = seq(30, 50, 5),
#'                          first_dose = 30, last_dose = 50,
#'                          min_dose = 30, max_dose = 50,
#'                          next_patient_cov = 2, min_cov = 1, max_cov = 4,
#'                          mtd_prior = matrix(1, nrow = 1, ncol = 2),
#'                          rho_prior = matrix(1, nrow = 2, ncol = 2))
#'summary(test)
#'
#'@references Babb JS, Rogatko A. Patient specific dosing in a cancer phase I clinical trial. Statistics in medicine. 2001 Jul 30;20(14):2079-90.
#'
#'@export
ewoc_d1continuous <- function(formula, theta, alpha,
                                mtd_prior, rho_prior,
                                next_patient_cov, min_cov, max_cov,
                                direction = c('positive', 'negative'),
                                type = c('continuous', 'discrete'),
                                first_dose = NULL, last_dose = NULL,
                                dose_set = NULL,
                                min_dose = NULL, max_dose = NULL,
                                rounding = c("down", "nearest"),
                                n_adapt = 5000, burn_in = 1000,
                                n_mcmc = 1000, n_thin = 1, n_chains = 1) {

  formula <- Formula(formula)
  if (class(formula)[2] != "formula")
    stop("Invalid formula! \n")

  data_base <- model.frame(formula, na.action = na.exclude,
                           drop.unused.levels = FALSE)

  dose_matrix <- model.matrix(formula, data_base, rhs = 1)

  if (length(formula)[2] == 1) {
    stop("This design requires a multinomial covariable.")
  } else {
    covariable_matrix <- model.matrix(formula, data_base, rhs = 2)
    covariable <- covariable_matrix[, 2]

    if (direction == "negative") {
      covariable <- -covariable_matrix
      min_cov <- -max_cov
      max_cov <- -min_cov
    }
  }

  design_matrix <- cbind(dose_matrix, covariable)
  colnames(design_matrix) <- c("intercept", "dose", "covariable")
  response <- model.response(data_base)

  if (!is.matrix(response))
    stop("The left side of the formula should be a matrix:
         number of DLT and number of patients for each dose!\n")

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

  if (covariable < min_cov | covariable > max_cov)
    stop("'covariable' in formula has to be between 'min_cov' and 'max_cov'.")

  if (next_patient_cov < min_cov | next_patient_cov > max_cov)
    stop("'covariable' in formula has to be between 'min_cov' and 'max_cov'.")

  if (is.null(first_dose) | is.null(last_dose)) {

    if (type == "continuous")
      stop("'first_dose' and the 'last_dose' should be informed since
           type = 'continuous'.")

    if (type == "discrete") {
      first_dose_np <- dose_set[1]
      last_dose_np <- dose_set[length(dose_set)]
      warning("'first_dose' and the 'last_dose' were defined as the first and last element of 'dose_set', respectively..")
    }

  } else {

    if (is.function(first_dose)){
      first_dose_np <- first_dose(next_patient_cov)

    } else {
      first_dose_np <- first_dose
    }

    if (is.function(last_dose)) {
      last_dose_np <- upper_dose(next_patient_cov)

    } else {
      last_dose_np <- last_dose
    }

    if (first_dose_np > last_dose_np)
      stop("'first_dose' should be smaller than the 'last_dose'.")
  }

  if (is.null(min_dose) | is.null(max_dose)) {
    min_dose_np <- first_dose_np
    min_dose_sp <- first_dose_np
    max_dose_np <- last_dose_np
    max_dose_sp <- last_dose_np

    if (type == "discrete")
      if (rounding == "down")
        max_dose_sp <- last_dose_np + 1
  } else {

    if (is.function(min_dose)){
      min_dose_sp <- min_dose(covariable)
      min_dose_np <- min_dose(next_patient_cov)

    } else {
      min_dose_sp <- min_dose
      min_dose_np <- min_dose
    }

    if (is.function(max_dose)) {
      max_dose_sp <- max_dose(covariable)
      max_dose_np <- max_dose(next_patient_cov)

    } else {
      max_dose_sp <- max_dose
      max_dose_np <- max_dose
    }

    if (any(min_dose_sp > max_dose_sp))
      stop("'min_dose' should be smaller than the 'max_dose'.")

    if (min_dose_np > max_dose_np)
      stop("'min_dose' should be smaller than the 'max_dose'.")
  }

  my_data <- list(response = response, design_matrix = design_matrix,
                  theta = theta, alpha = alpha,
                  first_dose_np = first_dose_np, last_dose_np = last_dose_np,
                  min_dose_np = min_dose_np, max_dose_np = max_dose_np,
                  min_dose_sp = min_dose_sp, max_dose_sp = max_dose_sp,
                  dose_set = dose_set,
                  mtd_prior = mtd_prior, rho_prior = rho_prior,
                  next_patient_cov = next_patient_cov,
                  max_cov = max_cov, min_cov = min_cov,
                  type = type, rounding = rounding)
  class(my_data) <- "continuous"

  out <- qmtd_jags(my_data, n_adapt, burn_in, n_mcmc, n_thin, n_chains)

  out$n <- nrow(response)
  out$alpha <- alpha
  out$theta <- theta
  out$max_dose <- max_dose_np
  out$min_dose <- min_dose_np

  class(out) <- "ewoc"
  return(out)
}

#'@export
ewoc_jags.d1continuous <- function(data, n_adapt, burn_in,
                                   n_mcmc, n_thin, n_chains) {

  theta <- data$theta
  min_dose <- data$min_dose_np
  max_dose <- data$max_dose_np
  lb <- - min_dose/(max_dose - min_dose)
  rho_prior <- data$rho_prior
  mtd_prior <- data$mtd_prior
  design_matrix <- data$design_matrix
  min_cov <- data$min_cov
  max_cov <- data$max_cov
  dlt <- data$response[, 1]
  npatients <- data$response[, 2]

  # JAGS model function
  jfun <- function() {

    for(i in 1:nobs) {
      dlt[i] ~ dbin(p[i], 1)
      p[i] <- 1/(1 + exp(-lp[i]))
      lp[i] <- inprod(design_matrix[i, ], beta)
    }

    beta[1] <- logit(rho[1]) - min_cov*(logit(rho[2]) - logit(rho[1]))/(max_cov - min_cov)
    beta[2] <- (logit(theta) - logit(rho[2]))/gamma
    beta[3] <- (logit(rho[2]) - logit(rho[1]))/(max_cov - min_cov)

    rho[1] <- theta*r[1]
    r[1] ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])
    rho[2] <- theta*r[2]
    r[2] ~ dbeta(rho_prior[2, 1], rho_prior[2, 2])
    gamma <- g[1]
    g[1] ~ dbeta(mtd_prior[1, 1], mtd_prior[1, 2])
  }


  tc1 <- textConnection("jmod", "w")
  write.model(jfun, tc1)
  close(tc1)

  data_base <- list('dlt' = dlt, 'design_matrix' = design_matrix,
                    'theta' = theta, 'nobs' = length(dlt),
                    'mtd_prior' = mtd_prior, 'rho_prior' = rho_prior,
                    'min_cov' = min_cov, 'max_cov' = max_cov)

  inits <- function() {
    out <- list(r = rbeta(nrow(rho_prior), rho_prior[, 1], rho_prior[, 2]),
                g = rbeta(nrow(mtd_prior), mtd_prior[, 1], mtd_prior[, 2]))
    return(out)
  }

  # Calling JAGS
  tc2 <- textConnection(jmod)
  j <- jags.model(tc2,
                  data = data_base,
                  inits = inits(),
                  n.chains = n_chains,
                  n.adapt = n_adapt)
  close(tc2)
  update(j, burn_in)
  sample <- coda.samples(j, variable.names = c("gamma", "rho"),
                         n.iter = n_mcmc, thin = n_thin,
                         n.chains = n_chains)[[1]]
  gamma <- sample[, 1]
  rho <- sample[, 2:3]


  out <- list(gamma = gamma, rho = rho)

  return(out)
}


response_continuous <- function(rho, dose, cov, min_cov, max_cov) {

  beta <- rep(NA, 3)
  beta[1] <- logit(rho[1]) - min_cov*(logit(rho[2]) - logit(rho[1]))/(max_cov - min_cov)
  beta[2] <- logit(rho[3]) - logit(rho[1])
  beta[3] <- (logit(rho[2]) - logit(rho[1]))/(max_cov - min_cov)

  design_matrix <- cbind(1, dose, cov)
  eta <- design_matrix%*%beta
  p <- plogis(eta)

  out <- rbinom(1, 1, p)

  return(out)
}

continuous_cov <- function(n, groups, prob_groups, cov_limits) {

  index <- rmultinom(n, 1, prob = prob_groups)
  index <- which(index == 1, arr.ind = TRUE)[, 1]
  covariate <- apply(cov_limits[index, ], 1,
                     function(x) out <- runif(1, x[1], x[2])) #
  group <- groups[index]

  out <- list(covariate = covariate, group = group)

  return(out)
}


