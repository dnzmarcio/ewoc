#'Escalation Over with Dose Control
#'
#'Finding the next dose for a phase I clinical trial based on Escalation Over
#'Dose Control (EWOC) design considering extended parametrization for binary
#'response with ordinal covariable.
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
#'@param mtd_prior a matrix of hyperparameters for the Beta prior distribution
#'associated with the MTD. Each row corresponds to a paramater.
#'@param rho_prior a matrix of hyperparameters for the Beta prior distribution
#'associated with each rho. Each row corresponds to a paramater.
#'@param next_patient_cov a character value indicating the covariable associated to
#'next patient.
#'@param levels_cov a character vector of the possible values for the ordinal covariable.
#'@param type a character describing the type of the Maximum Tolerable Dose
#'(MTD) variable.
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
#'@param rounding a character indicating how to round a continuous dose to the
#'one of elements of the dose set. It is only necessary if type = 'discrete'.
#'@param n_adapt the number of iterations for adaptation.
#'See \code{\link[rjags]{adapt}} for details.
#'@param burn_in the number of iterations before to start monitoring.
#'@param n_mcmc the number of iterations to monitor.
#'@param n_thin thinning interval for monitors.
#'@param n_chains the number of parallel chains for the model.
#'
#'@return \code{next_dose} the next recommend dose.
#'@return \code{hpd_dose} the 95\% HPD for the next dose.
#'@return \code{pdlt} the probability of DLT for the next dose.
#'@return \code{hpd_pdlt} the 95\% HPD for the probability of DLT for the next dose.
#'@return \code{mtd} the posterior MTD distribution considering the next patient covariable.
#'@return \code{rho} the posterior rho_01 distribution.
#'@return \code{gamma} the posterior standardized MTD distributions for all levels of the covariable.
#'@return \code{sample} a list of the MCMC chains distribution.
#'@return \code{trial} a list of trial conditions.
#'
#'@examples
#'DLT <- rep(0, 1)
#'npatients <- rep(1, 1)
#'group <- "B"
#'group <- factor(group, levels = c("A", "B", "C"))
#'dose <- rep(30, 1)
#'test <- ewoc1d_ordinal(cbind(DLT, npatients) ~ dose | group, type = 'continuous',
#'                       theta = 0.33, alpha = 0.25,
#'                       first_dose = 30, last_dose = 50,
#'                       min_dose = 30, max_dose = 50,
#'                       levels_cov = c("A", "B", "C"),
#'                       mtd_prior = matrix(1, nrow = 3, ncol = 2),
#'                       rho_prior = matrix(1, nrow = 1, ncol = 2))
#'summary(test)
#'plot(test)
#'
#'test <- ewoc1d_ordinal(cbind(DLT, npatients) ~ dose | group, type = 'discrete',
#'                       theta = 0.33, alpha = 0.25,
#'                       dose_set = seq(30, 50, 5),
#'                       first_dose = 30, last_dose = 50,
#'                       levels_cov = c("A", "B", "C"),
#'                       mtd_prior = matrix(1, nrow = 3, ncol = 2),
#'                       rho_prior = matrix(1, nrow = 1, ncol = 2))
#'summary(test)
#'plot(test)
#'
#'@references Tighiouart M, Cook-Wiens G, Rogatko A. Incorporating a patient dichotomous characteristic in cancer phase I clinical trials using escalation with overdose control. Journal of Probability and Statistics. 2012 Oct 2;2012.
#'
#'@export
ewoc_d1ordinal <- function(formula, theta, alpha,
                           rho_prior, mtd_prior,
                           levels_cov,
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
    stop("This design requires a ordinal covariable.")
  } else {
    nlc <- length(levels_cov)
    covariable_matrix <- model.matrix(formula, data_base, rhs = 2)
    covariable <- levels_cov[rowSums(covariable_matrix)]
  }

  design_matrix <- cbind(dose_matrix,
                         matrix(covariable_matrix[, -1], ncol = (nlc - 1)))
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

  if(nrow(mtd_prior) != nlg | ncol(mtd_prior) != 2)
    stop(paste0("'mtd_prior' should be a matrix with 2 columns and ", nlg, " rows."))

  if(nrow(rho_prior) != 1 | ncol(mtd_prior) != 2)
    stop(paste0("'rho_prior' should be a matrix with 2 columns and 1 row."))

  limits <- limits_d1cov(first_dose = first_dose, last_dose = last_dose,
                         min_dose = min_dose, max_dose = max_dose,
                         type = type, rounding = rounding,
                         dose_set = dose_set,
                         covariable = covariable)

  design_matrix[, 2] <-
    standard_dose(dose = design_matrix[, 2],
                  min_dose = limits$min_dose(covariable),
                  max_dose = limits$max_dose(covariable))

  my_data <- list(response = response, design_matrix = design_matrix,
                  theta = theta, alpha = alpha, limits = limits,
                  dose_set = dose_set,
                  rho_prior = rho_prior, mtd_prior = mtd_prior,
                  levels_cov = levels_cov, covariable = covariable,
                  type = type, rounding = rounding)
  class(my_data) <- "d1ordinal"

  out <- qmtd_jags(my_data, n_adapt, burn_in, n_mcmc, n_thin, n_chains)

  trial <- list(response = response, design_matrix = design_matrix,
                theta = theta, alpha = alpha,
                first_dose = limits$first_dose, last_dose = limits$last_dose,
                min_dose = limits$min_dose, max_dose = limits$max_dose,
                dose_set = dose_set,
                rho_prior = rho_prior, mtd_prior = mtd_prior,
                levels_cov = levels_cov, covariable = covariable_matrix,
                type = type, rounding = rounding,
                n_adapt = n_adapt, burn_in = burn_in, n_mcmc = n_mcmc,
                n_thin = n_thin, n_chains = n_chains)
  out$trial <- trial

  class(out) <- c("ewoc_d1ordinal", "d1ordinal")
  return(out)
}



#'@export
ewoc_jags.d1ordinal <- function(data, n_adapt, burn_in,
                              n_mcmc, n_thin, n_chains) {

  theta <- data$theta
  min_dose <- data$limits$min_dose
  max_dose <- data$limits$max_dose
  rho_prior <- data$rho_prior
  mtd_prior <- data$mtd_prior
  design_matrix <- data$design_matrix
  np <- ncol(design_matrix) - 1
  dose <- data$design_matrix
  dlt <- data$response[, 1]
  npatients <- data$response[, 2]

  # JAGS model function
  jfun <- function() {

    for(i in 1:nobs) {
      dlt[i] ~ dbin(p[i], 1)
      p[i] <- ifelse(1/(1 + exp(-eta[i])) == 1, 0.99, 1/(1 + exp(-eta[i])))
      eta[i] <- inprod(design_matrix[i, ], beta)
    }

    for (i in 3:(np+1)) {
      beta[i] <- logit(theta) - logit(rho[1]) - gamma[i-1]*(logit(theta) -
                                                              beta[1])/gamma[i-2]
    }
    beta[2] <- (logit(theta) - logit(rho[1]))/gamma[1]
    beta[1] <- logit(rho[1])

    rho[1] <- theta*r[1]
    r[1] ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])

    gamma[1] <- v[1]
    v[1] ~ dbeta(mtd_prior[1, 1], mtd_prior[1, 2])

    for (i in 2:np) {
      gamma[i] <- gamma[i-1]*v[i]
      v[i] ~ dbeta(mtd_prior[i, 1], mtd_prior[i, 2])
    }
  }

  tc1 <- textConnection("jmod", "w")
  write.model(jfun, tc1)
  close(tc1)

  data_base <- list('dlt' = dlt, 'design_matrix' = design_matrix, 'theta' = theta,
                    'nobs' = length(dlt), 'np' = np,
                    'rho_prior' = rho_prior, 'mtd_prior' = mtd_prior)

  inits <- function() {
    out <- list(r = rbeta(nrow(rho_prior),
                          rho_prior[, 1], rho_prior[, 2]),
                v = rbeta(nrow(mtd_prior),
                            mtd_prior[, 1], mtd_prior[, 2]))
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
  sample <- coda.samples(j, variable.names =  c("gamma", "rho"),
                         n.iter = n_mcmc, thin = n_thin,
                         n.chains = n_chains)[[1]]
  gamma <- sample[, 1:np]
  rho <- sample[, (np + 1)]

  out <- list(gamma = gamma, rho = rho)

  return(out)
}

#'@export
response_ordinal <- function(n, dose, rho, cov) {

  cov <- matrix(cov, nrow = 1)
  np <- 2 + ncol(cov)

  beta <- rep(NA, np)
  beta[1] <- logit(rho[1])
  beta[np] <- logit(rho[np]) - logit(rho[1])
  for (i in 2:(np-1))
    beta[i] <- logit(rho[i]) - beta[1]

  design_matrix <- cbind(1, cov, dose)
  lp <- design_matrix%*%beta
  p <- plogis(lp)

  out <- rbinom(n, 1, p)

  return(out)
}

#'@export
rho_ordinal <- function(rho11, true_mtd, min_dose, max_dose, theta) {

  aux <- function(rho11, true_mtd, theta) {
    rho0 <- rep(NA, length(true_mtd))
    rho0[1] <- plogis((logit(theta) - true_mtd[1]*
                         logit(rho11))/(1 - true_mtd[1]))
    for(i in 2:length(true_mtd)) {
      rho0[i] <- plogis(logit(theta) - true_mtd[i]*
                          (logit(rho11) - logit(rho0[1])))
    }
    return(rho0)
  }

  true_mtd <- (true_mtd - min_dose)/(max_dose - min_dose)

  rho0_Vec <- Vectorize(aux, vectorize.args = "rho11")
  out <-  rho0_Vec(rho11, true_mtd, theta)

  return(out)
}
