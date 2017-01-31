#'Escalation Over with Dose Control
#'
#'Finding the next dose for a phase I clinical trial based on Escalation Over
#'Dose Control (EWOC) design considering extended parametrization for time
#'to event response and single agent.
#'
#'@param formula an object of class \code{\link[Formula]{Formula}}: a symbolic
#'description of the model to be fitted with only one regressor term
#'corresponding to the dose for the right side and a matrix as a response
#'containing time and status for the left side.
#'@param theta a numerical value defining the proportion of expected patients
#'to experience a medically unacceptable, dose-limiting toxicity (DLT) if
#'administered the MTD.
#'@param alpha a numerical value defining the probability that the dose selected
#'by EWOC is higher than the MTD.
#'@param tau a numerical value defining the period of time for a possible
#'toxicity be observed.
#'@param rho_prior a matrix 2x2 of hyperparameters for the Beta prior
#'distribution associated with each rho. Each row corresponds to a paramater.
#'@param shape_prior a vector of hyperparameters for the Gamma prior
#'distribution associated with the shape parameter r for the Weibull
#'distribution.
#'It is only necessary if distribution = 'weibull'.
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
#'@param distribution a character establishing the distribution for the time of
#'events.
#'@param rounding a character indicating how to round a continuous dose to the
#'one of elements of the dose set. It is only necessary if type = discrete'.
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
#'@return \code{rho} the posterior rho_0 and rho_1 distributions.
#'@return \code{sample} a list of the MCMC chains distribution.
#'@return \code{trial} a list of trial conditions.
#'
#'@examples
#'time <- 9
#'status <- 0
#'dose <- 30
#'
#'test <- ewoc_d1ph(cbind(time, status) ~ dose, type = 'discrete',
#'                  theta = 0.33, alpha = 0.25, tau = 10,
#'                  dose_set = seq(30, 50, 5),
#'                  rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                  mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                  distribution = 'exponential',
#'                  rounding = 'down')
#'summary(test)
#'
#'test <- ewoc_d1ph(cbind(time, status) ~ dose, type = 'discrete',
#'                  theta = 0.33, alpha = 0.25, tau = 10,
#'                  dose_set = seq(30, 50, 5),
#'                  rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                  mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                  shape_prior = matrix(1, ncol = 2, nrow = 1),
#'                  distribution = 'weibull',
#'                  rounding = 'down')
#'summary(test)
#'
#'test <- ewoc_d1ph(cbind(time, status) ~ dose, type = 'continuous',
#'                  theta = 0.33, alpha = 0.25, tau = 10,
#'                  first_dose = 30, last_dose = 50,
#'                  rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                  mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                  distribution = 'exponential')
#'summary(test)
#'plot(test)
#'
#'test <- ewoc_d1ph(cbind(time, status) ~ dose, type = 'continuous',
#'                  theta = 0.33, alpha = 0.25, tau = 10,
#'                  first_dose = 30, last_dose = 50,
#'                  rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                  mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                  shape_prior = matrix(1, ncol = 2, nrow = 1),
#'                  distribution = 'weibull')
#'summary(test)
#'plot(test)
#'
#'@references Tighiouart M, Liu Y, Rogatko A. Escalation with overdose control using time to toxicity for cancer phase I clinical trials. PloS one. 2014 Mar 24;9(3):e93070.
#'
#'@export
ewoc_d1ph <- function(formula, theta, alpha, tau,
                      type = c('continuous', 'discrete'),
                      rho_prior, mtd_prior, shape_prior = NULL,
                      first_dose = NULL, last_dose = NULL,
                      dose_set = NULL,
                      min_dose = NULL, max_dose = NULL,
                      distribution = c('exponential', 'weibull'),
                      rounding = c("down", "nearest"),
                      n_adapt = 5000, burn_in = 1000,
                      n_mcmc = 1000, n_thin = 1, n_chains = 1) {

  formula <- Formula(formula)
  if (class(formula)[2] != "formula")
    stop("Invalid formula! \n")

  data_base <- model.frame(formula, na.action = na.exclude,
                           drop.unused.levels = FALSE)

  dose_matrix <- model.matrix(formula, data_base, rhs = 1)

  if (length(formula)[2] == 1){
    covariable_matrix <- NULL
    design_matrix <- dose_matrix
    colnames(design_matrix) <- c("intercept", "dose")
  } else {
    stop("This design cannot accommodate a covariable.")
  }

  response <- model.response(data_base)

  if (!is.matrix(response))
    stop("The left side of the formula should be a matrix: time and status!\n")

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

  if (nrow(rho_prior) != 1 | ncol(rho_prior) != 2)
    stop("'rho_prior' should be a matrix with 1 column and 2 rows.")

  if (nrow(mtd_prior) != 1 | ncol(mtd_prior) != 2)
    stop("'mtd_prior' should be a matrix with 1 column and 2 rows.")

  if (distribution == 'weibull')
    if (is.null(shape_prior)) {
      stop("'shape_prior' should be informed if 'distribution' = 'weibull'")
    } else {
      if (!(nrow(shape_prior) == 1 & ncol(shape_prior) == 2))
        stop("'shape_prior' should be a matrix with 2 columns and 1 row.")
    }

  limits <- limits_d1nocov(first_dose = first_dose, last_dose = last_dose,
                           min_dose = min_dose, max_dose = max_dose,
                           type = type, rounding = rounding,
                           dose_set = dose_set)

  design_matrix[, 2] <-
    standard_dose(dose = design_matrix[, 2],
                  min_dose = limits$min_dose(),
                  max_dose = limits$max_dose())

  my_data <- list(response = response, design_matrix = design_matrix,
                  theta = theta, alpha = alpha,
                  limits = limits,
                  dose_set = dose_set,
                  rho_prior = rho_prior, mtd_prior= mtd_prior,
                  shape_prior = shape_prior,
                  distribution = distribution, tau = tau,
                  type = type, rounding = rounding)
  class(my_data) <- "d1ph"

  out <- qmtd_jags(my_data, n_adapt, burn_in, n_mcmc, n_thin, n_chains)

  trial <- list(response = response, design_matrix = design_matrix,
                theta = theta, alpha = alpha,
                first_dose = limits$first_dose, last_dose = limits$last_dose,
                min_dose = limits$min_dose, max_dose = limits$max_dose,
                dose_set = dose_set,
                rho_prior = rho_prior, mtd_prior = mtd_prior,
                shape_prior = shape_prior,
                distribution = distribution, tau = tau,
                type = type, rounding = rounding,
                n_adapt = n_adapt, burn_in = burn_in, n_mcmc = n_mcmc,
                n_thin = n_thin, n_chains = n_chains)
  out$trial <- trial

  class(out) <- c("ewoc_d1ph", "d1ph")

  return(out)
}

#'@export
ewoc_jags.d1ph <- function(data, n_adapt, burn_in,
                         n_mcmc, n_thin, n_chains) {

  theta <- data$theta
  rho_prior <- data$rho_prior
  mtd_prior <- data$mtd_prior
  shape_prior <- data$shape_prior
  design_matrix <- data$design_matrix
  dose <- data$design_matrix
  time_cens <- data$response[, 1]
  status <- data$response[, 2]
  distribution <- data$distribution
  tau <- data$tau

  time_mod <- time_cens
  time_mod[status == 0] <- NA
  censored <- as.numeric(!status)

  # JAGS model function

  if (distribution == "weibull") {
    jfun <- function() {

      for(i in 1:nobs) {
        censored[i] ~ dinterval(time_mod[i], time_cens[i])
        time_mod[i] ~ dweib(shape, rate[i])
        rate[i] <- exp(inprod(design_matrix[i, ], beta))
      }

      beta[1] <- log(-log(1 - rho[1])/(tau^shape))
      beta[2] <- log(log(1 - theta)/log(1 - rho[1]))/gamma

      rho[1] <- theta*r
      r ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])
      gamma ~ dbeta(mtd_prior[1, 1], mtd_prior[1, 2])
      shape ~ dgamma(shape_prior[1, 1], shape_prior[1, 2])
    }

    inits <- function() {
      time_init <- rep(NA, length(time_mod))
      time_init[which(!status)] <- time_cens[which(!status)] + 1

      out <- list(r = rbeta(nrow(rho_prior),
                            rho_prior[, 1], rho_prior[, 2]),
                  gamma = rbeta(nrow(rho_prior),
                                mtd_prior[, 1], mtd_prior[, 2]),
                  shape = rgamma(1, shape_prior[1], shape_prior[2]),
                  time_mod = time_init)
      return(out)
    }

    data_base <- list('time_mod' = time_mod, 'time_cens' = time_cens,
                      'censored' = censored, 'tau' = tau,
                      'design_matrix' = design_matrix, 'theta' = theta,
                      'nobs' = length(time_cens[!is.na(time_cens)]),
                      'rho_prior' = rho_prior,
                      'mtd_prior' = mtd_prior,'shape_prior' = shape_prior)
  } else {
    jfun <- function() {

      for(i in 1:nobs) {
        censored[i] ~ dinterval(time_mod[i], time_cens[i])
        time_mod[i] ~ dexp(rate[i])
        rate[i] <- exp(inprod(design_matrix[i, ], beta))
      }

      beta[1] <- log(-log(1 - rho[1])/(tau^1))
      beta[2] <- log(log(1 - theta)/log(1 - rho[1]))/gamma

      rho[1] <- theta*r
      r ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])
      gamma ~ dbeta(mtd_prior[1, 1], mtd_prior[1, 2])
    }

    inits <- function() {
      time_init <- rep(NA, length(time_mod))
      time_init[which(!status)] <- time_cens[which(!status)] + 1

      out <- list(r = rbeta(nrow(rho_prior),
                            rho_prior[, 1], rho_prior[, 2]),
                  gamma = rbeta(nrow(rho_prior),
                                mtd_prior[, 1], mtd_prior[, 2]),
                  time_mod = time_init)
      return(out)
    }

    data_base <- list('time_mod' = time_mod, 'time_cens' = time_cens,
                      'censored' = censored, 'tau' = tau,
                      'design_matrix' = design_matrix, 'theta' = theta,
                      'nobs' = length(time_cens[!is.na(time_cens)]),
                      'rho_prior' = rho_prior, 'mtd_prior' = mtd_prior)

  }

  tc1 <- textConnection("jmod", "w")
  write.model(jfun, tc1)
  close(tc1)

  # Calling JAGS
  tc2 <- textConnection(jmod)
  j <- jags.model(tc2,
                  data = data_base,
                  inits = inits(),
                  n.chains = n_chains,
                  n.adapt = n_adapt)
  close(tc2)
  update(j, burn_in)

  if (distribution == "weibull"){
    result <- coda.samples(j, variable.names =  c("gamma", "rho", "shape"),
                           n.iter = n_mcmc, thin = n_thin, n.chains = n_chains)[[1]]

    gamma <- result[, 1]
    rho <- result[, 2]
    shape <- result[, 3]

    out <- list(gamma = gamma, rho = rho, shape = shape)
  } else {
    result <- coda.samples(j, variable.names =  c("gamma", "rho"),
                           n.iter = n_mcmc, thin = n_thin, n.chains = n_chains)[[1]]

    gamma <- result[, 1]
    rho <- result[, 2]

    out <- list(gamma = gamma, rho = rho)
  }
  return(out)
}



#'@export
rho_ph <- function(rho0, true_mtd, min_dose, max_dose, theta) {

  p0 <- (max_dose - min_dose)/(true_mtd - min_dose)
  p1 <- (log(1 - theta)/log(1 - rho0))^p0
  p2 <- log(1 - rho0)*p1
  out <- 1 - exp(p2)

  return(out)
}

