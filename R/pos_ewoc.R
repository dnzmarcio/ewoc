#'Escalation With Overdose Control
#'
#'Finding the next dose for a phase I clinical trial based on Escalation
#'with Overdose Control (EWOC) design using the parametrization for time
#'to event response with the proportional odds assumption and single agent.
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
#'@param mtd_prior a matrix 1x2 of hyperparameters for the Beta prior
#'distribution associated with the parameter MTD.
#'@param rho_prior a matrix 1x2 of hyperparameters for the Beta prior
#'distribution associated with the parameter rho.
#'@param shape_prior a matrix 1x2 of hyperparameters for the Gamma prior
#'distribution associated with the shape parameter r for the Weibull
#'distribution.
#'It is only necessary if distribution = 'weibull'.
#'@param type a character describing the type of the Maximum Tolerable Dose
#'(MTD) variable. It can be 'discrete' or 'continuous'.
#'@param min_dose a numerical value defining the lower bound of the support of
#'the MTD.
#'@param max_dose a numerical value defining the upper bound of the support of
#'the MTD.
#'@param first_dose a numerical value for the first allowable dose in the trial.
#'It is only necessary if type = 'continuous'.
#'@param last_dose a numerical value for the last allowable dose in the trial.
#'It is only necessary if type = 'continuous'.
#'@param dose_set a numerical vector of allowable doses in the trial. It is only
#'necessary if type = 'discrete'.
#'@param max_increment a numerical value indicating the maximum increment from the current dose to the next dose.
#'@param distribution a character establishing the distribution for the time of
#'events. It can be 'exponential' or 'weibull'.
#'@param rounding a character indicating how to round a continuous dose to the
#'one of elements of the dose set. It can be 'nearest' or 'down'.
#'It is only necessary if type = 'discrete'.
#'@param n_adapt the number of iterations for adaptation.
#'See \code{\link[rjags]{adapt}} for details.
#'@param burn_in the number of iterations before to start monitoring.
#'@param n_mcmc the number of iterations to monitor.
#'@param n_thin thinning interval for monitors.
#'@param n_chains the number of parallel chains for the model.
#'
#'@return \code{next_dose} the next recommend dose.
#'@return \code{mtd} the posterior MTD distribution.
#'@return \code{rho} the posterior rho_0 distribution.
#'@return \code{sample} a list of the MCMC chains distribution.
#'@return \code{trial} a list of the trial conditions.
#'
#'@references Tighiouart M, Liu Y, Rogatko A. Escalation with overdose control using time to toxicity for cancer phase I clinical trials. PloS one. 2014 Mar 24;9(3):e93070.
#'
#'@examples
#'time <- 9
#'status <- 0
#'dose <- 30
#'
#'test <- ewoc_d1pos(cbind(time, status) ~ dose, type = 'continuous',
#'                   theta = 0.33, alpha = 0.25, tau = 10,
#'                   min_dose = 30, max_dose = 50,
#'                   rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                   mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                   distribution = 'exponential',
#'                   rounding = 'nearest')
#'summary(test)
#'plot(test)
#'
#'@import stats
#'
#'@export
ewoc_d1pos <- function(formula, theta, alpha, tau,
                      type = c('continuous', 'discrete'),
                      rho_prior, mtd_prior, shape_prior = NULL,
                      min_dose, max_dose,
                      first_dose = NULL, last_dose = NULL,
                      dose_set = NULL, max_increment = NULL,
                      distribution = c('exponential', 'weibull'),
                      rounding = c('down', 'nearest'),
                      n_adapt = 5000, burn_in = 1000,
                      n_mcmc = 1000, n_thin = 1, n_chains = 1) {

  formula <- Formula::Formula(formula)
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

    if (is.null(max_increment))
      max_increment <- max(diff(dose_set))
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
                  dose_set = dose_set,
                  max_increment = max_increment, current_dose = current_dose,
                  rho_prior = rho_prior, mtd_prior= mtd_prior,
                  shape_prior = shape_prior,
                  distribution = distribution, tau = tau,
                  type = type, rounding = rounding)
  class(my_data) <- c("ewoc_d1pos", "d1pos")

  my_data$mcmc <- jags(my_data, n_adapt, burn_in, n_mcmc, n_thin, n_chains)
  out <- next_dose(my_data)

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

  class(out) <- c("ewoc_d1pos", "d1pos")

  return(out)
}

#'@importFrom rjags jags.model coda.samples
jags.d1pos <- function(data, n_adapt, burn_in,
                      n_mcmc, n_thin, n_chains) {

  # JAGS model function

  if (data$distribution == "weibull") {
    jfun <- "model {

    C <- 10000

    for (i in 1:nobs){
      zeros[i] ~ dpois(phi[i])

      phi[i] <- -l[i] + C

      l[i] <- status[i]*(log(h0[i] + 10^(-3)) - log(1 + s0[i]*(link[i] - 1)) + 10^(-3)) +
      log(link[i] + 10^(-3)) + log(s0[i] + 10^(-3)) -
      log(1 + s0[i]*(link[i] - 1) + 10^(-3))

      h0[i] <- exp(log(shape) + shape*log(lambda) + (shape - 1)*log(time_cens[i]))
      s0[i] <- exp(-(lambda*time_cens[i])^shape)
      link[i] <- exp(-beta*design_matrix[i, 2])
    }

      lambda <- (1/tau)*(-log(1-rho + 10^(-3)))^(1/shape)
      beta <- (log((1-rho)*theta  + 10^(-3)) -
      log((1-theta)*rho + 10^(-3)))*exp(-log(gamma + 10^(-3)))

      rho <- theta*r
      r ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])
      gamma ~ dbeta(mtd_prior[1, 1], mtd_prior[1, 2])
      shape ~ dgamma(shape_prior[1, 1], shape_prior[1, 2])
    }"

    inits <- function() {
      out <- list(r = rbeta(nrow(data$rho_prior),
                            data$rho_prior[, 1], data$rho_prior[, 2]),
                  gamma = rbeta(nrow(data$mtd_prior),
                            data$mtd_prior[, 1], data$mtd_prior[, 2]),
                  shape = rgamma(nrow(data$shape_prior),
                             data$shape_prior[, 1], data$shape_prior[, 2]))
      return(out)
    }

    data_base <- list('time_cens' = data$response[, 1],
                      'zeros' = rep(0, nrow(data$design_matrix)),
                      'status' = status, 'tau' = data$tau,
                      'design_matrix' = data$design_matrix,
                      'theta' = data$theta,
                      'nobs' = nrow(data$design_matrix),
                      'rho_prior' = data$rho_prior,
                      'mtd_prior' = data$mtd_prior,
                      'shape_prior' = data$shape_prior)
    } else {
      jfun <- "model {

    C <- 10000

    for (i in 1:nobs){
      zeros[i] ~ dpois(phi[i])

      phi[i] <- -l[i] + C

      l[i] <- status[i]*(log(h0[i] + 10^(-3)) -
                         log(1 + s0[i]*(link[i] - 1) + 10^(-3))
                         ) +
              log(link[i] + 10^(-3)) + log(s0[i] + 10^(-3)) -
              log(1 + s0[i]*(link[i] - 1) + 10^(-3))

      h0[i] <- lambda
      s0[i] <- exp(-lambda*time_cens[i])
      link[i] <- exp(-beta*design_matrix[i, 2])
    }

      lambda <- (1/tau)*(-log(1-rho + 10^(-3)))
      beta <- (log((1-rho)*theta + 10^(-3))-log((1-theta)*rho) +
           10^(-3))*exp(-log(gamma + 10^(-3)))

      rho <- theta*r
      gamma <- g

      r ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])
      g ~ dbeta(mtd_prior[1, 1], mtd_prior[1, 2])
    }"

      inits <- function() {
        out <- list(r = rbeta(nrow(data$rho_prior),
                              data$rho_prior[, 1], data$rho_prior[, 2]),
                    gamma = rbeta(nrow(data$mtd_prior),
                                  data$mtd_prior[, 1], data$mtd_prior[, 2]))
        return(out)
      }

      data_base <- list('time_cens' = data$response[, 1],
                        'zeros' = rep(0, nrow(data$design_matrix)),
                        'status' = data$response[, 2], 'tau' = data$tau,
                        'design_matrix' = data$design_matrix,
                        'theta' = data$theta,
                        'nobs' = nrow(data$design_matrix),
                        'rho_prior' = data$rho_prior,
                        'mtd_prior' = data$mtd_prior)
      }

  initial <- inits()
  # Calling JAGS
  j <- jags.model(textConnection(jfun),
                  data = data_base,
                  inits = initial,
                  n.chains = n_chains,
                  n.adapt = n_adapt)
  update(j, burn_in)

  if (data$distribution == "weibull"){
    sample <- coda.samples(j,
                           variable.names =
                             c("beta", "gamma", "lambda", "rho", "shape"),
                           n.iter = n_mcmc, thin = n_thin,
                           n.chains = n_chains)

    beta <- sample[[1]][, 1]
    gamma <- sample[[1]][, 2]
    lambda <- sample[[1]][, 3]
    rho <- sample[[1]][, 4]
    shape <- sample[[1]][, 5]

    out <- list(beta = beta, gamma = gamma, lambda = lambda,
                rho = rho, shape = shape,
                sample = sample)
  } else {
    sample <- coda.samples(j, variable.names =  c("beta", "gamma", "lambda", "rho"),
                           n.iter = n_mcmc, thin = n_thin,
                           n.chains = n_chains)

    beta <- sample[[1]][, 1]
    gamma <- sample[[1]][, 2]
    lambda <- sample[[1]][, 3]
    rho <- sample[[1]][, 4]
    out <- list(beta = beta, gamma = gamma, lambda = lambda,
                rho = rho, sample = sample)
  }
  return(out)
  }
