#'Escalation Over with Dose Control
#'
#'Finding the next dose for a phase I clinical trial based on Escalation Over
#'Dose Control (EWOC) design considering classical parametrization for binary
#'response with multinomial covariable.
#'
#'@param formula an object of class \code{\link[Formula]{Formula}}: a symbolic
#'description of the model to be fitted with two regressor parts separated by `|`
#'corresponding to the dose and covariable, respectively, for the right side and
#'a vector as a response containing number of DLT for the left side.
#'@param theta a numerical value defining the proportion of expected patients
#'to experience a medically unacceptable, dose-limiting toxicity (DLT) if
#'administered the MTD.
#'@param alpha a numerical value defining the probability that the dose selected
#'by EWOC is higher than the MTD.
#'@param order a logical value indicating if the MTD prior should be ordered.
#'@param mtd_prior a matrix of hyperparameters for the Beta prior distribution
#'associated with the MTD. Each row corresponds to a paramater.
#'@param rho_prior a matrix of hyperparameters for the Beta prior distribution
#'associated with each rho. Each row corresponds to a paramater.
#'@param next_patient_cov a character value indicating the covariable associated to
#'next patient.
#'@param levels_cov a character vector of the possible values for the ordinal covariable.
#'@param type a character describing the type of the Maximum Tolerable Dose
#'(MTD) variable.
#'@param min_dose either a numerical value or a function of the covariable
#'defining the lower bound of the support of the MTD used to standardize the doses.
#'@param max_dose either a numerical value or a function of the covariable
#'defining the upper bound of the support of the MTD used to standardize the doses.
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
#'@return \code{next_dose} the next recommend dose.
#'@return \code{mtd} the posterior MTD distribution for all levels of the covariable.
#'@return \code{rho} the posterior rho_01 distribution.
#'@return \code{gamma} the posterior standardized MTD distributions for all levels of the covariable.
#'@return \code{sample} a list of the MCMC chains distribution.
#'@return \code{trial} a list of the trial conditions.
#'
#'@references Tighiouart M, Cook-Wiens G, Rogatko A. Incorporating a patient dichotomous characteristic in cancer phase I clinical trials using escalation with overdose control. Journal of Probability and Statistics. 2012 Oct 2;2012.
#'
#'@examples
#'\dontrun{
#'library(ewoc)
#'
#'DLT <- rep(0, 2)
#'group <- c("B", "C")
#'group <- factor(group, levels = c("A", "B", "C"))
#'dose <- rep(30, 2)
#'test <- ewoc_d1multinomial(DLT ~ dose | group,
#'                           type = 'continuous',
#'                           theta = 0.33, alpha = 0.25,
#'                           min_dose = 30, max_dose = 50,
#'                           levels_cov = c("A", "B", "C"),
#'                           next_patient_cov = "A",
#'                           mtd_prior = matrix(1, nrow = 3, ncol = 2),
#'                           rho_prior = matrix(1, nrow = 1, ncol = 2))
#'summary(test)
#'plot(test)
#'}
#'
#'
#'@export
ewoc_d1multinomial <- function(formula, theta, alpha,
                               mtd_prior, rho_prior,
                               levels_cov, next_patient_cov,
                               min_dose, max_dose,
                               type = c('continuous', 'discrete'),
                               order = FALSE,
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
    stop("This design requires a multinomial covariable.")
  } else {
    nlc <- length(levels_cov)
    covariable_matrix <- model.matrix(formula, data_base, rhs = 2)
    index <- apply(covariable_matrix, 1, function(x) max(which(x != 0)))
    covariable <- levels_cov[index]
  }

  design_matrix <- cbind(dose_matrix,
                         matrix(covariable_matrix[, -1],
                                ncol = (length(levels_cov) - 1)))
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

  if(nrow(mtd_prior) != length(levels_cov) | ncol(mtd_prior) != 2)
    stop(paste0("'mtd_prior' should be a matrix with 2 columns and ", nlc, " rows."))

  if(nrow(rho_prior) != 1 | ncol(mtd_prior) != 2)
    stop(paste0("'rho_prior' should be a matrix with 2 columns and 1 row."))

  if(any(!(next_patient_cov %in% levels_cov)))
    stop(paste0("'next_patient_cov' should be choose from 'levels_cov'."))

  limits <- limits_d1cov(first_dose = first_dose, last_dose = last_dose,
                         min_dose = min_dose, max_dose = max_dose,
                         type = type, rounding = rounding,
                         dose_set = dose_set,
                         covariable = next_patient_cov)

  design_matrix[, 2] <- standard_dose(dose = design_matrix[, 2],
                                      min_dose = limits$min_dose(covariable),
                                      max_dose = limits$max_dose(covariable))

  my_data <- list(response = response, design_matrix = design_matrix,
                  theta = theta, alpha = alpha, limits = limits,
                  dose_set = dose_set, order = order,
                  rho_prior = rho_prior, mtd_prior = mtd_prior,
                  levels_cov = levels_cov,
                  type = type, rounding = rounding)
  class(my_data) <- "d1multinomial"

  out <- qmtd_jags(my_data, n_adapt, burn_in, n_mcmc, n_thin, n_chains)

  trial <- list(response = response, design_matrix = design_matrix,
                theta = theta, alpha = alpha,
                first_dose = limits$first_dose, last_dose = limits$last_dose,
                min_dose = limits$min_dose, max_dose = limits$max_dose,
                dose_set = dose_set, order = order,
                rho_prior = rho_prior, mtd_prior = mtd_prior,
                levels_cov = levels_cov, covariable = covariable,
                next_patient_cov = next_patient_cov,
                type = type, rounding = rounding,
                n_adapt = n_adapt, burn_in = burn_in, n_mcmc = n_mcmc,
                n_thin = n_thin, n_chains = n_chains)
  out$trial <- trial

  class(out) <- c("ewoc_d1multinomial", "d1multinomial")

  return(out)
}



#'@export
ewoc_jags.d1multinomial <- function(data, n_adapt, burn_in,
                                    n_mcmc, n_thin, n_chains) {

  if (!data$order){
  # JAGS model function
    jfun <- function() {

      for(i in 1:nobs) {
        dlt[i] ~ dbin(p[i], 1)
        p[i] <- ifelse(1/(1 + exp(-eta[i])) == 1, 0.99, 1/(1 + exp(-eta[i])))
        eta[i] <- inprod(design_matrix[i, ], beta)
      }

      for (i in 3:(np+1)) {
        beta[i] <- logit(theta) - logit(rho[1]) -
          gamma[i-1]*(logit(theta) - beta[1])/gamma[i-2]
      }
      beta[2] <- (logit(theta) - logit(rho[1]))/gamma[1]
      beta[1] <- logit(rho[1])

      rho[1] <- theta*r[1]
      r[1] ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])

      for (i in 1:np) {
        gamma[i] <- v[i]
        v[i] ~ dbeta(mtd_prior[i, 1], mtd_prior[i, 2])
      }
    }
  } else {
    jfun <- function() {

      for(i in 1:nobs) {
        dlt[i] ~ dbin(p[i], 1)
        p[i] <- ifelse(1/(1 + exp(-eta[i])) == 1, 0.99, 1/(1 + exp(-eta[i])))
        eta[i] <- inprod(design_matrix[i, ], beta)
      }

      for (i in 3:(np+1)) {
        beta[i] <- logit(theta) - logit(rho[1]) -
          gamma[i-1]*(logit(theta) - beta[1])/gamma[i-2]
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
  }

  tc1 <- textConnection("jmod", "w")
  R2WinBUGS::write.model(jfun, tc1)
  close(tc1)

  data_base <- list('dlt' = data$response,
                    'design_matrix' = data$design_matrix,
                    'theta' = data$theta,
                    'nobs' = length(data$response),
                    'np' = ncol(data$design_matrix) - 1,
                    'rho_prior' = data$rho_prior, 'mtd_prior' = data$mtd_prior)

  inits <- function() {
    out <- list(r = rbeta(nrow(data$rho_prior),
                          data$rho_prior[, 1], data$rho_prior[, 2]),
                v = rbeta(nrow(data$mtd_prior),
                          data$mtd_prior[, 1], data$mtd_prior[, 2]))
    return(out)
  }

  # Calling JAGS
  tc2 <- textConnection(jmod)
  j <- rjags::jags.model(tc2,
                         data = data_base,
                         inits = inits(),
                         n.chains = n_chains,
                         n.adapt = n_adapt)
  close(tc2)
  update(j, burn_in)
  sample <- rjags::coda.samples(j, variable.names =  c("gamma", "rho"),
                                n.iter = n_mcmc, thin = n_thin,
                                n.chains = n_chains)
  gamma <- sample[[1]][, 1:(ncol(data$design_matrix) - 1)]
  rho <- sample[[1]][, ncol(data$design_matrix)]

  out <- list(gamma = gamma, rho = rho, sample = sample)

  return(out)
}
