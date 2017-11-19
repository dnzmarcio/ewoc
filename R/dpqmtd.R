#'@importFrom rjags jags.model coda.samples
jags.d1classic <- function(data, n_adapt, burn_in,
                           n_mcmc, n_thin, n_chains) {

  # JAGS model function
  jfun <- "model {

  for(i in 1:nobs) {
  dlt[i] ~ dbin(p[i], 1)
  p[i] <- ifelse(1/(1 + exp(-lp[i])) == 1, 0.99, 1/(1 + exp(-lp[i])))
  lp[i] <- inprod(design_matrix[i, ], beta)
  }

  beta[1] <-logit(rho)
  beta[2] <- (logit(theta) - logit(rho))/gamma

  rho <- theta*v[1]
  v[1] ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])

  gamma <- v[2]
  v[2] ~ dbeta(mtd_prior[1, 1], mtd_prior[1, 2])
}"

  data_base <- list('dlt' = data$response,
                    'design_matrix' = data$design_matrix,
                    'theta' = data$theta,
                    'nobs' = length(data$response),
                    'rho_prior' = data$rho_prior,
                    'mtd_prior' = data$mtd_prior)

  inits <- function() {
    v <- rep(NA, 2)
    v[1] <- rbeta(1, data$rho_prior[1], data$rho_prior[2])
    v[2] <- rbeta(1, data$mtd_prior[1], data$mtd_prior[2])
    out <- list(v = v)
    return(v)
  }

  # Calling JAGS
  j <- jags.model(textConnection(jfun),
                  data = data_base,
                  inits = list(v = inits()),
                  n.chains = n_chains,
                  n.adapt = n_adapt)
  update(j, burn_in)
  sample <- coda.samples(j, variable.names = c("beta", "gamma", "rho"),
                         n.iter = n_mcmc, thin = n_thin,
                         n.chains = n_chains)

  beta <- sample[[1]][, 1:2]
  gamma <- sample[[1]][, 3]
  rho <- sample[[1]][, 4]
  out <- list(beta = beta, gamma = gamma, rho = rho, sample = sample)

  return(out)
}
#'@importFrom rjags jags.model coda.samples
jags.d1extended <- function(data, n_adapt, burn_in,
                            n_mcmc, n_thin, n_chains) {

  min_dose <- data$limits$min_dose
  max_dose <- data$limits$max_dose
  lb <- - min_dose/(max_dose - min_dose)

  # JAGS model function
  jfun <- "model {

  for(i in 1:nobs) {
  dlt[i] ~ dbin(p[i], 1)
  p[i] <- ifelse(1/(1 + exp(-lp[i])) == 1, 0.99, 1/(1 + exp(-lp[i])))
  lp[i] <- inprod(design_matrix[i, ], beta)
  }

  beta[1] <- logit(rho[1])
  beta[2] <- logit(rho[2]) - logit(rho[1])

  rho[1] <- min*v[1]
  min <- min(rho[2], limit)
  v[1] ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])

  limit <- plogis(numerator/denominator, 0, 1)
  numerator <- logit(theta) - lb*logit(rho[2])
  denominator <- 1 - lb

  rho[2] <- v[2]
  v[2] ~ dbeta(rho_prior[2, 1], rho_prior[2, 2])
}"

  data_base <- list('dlt' = data$response,
                    'design_matrix' = data$design_matrix,
                    'nobs' = length(data$response),
                    'rho_prior' = data$rho_prior,
                    'theta' = data$theta,
                    'lb' = lb)

  inits <- function() {
    v <- rep(NA, 2)
    v <- rbeta(2, data$rho_prior[, 1], data$rho_prior[, 2])
    out <- list(v = v)
    return(v)
  }

  # Calling JAGS
  j <- jags.model(textConnection(jfun),
                  data = data_base,
                  inits = list(v = inits()),
                  n.chains = n_chains,
                  n.adapt = n_adapt)
  update(j, burn_in)
  sample <- coda.samples(j, variable.names = c("beta", "rho"),
                         n.iter = n_mcmc, thin = n_thin,
                         n.chains = n_chains)

  beta <- sample[[1]][, 1:2]
  rho <- sample[[1]][, 3:4]

  out <- list(beta = beta, rho = rho, sample = sample)

  return(out)
}

#'@importFrom rjags jags.model coda.samples
jags.d1ph <- function(data, n_adapt, burn_in,
                      n_mcmc, n_thin, n_chains) {

  time_cens <- data$response[, 1]
  status <- data$response[, 2]
  time_mod <- time_cens
  time_mod[status == 0] <- NA
  censored <- as.numeric(!status)

  # JAGS model function

  if (data$distribution == "weibull") {
    jfun <- "model {

    for(i in 1:nobs) {
    censored[i] ~ dinterval(time_mod[i], time_cens[i])
    time_mod[i] ~ dweib(shape, rate[i])
    rate[i] <- exp(inprod(design_matrix[i, ], beta))
    }

    beta[1] <- log(-log(1 - rho)) - shape*log(tau)
    beta[2] <- (log(-log(1 - theta)) -
    log(-log(1 - rho)))*
    exp(-log(gamma))

    rho <- theta*r
    gamma <- g + 10^(-2)
    shape <- s + 10^(-2)
    r ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])
    g ~ dbeta(mtd_prior[1, 1], mtd_prior[1, 2])
    s ~ dgamma(shape_prior[1, 1], shape_prior[1, 2])
  }"

    inits <- function() {
      time_init <- rep(NA, length(time_mod))
      time_init[which(!status)] <- time_cens[which(!status)] + 1

      out <- list(r = rbeta(nrow(data$rho_prior),
                            data$rho_prior[, 1], data$rho_prior[, 2]),
                  g = rbeta(nrow(data$mtd_prior),
                            data$mtd_prior[, 1], data$mtd_prior[, 2]),
                  s = rgamma(nrow(data$shape_prior),
                             data$shape_prior[, 1], data$shape_prior[, 2]),
                  time_mod = time_init)
      return(out)
    }

    data_base <- list('time_mod' = time_mod, 'time_cens' = time_cens,
                      'censored' = censored, 'tau' = data$tau,
                      'design_matrix' = data$design_matrix,
                      'theta' = data$theta,
                      'nobs' = length(time_cens[!is.na(time_cens)]),
                      'rho_prior' = data$rho_prior,
                      'mtd_prior' = data$mtd_prior,
                      'shape_prior' = data$shape_prior)
  } else {
    jfun <- "model {

  for(i in 1:nobs) {
  censored[i] ~ dinterval(time_mod[i], time_cens[i])
  time_mod[i] ~ dexp(rate[i])
  rate[i] <- exp(inprod(design_matrix[i, ], beta) + 10^(-3))
  }

  beta[1] <- log(-log(1 - rho[1])) - log(tau)
  beta[2] <- (log(-log(1 - theta)) -
  log(-log(1 - rho[1])))*
  exp(-log(gamma + 10^(-2)))

  rho[1] <- theta*r
  r ~ dbeta(rho_prior[1, 1], rho_prior[1, 2])
  gamma ~ dbeta(mtd_prior[1, 1], mtd_prior[1, 2])
}"

    inits <- function() {
      time_init <- rep(NA, length(time_mod))
      time_init[which(!status)] <- time_cens[which(!status)] + 1

      out <- list(r = rbeta(nrow(data$rho_prior),
                            data$rho_prior[, 1], data$rho_prior[, 2]),
                  gamma = rbeta(nrow(data$mtd_prior),
                                data$mtd_prior[, 1], data$mtd_prior[, 2]),
                  time_mod = time_init)
      return(out)
    }

    data_base <- list('time_mod' = time_mod, 'time_cens' = time_cens,
                      'censored' = censored, 'tau' = data$tau,
                      'design_matrix' = data$design_matrix,
                      'theta' = data$theta,
                      'nobs' = length(time_cens[!is.na(time_cens)]),
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
                             c("beta", "gamma", "rho", "shape"),
                           n.iter = n_mcmc, thin = n_thin,
                           n.chains = n_chains)

    beta <- sample[[1]][, 1:2]
    gamma <- sample[[1]][, 3]
    rho <- sample[[1]][, 4]
    shape <- sample[[1]][, 5]

    out <- list(beta = beta, gamma = gamma, rho = rho, shape = shape,
                sample = sample)
  } else {
    sample <- coda.samples(j, variable.names =  c("beta", "gamma", "rho"),
                           n.iter = n_mcmc, thin = n_thin,
                           n.chains = n_chains)

    beta <- sample[[1]][, 1:2]
    gamma <- sample[[1]][, 3]
    rho <- sample[[1]][, 4]

    out <- list(beta = beta, gamma = gamma, rho = rho, sample = sample)
  }
  return(out)
}

