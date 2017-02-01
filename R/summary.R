#'@export
summary.ewoc_d1basic <- function(object, ..., pdlt = pdlt_d1basic, print = TRUE){

  tab00 <- data.frame(object$trial$min_dose(),
                    object$trial$max_dose(),
                    object$trial$theta, object$trial$alpha,
                    nrow(object$trial$response))
  colnames(tab00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  hpd_dose <- paste0("(", round(hpd_dose[1], 2), " ; ",
                     round(hpd_dose[2], 2), ")")
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose)

  prob_dlt <- pdlt(dose = object$next_dose, rho = object$rho,
                   gamma = object$gamma, theta = object$trial$theta,
                   min_dose = object$trial$min_dose(),
                   max_dose = object$trial$max_dose())
  hpd_pdlt <- coda::HPDinterval(as.mcmc(prob_dlt))

  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  hpd_pdlt <- paste0("(", round(hpd_pdlt[1], 2), " ; ",
                     round(hpd_pdlt[2], 2), ")")
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd_pdlt)

  if (print){
    cat("Conditions\n")
    print(tab)
    cat("\n")

    cat("Next Dose\n")
    colnames(tab01) <- c("Estimate", "95% HPD")
    print(tab01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    colnames(tab02) <- c("Estimate", "95% HPD")
    print(tab02)
  }

  out <- list(next_dose = next_dose, hpd_dose = hpd_dose,
              prob_dlt = prob_dlt, hpd_pdlt = hpd_pdlt)
}

#'@export
summary.ewoc_d1extended <- function(object, ..., pdlt = pdlt_d1extended,
                                    print = TRUE){

  n <- nrow(object$trial$response)
  alpha <- object$trial$alpha
  theta <- object$trial$theta
  max_dose <- object$trial$max_dose
  min_dose <- object$trial$min_dose

  tab00 <- data.frame(min_dose(), max_dose(), theta, alpha, n)
  colnames(tab00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  hpd_dose <- paste0("(", round(hpd_dose[1], 2), " ; ",
                     round(hpd_dose[2], 2), ")")
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose)

  prob_dlt <- pdlt(dose = object$next_dose, rho = object$rho,
                   min_dose = min_dose(), max_dose = max_dose())
  hpd_pdlt <- coda::HPDinterval(as.mcmc(prob_dlt))

  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  hpd_pdlt <- paste0("(", round(hpd_pdlt[1], 2), " ; ",
                     round(hpd_pdlt[2], 2), ")")
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd_pdlt)

  if (print){
    cat("Conditions\n")
    print(tab)
    cat("\n")

    cat("Next Dose\n")
    colnames(tab01) <- c("Estimate", "95% HPD")
    print(tab01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    colnames(tab02) <- c("Estimate", "95% HPD")
    print(tab02)
  }

  out <- list(next_dose = next_dose, hpd_dose = hpd_dose,
              prob_dlt = prob_dlt, hpd_pdlt = hpd_pdlt)
}

#'@export
summary.ewoc_d1ph <- function(object, ..., pdlt = pdlt_d1ph, print = TRUE){

  n <- nrow(object$trial$response)
  alpha <- object$trial$alpha
  theta <- object$trial$theta
  max_dose <- object$trial$max_dose
  min_dose <- object$trial$min_dose

  tab00 <- data.frame(min_dose(), max_dose(), theta, alpha, n)
  colnames(tab00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(as.mcmc(object$mtd))
  hpd <- round(as.numeric(hpd_dose), 2)
  hpd <- paste0("(", round(hpd[1], 2), " ; ", round(hpd[2], 2), ")")
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd)

  prob_dlt <- pdlt(dose = object$next_dose, rho = object$rho,
                   gamma = object$gamma, shape = object$shape,
                   theta = theta,
                   min_dose = min_dose(),
                   max_dose = max_dose(),
                   tau = object$trial$tau,
                   distribution = object$trial$distribution)
  hpd_pdlt <- coda::HPDinterval(as.mcmc(prob_dlt))

  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  hpd_pdlt <- paste0("(", round(hpd_pdlt[1], 2), " ; ",
                     round(hpd_pdlt[2], 2), ")")
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd_pdlt)

  if (print){
    cat("Conditions\n")
    print(tab)
    cat("\n")

    cat("Next Dose\n")
    colnames(tab01) <- c("Estimate", "95% HPD")
    print(tab01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    colnames(tab02) <- c("Estimate", "95% HPD")
    print(tab02)
  }

  out <- list(next_dose = next_dose, hpd_dose = hpd_dose,
              prob_dlt = prob_dlt, hpd_pdlt = hpd_pdlt)
}

#'@export
summary.ewoc_d1multinomial <- function(object, ..., next_covariable = NULL,
                                       pdlt = pdlt_d1multinomial,
                                       print = TRUE){

  npatients <- colSums(object$trial$covariable)
  npatients[1] <- npatients[1] - sum(npatients[2:length(npatients)])
  alpha <- object$trial$alpha
  theta <- object$trial$theta
  max_dose <- object$trial$max_dose
  min_dose <- object$trial$min_dose

  if (is.null(next_covariable))
    next_covariable <- object$trial$levels_cov

  tab00 <- data.frame(group = object$trial$levels_cov,
                      min_dose(object$trial$levels_cov),
                      max_dose(object$trial$levels_cov),
                      theta, alpha, npatients)
  colnames(tab00) <- c("Group", "Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  next_dose <- rep(NA, length(next_covariable))
  index <- which(object$trial$levels_cov ==  next_covariable)
  counter <- 1

  for (i in index){

    hpd_dose <- coda::HPDinterval(as.mcmc(object$mtd[, i]))
    hpd_dose <- round(as.numeric(hpd_dose), 2)
    hpd_dose <- paste0("(", round(hpd_dose[1], 2), " ; ",
                       round(hpd_dose[2], 2), ")")
    next_dose[i] <- round(as.numeric(object$next_dose[i]), 2)

    temp01 <- cbind(next_dose[i], hpd_dose)
    covariable <- rep(0, length(object$trial$levels_cov))
    covariable[c(1, i)] <- 1

    prob_dlt <- pdlt(dose = object$next_dose[i], rho = object$rho,
                     gamma = object$gamma, theta = theta,
                     min_dose = min_dose(object$trial$levels_cov)[i],
                     max_dose = max_dose(object$trial$levels_cov)[i],
                     cov = covariable)
    hpd_pdlt <- coda::HPDinterval(as.mcmc(prob_dlt))

    hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
    hpd_pdlt <- paste0("(", round(hpd_pdlt[1], 2), " ; ",
                       round(hpd_pdlt[2], 2), ")")
    prob_dlt <- round(median(prob_dlt), 2)

    temp02 <- cbind(prob_dlt, hpd_pdlt)

    if (counter > 1){
      temp01 <- rbind(tab01, temp01)
      temp02 <- rbind(tab02, temp02)
    }

    tab01 <- temp01
    tab02 <- temp02
    counter <- counter + 1
  }

  if (print){
    cat("Conditions\n")
    print(tab00)
    cat("\n")

    cat("Next Dose\n")
    tab01 <- data.frame(group = next_covariable, tab01)
    colnames(tab01) <- c("Group", "Estimate", "95% HPD")
    print(tab01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    tab02 <- data.frame(group = next_covariable, tab02)
    colnames(tab02) <- c("Group", "Estimate", "95% HPD")
    print(tab02)
  }

  out <- list(next_dose = next_dose, hpd_dose = hpd_dose,
              prob_dlt = prob_dlt, hpd_pdlt = hpd_pdlt)

}

#'@export
summary.ewoc_d1ordinal <- function(object, ..., next_covariable = NULL,
                                   pdlt = pdlt_d1ordinal,
                                   print = TRUE){

  npatients <- colSums(object$trial$covariable)
  npatients[1] <- npatients[1] - sum(npatients[2:length(npatients)])
  alpha <- object$trial$alpha
  theta <- object$trial$theta
  max_dose <- object$trial$max_dose
  min_dose <- object$trial$min_dose

  if (is.null(next_covariable))
    next_covariable <- object$trial$levels_cov

  tab00 <- data.frame(group = object$trial$levels_cov,
                      min_dose(object$trial$levels_cov),
                      max_dose(object$trial$levels_cov),
                      theta, alpha, npatients)
  colnames(tab00) <- c("Group", "Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  next_dose <- rep(NA, length(covariable))
  index <- which(object$trial$levels_cov ==  next_covariable)
  counter <- 1

  for (i in index){

    hpd_dose <- coda::HPDinterval(as.mcmc(object$mtd[, i]))
    hpd_dose <- round(as.numeric(hpd_dose), 2)
    hpd_dose <- paste0("(", round(hpd_dose[1], 2), " ; ",
                       round(hpd_dose[2], 2), ")")
    next_dose[i] <- round(as.numeric(object$next_dose[i]), 2)

    temp01 <- cbind(next_dose[i], hpd_dose)
    covariable <- rep(0, length(object$trial$levels_cov))
    covariable[c(1, i)] <- 1

    prob_dlt <- pdlt(dose = object$next_dose[i], rho = object$rho,
                     gamma = object$gamma, theta = theta,
                     min_dose = min_dose(object$trial$levels_cov)[i],
                     max_dose = max_dose(object$trial$levels_cov)[i],
                     cov = covariable)
    hpd_pdlt <- coda::HPDinterval(as.mcmc(prob_dlt))

    hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
    hpd_pdlt <- paste0("(", round(hpd_pdlt[1], 2), " ; ",
                       round(hpd_pdlt[2], 2), ")")
    prob_dlt <- round(median(prob_dlt), 2)

    temp02 <- cbind(prob_dlt, hpd_pdlt)

    if (counter > 1){
      temp01 <- rbind(tab01, temp01)
      temp02 <- rbind(tab02, temp02)
    }

    tab01 <- temp01
    tab02 <- temp02
    counter <- counter + 1
  }

  if (print){
    cat("Conditions\n")
    print(tab00)
    cat("\n")

    cat("Next Dose\n")
    tab01 <- data.frame(group = next_covariable, tab01)
    colnames(tab01) <- c("Group", "Estimate", "95% HPD")
    print(tab01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    tab02 <- data.frame(group = next_covariable, tab02)
    colnames(tab02) <- c("Group", "Estimate", "95% HPD")
    print(tab02)
  }

  out <- list(next_dose = next_dose, hpd_dose = hpd_dose,
              prob_dlt = prob_dlt, hpd_pdlt = hpd_pdlt)
}
