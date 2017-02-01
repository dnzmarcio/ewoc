#'@export
summary.ewoc_d1basic <- function(object, ..., pdlt = pdlt_d1basic, print = TRUE){

  p00 <- data.frame(object$trial$min_dose(),
                    object$trial$max_dose(),
                    object$trial$theta, object$trial$alpha,
                    nrow(object$trial$response))
  colnames(p00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(coda::as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose[1], hpd_dose[2])

  prob_dlt <- pdlt(dose = next_dose, rho = object$rho,
                   gamma = object$gamma, theta = object$trial$theta,
                   min_dose = object$trial$min_dose(),
                   max_dose = object$trial$max_dose())
  hpd_pdlt <- coda::HPDinterval(coda::as.mcmc(prob_dlt))

  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd_pdlt[1], hpd_pdlt[2])

  out <- list(next_dose = next_dose, hpd_dose = hpd_dose,
              prob_dlt = prob_dlt, hpd_pdlt = hpd_pdlt)

  if (print){
    cat("Conditions\n")
    print(p00)
    cat("\n")

    cat("Next Dose\n")
    p01 <- data.frame(estimate = tab01[, 1],
                        hpd = paste0("(", tab01[, 2], " ; ", tab01[, 3], ")"))
    colnames(p01) <- c("Estimate", "95% HPD")
    print(p01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    p02 <- data.frame(estimate = tab01[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Estimate", "95% HPD")
    print(p02)
  }

}

#'@export
summary.ewoc_d1extended <- function(object, ..., pdlt = pdlt_d1extended,
                                    print = TRUE){

  p00 <- data.frame(object$trial$min_dose(), object$trial$max_dose(),
                    object$trial$theta, object$trial$alpha,
                    nrow(object$trial$response))
  colnames(p00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(coda::as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose[1], hpd_dose[2])

  prob_dlt <- pdlt(dose = next_dose, rho = object$rho,
                   min_dose = object$trial$min_dose(),
                   max_dose = object$trial$max_dose())
  hpd_pdlt <- coda::HPDinterval(coda::as.mcmc(prob_dlt))

  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd_pdlt[1], hpd_pdlt[2])

  out <- list(next_dose = next_dose, hpd_dose = hpd_dose,
              prob_dlt = prob_dlt, hpd_pdlt = hpd_pdlt)

  if (print){
    cat("Conditions\n")
    print(p00)
    cat("\n")

    cat("Next Dose\n")
    p01 <- data.frame(estimate = tab01[, 1],
                        hpd = paste0("(", tab01[, 2], " ; ", tab01[, 3], ")"))
    colnames(p01) <- c("Estimate", "95% HPD")
    print(p01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    p02 <- data.frame(estimate = tab01[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Estimate", "95% HPD")
    print(p02)
  }
}

#'@export
summary.ewoc_d1ph <- function(object, ..., pdlt = pdlt_d1ph, print = TRUE){

  p00 <- data.frame(object$trial$min_dose(), object$trial$max_dose(),
                    object$trial$theta, object$trial$alpha,
                    nrow(object$trial$response))
  colnames(p00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(coda::as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose[1], hpd_dose[2])

  prob_dlt <- pdlt(dose = next_dose, rho = object$rho,
                   gamma = object$gamma, shape = object$shape,
                   theta = object$trial$theta,
                   min_dose = object$trial$min_dose(),
                   max_dose = object$trial$max_dose(),
                   tau = object$trial$tau,
                   distribution = object$trial$distribution)
  hpd_pdlt <- coda::HPDinterval(coda::as.mcmc(prob_dlt))

  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd_pdlt[1], hpd_pdlt[2])

  if (print){
    cat("Conditions\n")
    print(p00)
    cat("\n")

    cat("Next Dose\n")
    p01 <- data.frame(estimate = tab01[, 1],
                        hpd = paste0("(", tab01[, 2], " ; ", tab01[, 3], ")"))
    colnames(p01) <- c("Estimate", "95% HPD")
    print(p01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    p02 <- data.frame(estimate = tab01[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Estimate", "95% HPD")
    print(p02)
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

  if (is.null(next_covariable))
    next_covariable <- object$trial$levels_cov

  p00 <- data.frame(group = object$trial$levels_cov,
                    object$trial$min_dose(object$trial$levels_cov),
                    object$trial$max_dose(object$trial$levels_cov),
                      object$trial$theta, object$trial$alpha, npatients)
  colnames(p00) <- c("Group", "Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  index <- which(object$trial$levels_cov ==  next_covariable)
  counter <- 1

  for (i in index){
    hpd_dose <- coda::HPDinterval(coda::as.mcmc(object$mtd[, i]))
    hpd_dose <- round(as.numeric(hpd_dose), 2)
    next_dose <- round(as.numeric(object$next_dose[i]), 2)

    temp01 <- cbind(next_dose, hpd_dose[1], hpd_dose[2])
    covariable <- rep(0, length(object$trial$levels_cov))
    covariable[c(1, i)] <- 1

    prob_dlt <- pdlt(dose = next_dose, rho = object$rho,
                     gamma = object$gamma, theta = object$trial$theta,
                     min_dose = object$trial$min_dose(object$trial$levels_cov)[i],
                     max_dose = object$trial$max_dose(object$trial$levels_cov)[i],
                     cov = covariable)
    hpd_pdlt <- coda::HPDinterval(coda::as.mcmc(prob_dlt))

    hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
    prob_dlt <- round(median(prob_dlt), 2)

    temp02 <- cbind(prob_dlt, hpd_pdlt[1], hpd_pdlt[2])

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
    print(p00)
    cat("\n")

    cat("Next Dose\n")
    p01 <- data.frame(group = next_covariable, estimate = tab01[, 1],
                        hpd = paste0("(", tab01[, 2], " ; ", tab01[, 3], ")"))
    colnames(p01) <- c("Group", "Estimate", "95% HPD")
    print(p01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    p02 <- data.frame(group = next_covariable, estimate = tab02[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Group", "Estimate", "95% HPD")
    print(p02)
  }

  out <- list(next_dose = tab01[, 1], hpd_dose = tab01[, 2:3],
              prob_dlt = tab02[, 1], hpd_pdlt = tab02[, 2:3])
}

#'@export
summary.ewoc_d1ordinal <- function(object, ..., next_covariable = NULL,
                                   pdlt = pdlt_d1ordinal,
                                   print = TRUE){

  npatients <- colSums(object$trial$covariable)
  npatients[1] <- npatients[1] - sum(npatients[2:length(npatients)])

  if (is.null(next_covariable))
    next_covariable <- object$trial$levels_cov

  p00 <- data.frame(group = object$trial$levels_cov,
                    object$trial$min_dose(object$trial$levels_cov),
                    object$trial$max_dose(object$trial$levels_cov),
                      object$trial$theta, object$trial$alpha, npatients)
  colnames(p00) <- c("Group", "Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  index <- which(object$trial$levels_cov ==  next_covariable)
  counter <- 1

  for (i in index){

    hpd_dose <- coda::HPDinterval(coda::as.mcmc(object$mtd[, i]))
    hpd_dose <- round(as.numeric(hpd_dose), 2)
    next_dose <- round(as.numeric(object$next_dose[i]), 2)

    temp01 <- cbind(next_dose, hpd_dose[1], hpd_dose[2])
    covariable <- rep(0, length(object$trial$levels_cov))
    covariable[c(1, i)] <- 1

    prob_dlt <- pdlt(dose = object$next_dose, rho = object$rho,
                     gamma = object$gamma, theta = object$trial$theta,
                     min_dose = object$trial$min_dose(object$trial$levels_cov)[i],
                     max_dose = object$trial$max_dose(object$trial$levels_cov)[i],
                     cov = covariable)
    hpd_pdlt <- coda::HPDinterval(coda::as.mcmc(prob_dlt))
    hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
    prob_dlt <- round(median(prob_dlt), 2)

    temp02 <- cbind(prob_dlt, hpd_pdlt[1], hpd_pdlt[2])

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
    print(p00)
    cat("\n")

    cat("Next Dose\n")
    p01 <- data.frame(group = next_covariable, estimate = tab01[, 1],
                        hpd = paste0("(", tab01[, 2], " ; ", tab01[, 3], ")"))
    colnames(p01) <- c("Group", "Estimate", "95% HPD")
    print(p01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    p02 <- data.frame(group = next_covariable, estimate = tab02[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Group", "Estimate", "95% HPD")
    print(p02)
  }

  out <- list(next_dose = tab01[, 1], hpd_dose = tab01[, 2:3],
              prob_dlt = tab02[, 1], hpd_pdlt = tab02[, 2:3])
}
