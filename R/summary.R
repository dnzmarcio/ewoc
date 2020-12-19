#'@export
summary.ewoc_d1classical <- function(object, ..., print = TRUE){

  p00 <- data.frame(min_dose = object$trial$min_dose,
                    max_dose = object$trial$max_dose,
                    theta = object$trial$theta,
                    alpha = object$trial$alpha,
                    n = length(object$trial$response))
  colnames(p00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(coda::as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose[1], hpd_dose[2])

  hpd_pdlt <- coda::HPDinterval(coda::as.mcmc(object$pdlt))
  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  prob_dlt <- round(median(object$pdlt), 2)
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
    p02 <- data.frame(estimate = tab02[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Estimate", "95% HPD")
    print(p02)

  } else {

    out <- list(next_dose = next_dose, hpd_dose = hpd_dose,
                prob_dlt = prob_dlt, hpd_pdlt = hpd_pdlt)
    return(out)
  }
}

#'@export
summary.ewoc_d1extended <- function(object, ..., print = TRUE){

  p00 <- data.frame(min_dose = object$trial$min_dose,
                    max_dose = object$trial$max_dose,
                    theta = object$trial$theta,
                    alpha = object$trial$alpha,
                    n = length(object$trial$response))
  colnames(p00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(coda::as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose[1], hpd_dose[2])

  hpd_pdlt <- coda::HPDinterval(coda::as.mcmc(object$pdlt))
  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  prob_dlt <- round(median(object$pdlt), 2)
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
    p02 <- data.frame(estimate = tab02[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Estimate", "95% HPD")
    print(p02)

  } else {

    out <- list(next_dose = next_dose, hpd_dose = hpd_dose,
                prob_dlt = prob_dlt, hpd_pdlt = hpd_pdlt)
    return(out)
  }
}

#'@export
summary.ewoc_d1ph <- function(object, ..., print = TRUE){

  p00 <- data.frame(min_dose = object$trial$min_dose,
                    max_dose = object$trial$max_dose,
                    theta = object$trial$theta,
                    alpha = object$trial$alpha,
                    n = nrow(object$trial$response))
  colnames(p00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(coda::as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose[1], hpd_dose[2])

  hpd_pdlt <- coda::HPDinterval(coda::as.mcmc(object$pdlt))
  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  prob_dlt <- round(median(object$pdlt), 2)
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
    p02 <- data.frame(estimate = tab02[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Estimate", "95% HPD")
    print(p02)

  } else {

    out <- list(next_dose = next_dose, hpd_dose = hpd_dose,
                prob_dlt = prob_dlt, hpd_pdlt = hpd_pdlt)
    return(out)
  }

}

#'@export
summary.ewoc_d1dicov <- function(object, ..., print = TRUE){

  covariable <- factor(object$trial$covariable,
                       levels = object$trial$levels_cov)

  p00 <- data.frame(group = object$trial$levels_cov,
                    min_dose = object$trial$min_dose(object$trial$levels_cov),
                    max_dose = object$trial$max_dose(object$trial$levels_cov),
                    theta = object$trial$theta, alpha = object$trial$alpha,
                    n = as.numeric(table(covariable)))
  colnames(p00) <- c("Group", "Minimum Dose", "Maximum Dose", "Theta",
                     "Alpha", "Number of patients")
  rownames(p00) <- NULL

  hpd_dose <- matrix(HPDinterval(as.mcmc(object$mtd)), ncol = 2)
  hpd_dose <- round(hpd_dose, 2)
  next_dose <- round(as.numeric(object$next_dose), 2)

  tab01 <- data.frame(covariate = object$trial$levels_cov,
                      next_dose = next_dose,
                      lower = hpd_dose[, 1], upper = hpd_dose[, 2])

  hpd_pdlt <- matrix(HPDinterval(as.mcmc(object$pdlt)), ncol = 2)
  hpd_pdlt <- round(hpd_pdlt, 2)
  prob_dlt <- round(median(object$pdlt), 2)

  tab02 <- data.frame(covariate = object$trial$levels_cov,
                      prob_dlt = prob_dlt,
                      lower = hpd_pdlt[, 1],
                      upper = hpd_pdlt[, 2])

  if (print){
    cat("Conditions\n")
    print(p00)
    cat("\n")

    index <- which(object$trial$levels_cov ==
                     object$trial$next_patient_cov)

    cat("Next Dose\n")
    p01 <- data.frame(covariate = object$trial$next_patient_cov,
                      next_dose = next_dose[index],
                      hpd = paste0("(", hpd_dose[index, 1], " ; ",
                                   hpd_dose[index, 2], ")"))
    colnames(p01) <- c("Group", "Estimate", "95% HPD")
    rownames(p01) <- NULL
    print(p01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    p02 <- data.frame(covariate = object$trial$next_patient_cov,
                      prob_dlt = prob_dlt,
                      hpd = paste0("(", hpd_pdlt[index, 1], " ; ",
                                   hpd_pdlt[index, 2], ")"))
    colnames(p02) <- c("Group", "Estimate", "95% HPD")
    rownames(p02) <- NULL
    print(p02)
  } else {

    out <- list(next_dose = tab01[, 1:2], hpd_dose = tab01[, c(1,3,4)],
                prob_dlt = tab02[, 1:2], hpd_pdlt = tab02[, c(1,3,4)])
    return(out)
  }
}

