#'@export
summary.ewoc_d1basic <- function(object, ..., pdlt = NULL, print = TRUE){

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

  if (is.null(pdlt))
    pdlt <- pdlt_d1basic(rho = object$rho, mtd = object$mtd,
                         theta = object$trial$theta,
                         min_dose = object$trial$min_dose,
                         max_dose = object$trial$max_dose)

  prob_dlt <- pdlt(object$next_dose)
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
summary.ewoc_d1extended <- function(object, ..., pdlt = NULL, print = TRUE){

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

  if (is.null(pdlt))
    pdlt <- pdlt_d1extended(rho = object$rho,
                            min_dose = object$trial$min_dose,
                            max_dose = object$trial$max_dose)

  prob_dlt <- pdlt(dose = object$next_dose)
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
summary.ewoc_d1ph <- function(object, ..., pdlt = NULL, print = TRUE){

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

  if (is.null(pdlt))
    pdlt <- pdlt_d1ph(rho = object$rho, mtd = object$mtd, shape = object$shape,
                      theta = object$trial$theta,
                      min_dose = object$trial$min_dose,
                      max_dose = object$trial$max_dose,
                      tau = object$trial$tau,
                      distribution = object$trial$distribution)

  prob_dlt <- pdlt(object$next_dose)
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
summary.ewoc_d1multinomial <- function(object, ..., pdlt = NULL, print = TRUE){

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

  covariable <- factor(object$trial$next_patient_cov,
                       levels = object$trial$levels_cov)


  for (i in 1:length(covariable)){
    index <- which(object$trial$next_patient_cov[i] == object$trial$levels_cov)

    hpd_dose <- matrix(HPDinterval(as.mcmc(object$mtd)), ncol = 2)
    hpd_dose <- round(hpd_dose[index, ], 2)
    next_dose <- round(as.numeric(object$next_dose[index]), 2)

    temp <- cbind(next_dose, hpd_dose[1], hpd_dose[2])

    if (i > 1)
      temp <- rbind(tab01, temp)

    tab01 <- temp

    if (is.null(pdlt))
      pdlt <- pdlt_d1multinomial(rho = object$rho, mtd = object$mtd,
                                 theta = object$trial$theta,
                                 min_dose = object$trial$min_dose,
                                 max_dose = object$trial$max_dose,
                                 levels_cov = object$trial$levels_cov)

    prob_dlt <- pdlt(dose = object$next_dose[index], cov = covariable[i])
    hpd_pdlt <- matrix(HPDinterval(coda::as.mcmc(prob_dlt)), ncol = 2)
    hpd_pdlt <- round(hpd_pdlt, 2)
    prob_dlt <- round(median(prob_dlt), 2)

    temp <- cbind(prob_dlt, hpd_pdlt[1], hpd_pdlt[2])

    if (i > 1)
      temp <- rbind(tab02, temp)

    tab02 <- temp
  }

  if (print){
    cat("Conditions\n")
    print(p00)
    cat("\n")

    cat("Next Dose\n")
    p01 <- data.frame(group = object$trial$next_patient_cov,
                      estimate = tab01[, 1],
                      hpd = paste0("(", tab01[, 2], " ; ", tab01[, 3], ")"))
    colnames(p01) <- c("Group", "Estimate", "95% HPD")
    rownames(p01) <- NULL
    print(p01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    p02 <- data.frame(group = object$trial$next_patient_cov,
                      estimate = tab02[, 1],
                      hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Group", "Estimate", "95% HPD")
    rownames(p02) <- NULL
    print(p02)
  } else {

    out <- list(next_dose = tab01[, 1, drop=F], hpd_dose = tab01[, 2:3, drop=F],
                prob_dlt = tab02[, 1, drop=F], hpd_pdlt = tab02[, 2:3, drop=F])
    return(out)
  }
}

#'@export
summary.ewoc_d1ordinal <- function(object, ..., pdlt = NULL, print = TRUE){

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

  covariable <- factor(object$trial$next_patient_cov,
                       levels = object$trial$levels_cov)


  for (i in 1:length(covariable)){
    index <- which(object$trial$next_patient_cov[i] == object$trial$levels_cov)

    hpd_dose <- matrix(HPDinterval(as.mcmc(object$mtd)), ncol = 2)
    hpd_dose <- round(hpd_dose[index, ], 2)
    next_dose <- round(as.numeric(object$next_dose[index]), 2)

    temp <- cbind(next_dose, hpd_dose[1], hpd_dose[2])

    if (i > 1)
      temp <- rbind(tab01, temp)

    tab01 <- temp

    if (is.null(pdlt))
      pdlt <- pdlt_d1ordinal(rho = object$rho, mtd = object$mtd,
                             theta = object$trial$theta,
                             min_dose = object$trial$min_dose,
                             max_dose = object$trial$max_dose,
                             levels_cov = object$trial$levels_cov)

    prob_dlt <- pdlt(dose = object$next_dose[index], cov = covariable[i])
    hpd_pdlt <- matrix(HPDinterval(coda::as.mcmc(prob_dlt)), ncol = 2)
    hpd_pdlt <- round(hpd_pdlt, 2)
    prob_dlt <- round(median(prob_dlt), 2)

    temp <- cbind(prob_dlt, hpd_pdlt[1], hpd_pdlt[2])

    if (i > 1)
      temp <- rbind(tab02, temp)

    tab02 <- temp
  }

  if (print){
    cat("Conditions\n")
    print(p00)
    cat("\n")

    cat("Next Dose\n")
    p01 <- data.frame(group = object$trial$next_patient_cov,
                      estimate = tab01[, 1],
                      hpd = paste0("(", tab01[, 2], " ; ", tab01[, 3], ")"))
    colnames(p01) <- c("Group", "Estimate", "95% HPD")
    rownames(p01) <- NULL
    print(p01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    p02 <- data.frame(group = object$trial$next_patient_cov,
                      estimate = tab02[, 1],
                      hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Group", "Estimate", "95% HPD")
    rownames(p02) <- NULL
    print(p02)
  } else {

    out <- list(next_dose = tab01[, 1, drop=F], hpd_dose = tab01[, 2:3, drop=F],
                prob_dlt = tab02[, 1, drop=F], hpd_pdlt = tab02[, 2:3, drop=F])
    return(out)
  }
}

#'@export
summary.ewoc_d1continuous <- function(object, ..., pdlt = NULL, print = TRUE){

  p00 <- data.frame(covariable = unique(object$trial$covariable),
                    min_dose = object$trial$min_dose(unique(object$trial$covariable)),
                    max_dose = object$trial$max_dose(unique(object$trial$covariable)),
                    theta = object$trial$theta, alpha = object$trial$alpha,
                    n = as.numeric(table(object$trial$covariable)))
  colnames(p00) <- c("Covariable", "Minimum Dose", "Maximum Dose", "Theta",
                     "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(coda::as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  tab01 <- data.frame(round(object$next_dose, 2), hpd_dose[1], hpd_dose[2])

  if (is.null(pdlt))
    prob_dlt <- pdlt_d1continuous(rho = object$rho,
                                  mtd = object$mtd,
                                  theta = object$trial$theta,
                                  min_dose = object$trial$min_dose(object$trial$next_patient_cov),
                                  max_dose = object$trial$max_dose(object$trial$next_patient_cov),
                                  min_cov = object$trial$min_cov,
                                  max_cov = object$trial$max_cov,
                                  cov = object$trial$next_patient_cov,
                                  direction = object$trial$direction)

  prob_dlt <- pdlt(dose = object$next_dose)
  hpd_pdlt <- coda::HPDinterval(coda::as.mcmc(prob_dlt))

  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd_pdlt[1], hpd_pdlt[2])

  if (print){
    cat("Conditions\n")
    print(p00)
    cat("\n")

    cat("Next Dose\n")
    p01 <- data.frame(covariable = object$trial$next_patient_cov,
                      estimate = tab01[, 1],
                      hpd = paste0("(", tab01[, 2], " ; ", tab01[, 3], ")"))
    colnames(p01) <- c("Covariable", "Estimate", "95% HPD")
    print(p01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    p02 <- data.frame(covariable = object$trial$next_patient_cov,
                      estimate = tab02[, 1],
                      hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Covariable", "Estimate", "95% HPD")
    print(p02)
  } else {
    out <- list(next_dose = tab01[, 1, drop=F], hpd_dose = tab01[, 2:3, drop=F],
                prob_dlt = tab02[, 1, drop=F], hpd_pdlt = tab02[, 2:3, drop=F])

    return(out)
  }
}


#'@export
summary.ewoc_d1excontinuous <- function(object, ..., pdlt = NULL, print = TRUE){

  p00 <- data.frame(covariable = unique(object$trial$covariable),
                    min_dose = object$trial$min_dose(unique(object$trial$covariable)),
                    max_dose = object$trial$max_dose(unique(object$trial$covariable)),
                    theta = object$trial$theta, alpha = object$trial$alpha,
                    n = as.numeric(table(object$trial$covariable)))
  colnames(p00) <- c("Covariable", "Minimum Dose", "Maximum Dose", "Theta",
                     "Alpha", "Number of patients")

  hpd_dose <- coda::HPDinterval(coda::as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  tab01 <- data.frame(round(object$next_dose, 2), hpd_dose[1], hpd_dose[2])

  if (is.null(pdlt))
    pdlt <- pdlt_d1excontinuous(rho = object$rho,
                                theta = object$trial$theta,
                                min_dose = object$trial$min_dose(object$trial$next_patient_cov),
                                max_dose = object$trial$max_dose(object$trial$next_patient_cov),
                                min_cov = object$trial$min_cov,
                                max_cov = object$trial$max_cov,
                                cov = object$trial$next_patient_cov,
                                direction = object$trial$direction)

  prob_dlt <- pdlt(dose = object$next_dose)
  hpd_pdlt <- coda::HPDinterval(coda::as.mcmc(prob_dlt))

  hpd_pdlt <- round(as.numeric(hpd_pdlt), 2)
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd_pdlt[1], hpd_pdlt[2])

  if (print){
    cat("Conditions\n")
    print(p00)
    cat("\n")

    cat("Next Dose\n")
    p01 <- data.frame(covariable = object$trial$next_patient_cov,
                      estimate = tab01[, 1],
                      hpd = paste0("(", tab01[, 2], " ; ", tab01[, 3], ")"))
    colnames(p01) <- c("Covariable", "Estimate", "95% HPD")
    print(p01)
    cat("\n")

    cat("P(DLT| next dose)\n")
    p02 <- data.frame(covariable = object$trial$next_patient_cov,
                      estimate = tab02[, 1],
                      hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
    colnames(p02) <- c("Covariable", "Estimate", "95% HPD")
    print(p02)
  } else {
    out <- list(next_dose = tab01[, 1], hpd_dose = tab01[, 2:3],
                prob_dlt = tab02[, 1], hpd_pdlt = tab02[, 2:3])

    return(out)
  }
}

