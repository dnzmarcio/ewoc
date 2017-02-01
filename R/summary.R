#'@export
summary.ewoc_d1basic <- function(object, ..., pdlt = pdlt_d1basic){

  n <- nrow(object$trial$response)
  alpha <- object$trial$alpha
  theta <- object$trial$theta
  max_dose <- object$trial$max_dose
  min_dose <- object$trial$min_dose

  tab <- data.frame(min_dose(), max_dose(), theta, alpha, n)
  colnames(tab) <- c("Minimum Dose", "Maximum Dose", "Theta",
                     "Alpha", "Number of patients")
  cat("Conditions\n")
  print(tab)
  cat("\n")

  hpd_dose <- HPDinterval(as.mcmc(object$mtd))
  hpd <- round(as.numeric(hpd_dose), 2)
  hpd <- paste0("(", round(hpd[1], 2), " ; ", round(hpd[2], 2), ")")
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd)

  prob_dlt <- pdlt(dose = object$next_dose, rho = object$rho,
                   gamma = gamma, theta = theta)
  hpd_pdlt <- HPDinterval(as.mcmc(prob_dlt))

  hpd <- round(as.numeric(hpd_pdlt), 2)
  hpd <- paste0("(", round(hpd[1], 2), " ; ", round(hpd[2], 2), ")")
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd)

  cat("Next Dose\n")
  colnames(tab01) <- c("Group", "Estimate", "95% HPD")
  print(tab01)
  cat("\n")

  cat("P(DLT| next dose)\n")
  colnames(tab02) <- c("Group", "Estimate", "95% HPD")
  print(tab02)
}

#'@export
summary.ewoc_d1extended <- function(object, ..., pdlt = pdlt_d1extended){

  n <- nrow(object$trial$response)
  alpha <- object$trial$alpha
  theta <- object$trial$theta
  max_dose <- object$trial$max_dose
  min_dose <- object$trial$min_dose

  tab <- data.frame(min_dose(), max_dose(), theta, alpha, n)
  colnames(tab) <- c("Minimum Dose", "Maximum Dose", "Theta",
                     "Alpha", "Number of patients")
  cat("Conditions\n")
  print(tab)
  cat("\n")

  hpd_dose <- HPDinterval(as.mcmc(object$mtd))
  hpd <- round(as.numeric(hpd_dose), 2)
  hpd <- paste0("(", round(hpd[1], 2), " ; ", round(hpd[2], 2), ")")
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd)

  prob_dlt <- pdlt(dose = object$next_dose, rho = object$rho,
                   min_dose = min_dose(), max_dose = max_dose())
  hpd_pdlt <- HPDinterval(as.mcmc(prob_dlt))

  hpd <- round(as.numeric(hpd_pdlt), 2)
  hpd <- paste0("(", round(hpd[1], 2), " ; ", round(hpd[2], 2), ")")
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd)

  cat("Next Dose\n")
  colnames(tab01) <- c("Group", "Estimate", "95% HPD")
  print(tab01)
  cat("\n")

  cat("P(DLT| next dose)\n")
  colnames(tab02) <- c("Group", "Estimate", "95% HPD")
  print(tab02)
}

#'@export
summary.ewoc_d1ph <- function(object, ..., pdlt = pdlt_d1ph){

  n <- nrow(object$trial$response)
  alpha <- object$trial$alpha
  theta <- object$trial$theta
  max_dose <- object$trial$max_dose
  min_dose <- object$trial$min_dose

  tab <- data.frame(min_dose(), max_dose(), theta, alpha, n)
  colnames(tab) <- c("Minimum Dose", "Maximum Dose", "Theta",
                     "Alpha", "Number of patients")
  cat("Conditions\n")
  print(tab)
  cat("\n")

  hpd_dose <- HPDinterval(as.mcmc(object$mtd))
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
  hpd_pdlt <- HPDinterval(as.mcmc(prob_dlt))

  hpd <- round(as.numeric(hpd_pdlt), 2)
  hpd <- paste0("(", round(hpd[1], 2), " ; ", round(hpd[2], 2), ")")
  prob_dlt <- round(median(prob_dlt), 2)
  tab02 <- data.frame(prob_dlt, hpd)

  cat("Next Dose\n")
  colnames(tab01) <- c("Estimate", "95% HPD")
  print(tab01)
  cat("\n")

  cat("P(DLT| next dose)\n")
  colnames(tab02) <- c("Estimate", "95% HPD")
  print(tab02)
}

#'@export
summary.ewoc_d1ordinal <- function(object, ..., pdlt = pdlt_d1ordinal){

  npat <- colSums(object$trial$covariable)
  npat[1] <- npat[1] - sum(npat[2:length(npat)])
  alpha <- object$trial$alpha
  theta <- object$trial$theta
  max_dose <- object$trial$max_dose
  min_dose <- object$trial$min_dose

  tab <- data.frame(group = object$trial$levels_cov,
                    min_dose(object$trial$levels_cov),
                    max_dose(object$trial$levels_cov),
                    theta, alpha, npat)
  colnames(tab) <- c("Group", "Minimum Dose", "Maximum Dose", "Theta",
                     "Alpha", "Number of patients")
  cat("Conditions\n")
  print(tab)
  cat("\n")

  next_dose <- rep(NA, length(object$trial$levels_cov))

  for (i in 1:length(object$trial$levels_cov)){

    hpd_dose <- HPDinterval(as.mcmc(object$mtd[, i]))
    hpd <- round(as.numeric(hpd_dose), 2)
    hpd <- paste0("(", round(hpd[1], 2), " ; ", round(hpd[2], 2), ")")
    next_dose[i] <- round(as.numeric(object$next_dose[i]), 2)

    temp01 <- cbind(next_dose[i], hpd)
    covariable <- rep(0, length(object$trial$levels_cov))
    covariable[c(1, i)] <- 1

    prob_dlt <- pdlt(dose = object$next_dose[i], rho = object$rho,
                     gamma = object$gamma, theta = theta,
                     min_dose = min_dose(object$trial$levels_cov),
                     max_dose = max_dose(object$trial$levels_cov),
                     cov = covariable)
    hpd_pdlt <- HPDinterval(as.mcmc(prob_dlt))

    hpd <- round(as.numeric(hpd_pdlt), 2)
    hpd <- paste0("(", round(hpd[1], 2), " ; ", round(hpd[2], 2), ")")
    prob_dlt <- round(median(prob_dlt), 2)

    temp02 <- cbind(prob_dlt, hpd)

    if (i > 1){
      temp01 <- rbind(tab01, temp01)
      temp02 <- rbind(tab02, temp02)
    }

    tab01 <- temp01
    tab02 <- temp02
  }

  cat("Next Dose\n")
  tab01 <- data.frame(group = object$trial$levels_cov, tab01)
  colnames(tab01) <- c("Group", "Estimate", "95% HPD")
  print(tab01)
  cat("\n")

  cat("P(DLT| next dose)\n")
  tab02 <- data.frame(group = object$trial$levels_cov, tab02)
  colnames(tab02) <- c("Group", "Estimate", "95% HPD")
  print(tab02)

}


#'@export
summary.ewoc_d1multinomial <- function(object, ..., pdlt = pdlt_d1multinomial){

  npat <- colSums(object$trial$covariable)
  npat[1] <- npat[1] - sum(npat[2:length(npat)])
  alpha <- object$trial$alpha
  theta <- object$trial$theta
  max_dose <- object$trial$max_dose
  min_dose <- object$trial$min_dose

  tab <- data.frame(group = object$trial$levels_cov,
                    min_dose(object$trial$levels_cov),
                    max_dose(object$trial$levels_cov),
                    theta, alpha, npat)
  colnames(tab) <- c("Group", "Minimum Dose", "Maximum Dose", "Theta",
                     "Alpha", "Number of patients")
  cat("Conditions\n")
  print(tab)
  cat("\n")

  next_dose <- rep(NA, length(object$trial$levels_cov))

  for (i in 1:length(object$trial$levels_cov)){

    hpd_dose <- HPDinterval(as.mcmc(object$mtd[, i]))
    hpd <- round(as.numeric(hpd_dose), 2)
    hpd <- paste0("(", round(hpd[1], 2), " ; ", round(hpd[2], 2), ")")
    next_dose[i] <- round(as.numeric(object$next_dose[i]), 2)

    temp01 <- cbind(next_dose[i], hpd)
    covariable <- rep(0, length(object$trial$levels_cov))
    covariable[c(1, i)] <- 1

    prob_dlt <- pdlt(dose = object$next_dose[i], rho = object$rho,
                     gamma = object$gamma, theta = theta,
                     min_dose = min_dose(object$trial$levels_cov),
                     max_dose = max_dose(object$trial$levels_cov),
                     cov = covariable)
    hpd_pdlt <- HPDinterval(as.mcmc(prob_dlt))

    hpd <- round(as.numeric(hpd_pdlt), 2)
    hpd <- paste0("(", round(hpd[1], 2), " ; ", round(hpd[2], 2), ")")
    prob_dlt <- round(median(prob_dlt), 2)

    temp02 <- cbind(prob_dlt, hpd)

    if (i > 1){
      temp01 <- rbind(tab01, temp01)
      temp02 <- rbind(tab02, temp02)
    }

    tab01 <- temp01
    tab02 <- temp02
  }

  cat("Next Dose\n")
  tab01 <- data.frame(group = object$trial$levels_cov, tab01)
  colnames(tab01) <- c("Group", "Estimate", "95% HPD")
  print(tab01)
  cat("\n")

  cat("P(DLT| next dose)\n")
  tab02 <- data.frame(group = object$trial$levels_cov, tab02)
  colnames(tab02) <- c("Group", "Estimate", "95% HPD")
  print(tab02)

}
