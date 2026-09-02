#'@importFrom coda as.mcmc HPDinterval
#'@export
summary.ewoc_d1classical <- function(object, ..., print = TRUE){

  p00 <- data.frame(min_dose = object$trial$min_dose,
                    max_dose = object$trial$max_dose,
                    theta = object$trial$theta,
                    alpha = object$trial$alpha,
                    n = length(object$trial$response))
  colnames(p00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- HPDinterval(as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose[1], hpd_dose[2])

  hpd_pdlt_next_dose <- HPDinterval(as.mcmc(object$pdlt))
  hpd_pdlt_next_dose <- round(as.numeric(hpd_pdlt_next_dose), 2)
  prob_dlt_next_dose <- round(median(object$pdlt), 2)
  tab02 <- data.frame(prob_dlt_next_dose, hpd_pdlt_next_dose[1], hpd_pdlt_next_dose[2])

  if(object$trial$type == "discrete"){
    pdlt_dose <- matrix(NA,
                        ncol = length(object$trial$dose_set),
                        nrow = length(object$mtd))

    for(i in 1:length(object$mtd)){
      pdlt_function <- pdlt_d1classical(rho = object$rho[i],
                                        mtd = object$mtd[i],
                                        theta = object$trial$theta,
                                        min_dose = object$trial$min_dose,
                                        max_dose = object$trial$max_dose)
      pdlt_dose[i, ] <- sapply(object$trial$dose_set, pdlt_function)
    }

    prob_dlt_dose <- apply(pdlt_dose, 2, median)
    hpd_pdlt_dose <- apply(pdlt_dose, 2, function(x) HPDinterval(as.mcmc(x)))
    tab03 <- data.frame(dose = object$trial$dose_set, prob_dlt_dose,
                        hpd_pdlt_dose[1, ], hpd_pdlt_dose[2, ])
  } else {
    prob_dlt_dose <- NULL
    hpd_pdlt_dose <- NULL
  }


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

    if(object$trial$type == "continuous"){
      cat("P(DLT| next dose)\n")
      p02 <- data.frame(estimate = tab02[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
      colnames(p02) <- c("Estimate", "95% HPD")
      print(p02)
      cat("\n")

    } else {
      cat("P(DLT| dose)\n")
      p03 <- data.frame(dose = tab03[, 1], estimate = round(tab03[, 2], 2),
                        hpd = paste0("(", round(tab03[, 3], 2), " ; ",
                                     round(tab03[, 4], 2), ")"))
      colnames(p03) <- c("Dose", "Estimate", "95% HPD")
      print(p03)
    }

  } else {

    out <- list(next_dose = next_dose,
                hpd_dose = hpd_dose,
                prob_dlt_next_dose = prob_dlt_next_dose,
                hpd_pdlt_next_dose = hpd_pdlt_next_dose,
                prob_dlt_dose = prob_dlt_dose,
                hpd_pdlt_dose = hpd_pdlt_dose)
    return(out)
  }
}

#'@importFrom coda as.mcmc HPDinterval
#'@export
summary.ewoc_d1extended <- function(object, ..., print = TRUE){

  p00 <- data.frame(min_dose = object$trial$min_dose,
                    max_dose = object$trial$max_dose,
                    theta = object$trial$theta,
                    alpha = object$trial$alpha,
                    n = length(object$trial$response))
  colnames(p00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- HPDinterval(as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose[1], hpd_dose[2])

  hpd_pdlt_next_dose <- HPDinterval(as.mcmc(object$pdlt))
  hpd_pdlt_next_dose <- round(as.numeric(hpd_pdlt_next_dose), 2)
  prob_dlt_next_dose <- round(median(object$pdlt), 2)
  tab02 <- data.frame(prob_dlt_next_dose, hpd_pdlt_next_dose[1], hpd_pdlt_next_dose[2])

  if(object$trial$type == "discrete"){
    pdlt_dose <- matrix(NA,
                            ncol = length(object$trial$dose_set),
                            nrow = length(object$mtd))

    for(i in 1:length(object$mtd)){
      pdlt_function <- pdlt_d1extended(rho = object$rho[i, ],
                                       min_dose = object$trial$min_dose,
                                       max_dose = object$trial$max_dose)
      pdlt_dose[i, ] <- sapply(object$trial$dose_set, pdlt_function)
    }

    prob_dlt_dose <- apply(pdlt_dose, 2, median)
    hpd_pdlt_dose <- apply(pdlt_dose, 2, function(x) HPDinterval(as.mcmc(x)))
    tab03 <- data.frame(object$trial$dose_set, prob_dlt_dose,
                        hpd_pdlt_dose[1, ], hpd_pdlt_dose[2, ])
  } else {
    prob_dlt_dose <- NULL
    hpd_pdlt_dose <- NULL
  }

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

    if(object$trial$type == "continuous"){
      cat("P(DLT| next dose)\n")
      p02 <- data.frame(estimate = tab02[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
      colnames(p02) <- c("Estimate", "95% HPD")
      print(p02)
      cat("\n")

    } else {
      cat("P(DLT| dose)\n")
      p03 <- data.frame(dose = tab03[, 1], estimate = round(tab03[, 2], 2),
                        hpd = paste0("(", round(tab03[, 3], 2), " ; ",
                                     round(tab03[, 4], 2), ")"))
      colnames(p03) <- c("Dose", "Estimate", "95% HPD")
      print(p03)
    }

  } else {

    out <- list(next_dose = next_dose,
                hpd_dose = hpd_dose,
                prob_dlt_next_dose = prob_dlt_next_dose,
                hpd_pdlt_next_dose = hpd_pdlt_next_dose,
                prob_dlt_dose = prob_dlt_dose,
                hpd_pdlt_dose = hpd_pdlt_dose)
    return(out)
  }
}

#'@importFrom coda as.mcmc HPDinterval
#'@export
summary.ewoc_d1ph <- function(object, ..., print = TRUE){

  p00 <- data.frame(min_dose = object$trial$min_dose,
                    max_dose = object$trial$max_dose,
                    theta = object$trial$theta,
                    alpha = object$trial$alpha,
                    n = nrow(object$trial$response))
  colnames(p00) <- c("Minimum Dose", "Maximum Dose", "Theta",
                       "Alpha", "Number of patients")

  hpd_dose <- HPDinterval(as.mcmc(object$mtd))
  hpd_dose <- round(as.numeric(hpd_dose), 2)
  next_dose <- round(as.numeric(object$next_dose), 2)
  tab01 <- data.frame(next_dose, hpd_dose[1], hpd_dose[2])

  hpd_pdlt_next_dose <- HPDinterval(as.mcmc(object$pdlt))
  hpd_pdlt_next_dose <- round(as.numeric(hpd_pdlt_next_dose), 2)
  prob_dlt_next_dose <- round(median(object$pdlt), 2)
  tab02 <- data.frame(prob_dlt_next_dose, hpd_pdlt_next_dose[1], hpd_pdlt_next_dose[2])

  if(object$trial$type == "discrete"){
    pdlt_dose <- matrix(NA,
                            ncol = length(object$trial$dose_set),
                            nrow = length(object$mtd))

    for(i in 1:length(object$mtd)){

      if (object$trial$distribution == "exponential"){
        pdlt_function <- pdlt_d1ph(rho = object$rho[i],
                                   mtd = object$mtd[i],
                                   shape = NULL,
                                   theta = object$trial$theta,
                                   min_dose = object$trial$min_dose,
                                   max_dose = object$trial$max_dose,
                                   tau = object$trial$tau,
                                   distribution = object$trial$distribution)
      } else {
        pdlt_function <- pdlt_d1ph(rho = object$rho[i],
                                   mtd = object$mtd[i],
                                   shape = object$shape[i],
                                   theta = object$trial$theta,
                                   min_dose = object$trial$min_dose,
                                   max_dose = object$trial$max_dose,
                                   tau = object$trial$tau,
                                   distribution = object$trial$distribution)
      }

      pdlt_dose[i, ] <- sapply(object$trial$dose_set, pdlt_function)
    }

    prob_dlt_dose <- apply(pdlt_dose, 2, median)
    hpd_pdlt_dose <- apply(pdlt_dose, 2, function(x) HPDinterval(as.mcmc(x)))
    tab03 <- data.frame(object$trial$dose_set, prob_dlt_dose,
                        hpd_pdlt_dose[1, ], hpd_pdlt_dose[2, ])
  } else {
    prob_dlt_dose <- NULL
    hpd_pdlt_dose <- NULL
  }

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

    if(object$trial$type == "continuous"){
      cat("P(DLT| next dose)\n")
      p02 <- data.frame(estimate = tab02[, 1],
                        hpd = paste0("(", tab02[, 2], " ; ", tab02[, 3], ")"))
      colnames(p02) <- c("Estimate", "95% HPD")
      print(p02)
      cat("\n")

    } else {
      cat("P(DLT| dose)\n")
      p03 <- data.frame(dose = tab03[, 1], estimate = round(tab03[, 2], 2),
                        hpd = paste0("(", round(tab03[, 3], 2), " ; ",
                                     round(tab03[, 4], 2), ")"))
      colnames(p03) <- c("Dose", "Estimate", "95% HPD")
      print(p03)
    }

  } else {

    out <- list(next_dose = next_dose,
                hpd_dose = hpd_dose,
                prob_dlt_next_dose = prob_dlt_next_dose,
                hpd_pdlt_next_dose = hpd_pdlt_next_dose,
                prob_dlt_dose = prob_dlt_dose,
                hpd_pdlt_dose = hpd_pdlt_dose)
    return(out)
  }

}

