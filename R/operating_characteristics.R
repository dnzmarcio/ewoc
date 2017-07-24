#'Overdose loss function
#'
#'Calculate the overdose loss function given by \eqn{alpha(true_mtd -
#'mtd_estimate)} if \eqn{true_mtd > mtd_estimate} and \eqn{(1 - alpha)
#'(mtd_estimate - true_mtd)} if \eqn{true_mtd < mtd_estimate}.
#'
#'@param mtd_estimate a numerical value of the MTD estimate.
#'@param true_MTD a numerical value of the true MTD.
#'@param alpha Value of the loss for an underdose while the loss for
#'an overdose is (1 - \code{alpha}).
#'
#'@return value of the evaluated loss function.
#'
#'@export
overdose_loss <- function (mtd_estimate, true_mtd, alpha) {

  out <- ifelse(mtd_estimate < true_mtd, alpha*(true_mtd - mtd_estimate),
                (1 - alpha)*(mtd_estimate - true_mtd))
  return(out)
}

#'Evaluation of the DLT rate
#'
#'Calculate the DLT rate for each trial, the average DLT rate, the percent
#'of trials which have \eqn{DLT rate > target_rate + margin}, the percent
#'of trials which have \eqn{DLT rate < target_rate - margin} and the percent
#'of trials which have \eqn{target_rate - margin < DLT rate < target_rate + margin}.
#'
#'@param dlt_matrix a matrix of the number of DLT for each step of the trial (column)
#'and for each trial (row).
#'@param trial a logical value indicating if the DLT rate for each trial should be returned.
#'@param target_rate a numerical value of the target rate of DLT.
#'@param margin a numerical value of the acceptable distance from the \code{target_rate}.
#'
#'@return \code{trial} a numerical vector of the DLT rate for each trial.
#'@return \code{average} a numerical value of the average of DLT rate considering a batch of trials.
#'@return \code{upper} the percent of trials which the
#'\code{DLT rate > target_rate + margin} if \code{margin != NULL} and
#'\code{target_rate != NULL}.
#'@return \code{lower} the percent of trials which the
#'\code{DLT rate < target_rate - margin} if \code{margin != NULL} and
#'\code{target_rate != NULL}.
#'@return \code{interval} the percent of trials which the
#'\code{target_rate - margin < DLT rate < target_rate + margin} if \code{margin != NULL} and
#'\code{target_rate != NULL}.
#'
#'@export
dlt_rate <- function(dlt_matrix, trial = FALSE,
                     target_rate = NULL, margin = NULL, digits = 2) {

  dlt_matrix <- as.matrix(dlt_matrix)

  aux_upper <- function(dlt, target_rate, margin) {
    out <- ifelse(mean(dlt, na.rm = TRUE) > target_rate + margin, 1, 0)
    return(out)
  }

  aux_lower <- function(dlt, target_rate, margin) {
    out <- ifelse(mean(dlt, na.rm = TRUE) < target_rate - margin, 1, 0)
    return(out)
  }

  dlt_trial <- round(rowMeans(dlt_matrix, na.rm = TRUE), digits)
  dlt_average <- round(mean(dlt_trial, na.rm = TRUE), digits)

  if (!is.null(margin) & !is.null(target_rate)){
    dlt_upper <- apply(dlt_matrix, 1, aux_upper,
                    target_rate = target_rate, margin = margin)
    dlt_upper <- round(mean(dlt_upper, na.rm = TRUE)*100, digits)

    dlt_lower <- apply(dlt_matrix, 1, aux_lower,
                         target_rate = target_rate, margin = margin)
    dlt_lower <- round(mean(dlt_lower, na.rm = TRUE)*100, digits)

    dlt_interval <- dlt_upper + dlt_lower

    out <- list(trial = dlt_trial, average = dlt_average,
                upper = dlt_upper, lower = dlt_lower, interval = dlt_interval)

    if (trial)
      out <- list(average = dlt_average,
                  upper = dlt_upper, lower = dlt_lower, interval = dlt_interval)

  } else {
    out <- list(trial = dlt_trial, average = dlt_average)

    if (trial)
      out <- list(average = dlt_average)
  }
  return(out)
}

#'Evaluation of the stop rule
#'
#'Calculate the average, minimum, maximum number of patients to stop a trial and
#'the percent of stopped trials. Stopped trials contain NA after the last
#'assigned dose.
#'
#'@param dose Matrix of the number of DLT for each step of the trial (column)
#'and for each trial (row).
#'@return A list consisting of
#'\itemize{
#' \item{\code{average}: }{Average number of patients to stop a trial.}
#' \item{\code{min}: }{Minimum number of patients to stop a trial.}
#' \item{\code{max}: }{Maximum number of patients to stop a trial.}
#' \item{\code{nstop}: }{Percent of stopped trials}.
#'}
#'@export
stop_rule <- function(dlt_matrix, sample_size, digits = 2) {

  dlt_matrix <- as.matrix(dlt_matrix)

  index <- which(rowSums(!is.na(dlt_matrix)) < sample_size, arr.ind = TRUE)

  if(length(index) > 0) {
    result <- apply(dlt_matrix[index, ], 1, function(x) sum(!is.na(x)))
  } else {
    result <- 0
  }

  out <- list(average = mean(result),
              min = min(result),
              max = max(result),
              nstop = round(100*length(index)/
                              nrow(dlt_matrix), digits))
  return(out)
}


#'Percent of dose greater than an upper bound
#'
#'Calculate the percent of dose which are greater than \code{true_mtd + margin}.
#'
#'@param dose a numerical matrix of assigned doses for each step of the trial (column)
#'and for each trial (row).
#'@param true_MTD a numerical value of the true Maximum Tolerable Dose.
#'@param margin a numerical value of the acceptable upper margin of distance from the
#'\code{true_MTD}.
#'
#'@return \code{percent} the percent of the dose assigned which are greater than
#' \code{true_MTD + margin} for each trial.
#'@return \code{average} the average of the percent the dose assigned which are
#' greater than \code{true_MTD + margin}.
#'
#'@export
overdose_percent <- function(dose_matrix, true_mtd,  margin = NULL, digits = 2) {

  dose_matrix <- as.matrix(dose_matrix)

  aux <- function(dose, true_mtd,  margin) {
    observed_number <- sum(dose > true_mtd + margin)
    out <- round(100*observed_number/length(dose), digits)
    return(out)
  }

  percent <- apply(dose_matrix, 1, aux,
               true_mtd = true_mtd,  margin = margin)
  average <- round(mean(percent, na.rm = TRUE), digits)
  out <- list(trial = percent, average = average)

  return(out)
}

#'Percent of dose smaller than a lower bound
#'
#'Calculate the percent of dose which are smaller than \code{true_MTD - margin}.
#'
#'@param dose a numerical matrix of assigned for each step of the trial (column)
#'and for each trial (row).
#'@param true_MTD a numerical value of the true Maximum Tolerable Dose.
#'@param margin a numerical Value of the acceptable lower margin of distance from the
#'\code{true_MTD}.
#'
#'@return \code{percent} the percent of the dose assigned which are smaller than
#' \code{true_MTD - margin} for each trial.
#'@return \code{average} the average of the percent the dose assigned which are
#' greater than \code{true_MTD - margin}.
#'
#'@export
underdose_percent <- function(dose_matrix, true_mtd,  margin, digits = 2) {

  dose_matrix <- as.matrix(dose_matrix)

  aux <- function(dose, true_mtd,  margin) {
    observed_number <- sum(dose < true_mtd - margin)
    out <- round(100*observed_number/length(dose), digits)
    return(out)
  }

  percent <- apply(dose_matrix, 1, aux,
                   true_mtd = true_mtd,  margin = margin)
  average <- round(mean(percent, na.rm = TRUE), digits)
  out <- list(trial = percent, average = average)

  return(out)
}

#'Percent of optimal doses
#'
#'Calculate the percent of dose which are inside the interval \code{[true_MTD -
#'margin ; true_MTD + margin]}.
#'
#'@param dose_matrix a numerical matrix of assigned doses for each step of the trial (column)
#'and for each trial (row).
#'@param true_MTD a numerical value of the true Maximum Tolerable Dose.
#'@param margin a numerical value of the acceptable margin of distance from the
#'\code{true_MTD}.
#'
#'@return \code{trial} the percent of optimal doses for each trial.
#'@return \code{average} the average percent of optimal doses.
#'
#'@export
optimal_mtd_interval <- function(dose_matrix, true_mtd, margin, digits = 2) {

  dose_matrix <- as.matrix(dose_matrix)

  aux <- function(dose, true_mtd,  margin) {
    observed_number <- sum(dose > true_mtd - margin & dose < true_mtd + margin)
    out <- round(100*observed_number/length(dose), digits)
    return(out)
  }

  percent <- apply(dose_matrix, 1, aux,
                   true_mtd = true_mtd,  margin = margin)
  average <- round(mean(percent, na.rm = TRUE), digits)
  out <- list(trial = percent, average = average)

  return(out)
}


#'@export
optimal_toxicity_interval <- function(dose_matrix, theta, margin, pdlt, digits = 2) {

  dose_matrix <- as.matrix(dose_matrix)

  aux <- function(dose, theta,  margin) {
    prob <- pdlt(dose)
    observed_number <- sum(prob > theta - margin & prob < theta + margin)
    out <- round(100*observed_number/length(dose), digits)
    return(out)
  }

  percent <- apply(dose_matrix, 1, aux,
                   theta = theta,  margin = margin)
  average <- round(mean(percent, na.rm = TRUE), digits)
  out <- list(trial = percent, average = average)

  return(out)
}

#'Bias of the MTD estimates
#'
#'Calculate the bias.
#'
#'@param mtd_estimate a numerical vector of the MTD estimates.
#'@param true_mtd a numerical value of the true Maximum Tolerable Dose.
#'
#'@return \code{bias} bias of the MTD estimates.
#'
#'@export
mtd_bias <- function(mtd_estimate, true_mtd) {
  out <- mean(mtd_estimate - true_mtd, na.rm = TRUE)
  return(out)
}

#'Mean Square Error of the MTD estimates
#'
#'Calculate the Mean Square Error (MSE).
#'
#'@param mtd_estimate a numerical vector of the MTD estimates.
#'@param true_mtd a numerical value of the true Maximum Tolerable Dose.
#'
#'@return \code{mse} MSE of the MTD estimates.
#'
#'@export
mtd_mse <- function(mtd_estimate, true_mtd) {
  out <- mean((mtd_estimate - true_mtd)^2, na.rm = TRUE)
  return(out)
}

#'@export
accuracy_index <- function (mtd_estimate, dose_set, true_prob, theta,
                            loss = c("squared", "absolute", "classification",
                                     "overdose"), alpha = NULL, digits = 5) {

  mtd_estimate <- factor(mtd_estimate, levels = dose_set)
  estimate_prob <- prop.table(table(mtd_estimate))

  if (loss == "squared")
    dist <- (true_prob -  theta)^2
  if (loss == "absolute")
    dist <- abs(true_prob -  theta)
  if (loss == "classification")
    dist <- as.numeric(round(true_prob, digits) !=  theta)
  if (loss == "overdose") {
    if (!is.null(alpha)) {
      dist <- overdose_loss(true_prob, theta, alpha)
    } else {
      stop("loss = 'overdose' requires a value for alpha.")
    }
  }
  out <- 1 - length(dose_set)*sum(dist*estimate_prob)/sum(dist)

  return(out)
}

#'@export
average_toxicity <- function (dose, dose_set, true_prob, theta) {

  aux_toxicity <- function (x, dose_set, true_prob) {
    x <- factor(x, levels = dose_set)
    freq <- as.numeric(table(x))
    observed_average <- sum(true_prob*freq)
  }

  observed_trial <- apply(dose, 1, FUN = aux_toxicity,
                          dose_set = dose_set, true_prob = true_prob)

  expected_average <- ncol(dose)*theta

  out <- list(observed_average = mean(observed_trial),
              expected_average = expected_average)

  return(out)
}






