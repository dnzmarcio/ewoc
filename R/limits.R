#'@export
limits_d1nocov <- function(first_dose, last_dose, min_dose, max_dose, type,
                           rounding, dose_set){

  if (is.null(min_dose) | is.null(max_dose))
    stop ("'min_dose' and 'max_dose have to be defined'")

  if (is.function(min_dose))
    stop("'min_dose' cannot be a function for this design.")

  if (is.function(max_dose))
    stop("'max_dose' cannot be a function for this design.")

  if (type == "discrete")
    if (rounding == "down")
      max_dose <- max_dose + 1

  if (min_dose > max_dose)
    stop("'min_dose' should be smaller than the 'max_dose'.")

  if (is.null(first_dose) | is.null(last_dose)) {
    if (type == "continuous"){
      first_dose <- min_dose
      last_dose <- max_dose
      warning("'first_dose' and 'last_dose' were defined as the minimum and maximum doses, respectively.")
    }
    if (type == "discrete") {
      first_dose <- dose_set[1]
      last_dose <- dose_set[length(dose_set)]
      warning("'first_dose' and  'last_dose' were defined as the first and last elements of 'dose_set', respectively.")
    }
  }

  if (is.function(first_dose))
    stop("'first_dose' cannot be a function for this design.")


  if (is.function(last_dose))
    stop("'last_dose' cannot be a function for this design.")

  out <- list(first_dose = first_dose, last_dose = last_dose,
              min_dose = min_dose, max_dose = max_dose)
  return(out)
}

#'@export
limits_d1cov <- function(first_dose, last_dose, min_dose, max_dose, type,
                          rounding, dose_set, covariate){

  if (is.null(min_dose) | is.null(max_dose))
    stop ("'min_dose' and 'max_dose' have to be defined.")

  if (!is.function(min_dose)){
    min_dose_value <- min_dose
    min_dose <- function(covariate = NULL){
      out <- min_dose_value
      return(out)
    }

    min_dose <- Vectorize(min_dose)
  }

  if (!is.function(max_dose)) {
    max_dose_value <- max_dose

    if (type == "discrete" & rounding == "down"){
      max_dose <- function(covariate = NULL){
        out <- max_dose_value + 1
        return(out)
      }
    } else {
      max_dose <- function(covariate = NULL){
        out <- max_dose_value
        return(out)
      }
    }

    max_dose <- Vectorize(max_dose)
  }

  if (any(min_dose(covariate) > max_dose(covariate)))
    stop("'min_dose' should be smaller than the 'max_dose'.")

  if (is.null(first_dose) | is.null(last_dose)) {
    if (type == "continuous"){
      first_dose <- min_dose
      last_dose <- max_dose
      warning("'first_dose' and 'last_dose' were defined as the minimum and maximum doses, respectively.")
    }
    if (type == "discrete") {
      first_dose <- dose_set[1]
      last_dose <- dose_set[length(dose_set)]
      warning("'first_dose' and  'last_dose' were defined as the first and last element of 'dose_set', respectively.")
    }
  }

  if (any(first_dose(covariate) > last_dose(covariate)))
    stop("'first_dose' should be smaller than the 'last_dose'.")

  out <- list(first_dose = first_dose, last_dose = last_dose,
              min_dose = min_dose, max_dose = max_dose)
  return(out)
}


