logit <- function(p) {
  out <- log(p/(1 - p))
  return(out)
}

comb <- function(...) {
  mapply('rbind', ..., SIMPLIFY=FALSE)
}


standard_dose <- function(dose, min_dose, max_dose) {

  out <- (dose - min_dose)/(max_dose - min_dose)
  return(out)
}


inv_standard_dose <- function(dose, min_dose, max_dose) {

  out <- dose*(max_dose - min_dose) + min_dose
  return(out)
}

round_down <- function(dose, grid){

  dif <- dose - grid
  index <- ifelse(length(which(dif >= 0)) == 0, 1, max(which(dif >= 0)))
  out <- grid[index]
  return(out)
}

round_nearest <- function(dose, grid){

  dif <- abs(dose - grid)
  index <- ifelse(length(which.min(dif)) == 0, 1, which.min(dif))
  out <- grid[index]
  return(out)
}

rounding_system <- function(dose, grid, rounding) {

  dose <- matrix(dose, ncol = 1)

  if (rounding == "down")
    out <- apply(dose, 1, round_down, grid =  grid)
  if (rounding == "nearest")
    out <- apply(dose, 1, round_nearest, grid = grid)

  return(out)

}

next_dose <- function(data) {
  UseMethod("next_dose")
}

ewoc_jags <- function(data, n_adapt, burn_in, n_mcmc, n_thin, n_chains) {
  UseMethod("ewoc_jags")
}

feasibility <- function(alpha, strategy, rate, dlt, resolution){
  if (strategy == "constant")
    next_alpha <- alpha[1]
  if (strategy == "increasing")
    next_alpha <- ifelse(round(alpha[length(alpha)], 3) >= 0.5, 0.5,
                         round((alpha[length(alpha)] + rate), 3))
  if (strategy == "conditional")
    next_alpha <-
      ifelse(round(alpha[length(alpha)], 3) >= 0.5, 0.5,
             ifelse(round((alpha[1] + sum(resolution*(1 - dlt))*rate), 1) >= 0.5,
                    0.5, round((alpha[1] + sum(resolution*(1 - dlt))*rate), 3)))

  out <- next_alpha
  return(out)
}

