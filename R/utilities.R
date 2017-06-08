#'Transform a probability into logit scale
#'
#'@export
logit <- function(p) {
  out <- log(p/(1 - p))
  return(out)
}

#'@export
standard_dose <- function(dose, min_dose, max_dose) {

  out <- (dose - min_dose)/(max_dose - min_dose)
  return(out)
}

#'@export
inv_standard_dose <- function(dose, min_dose, max_dose) {

  out <- dose*(max_dose - min_dose) + min_dose
  return(out)
}

#'@export
round_down <- function(dose, grid){

  dif <- dose - grid
  index <- ifelse(length(which(dif >= 0)) == 0, 1, max(which(dif >= 0)))
  out <- grid[index]
  return(out)
}

#'@export
round_nearest <- function(dose, grid){

  dif <- abs(dose - grid)
  index <- ifelse(length(which.min(dif)) == 0, 1, which.min(dif))
  out <- grid[index]
  return(out)
}

#'@export
rounding_system <- function(dose, grid, rounding) {

  dose <- matrix(dose, ncol = 1)

  if (rounding == "down")
    out <- apply(dose, 1, round_down, grid =  grid)
  if (rounding == "nearest")
    out <- apply(dose, 1, round_nearest, grid = grid)

  return(out)

}

#'Simulation of a trial using Escalation Over with Dose Control
#'
#'Performing a simulation of several phase I clinical trial based on the
#'Escalation Over Dose Control (EWOC).
#'
#'@param step_zero an object from the classes 'ewoc_d1basic', 'ewoc_d1extended',
#''ewoc_d1ph', 'ewoc_d1multinomial', 'ewoc_d1ordinal', and 'ewoc_d1continuous'.
#'@param n_sim a number indicating the number of phase I clinical trials
#'to be simulated.
#'@param sample_size a number indicating the number of patients enrolled for
#'each clinical trial.
#'@param alpha_strategy a character indicating the strategy to apply for the
#'feasibility value. Default is "fixed". Options are "increasing" and
#'"conditional".
#'@param alpha_rate a numerical value indicating the rate of the
#'feasibility strategy. Only necessary if alpha_strategy is either
#''increasing' or 'conditional'.
#'@param response_sim a function which is self-contained and will beused
#'as a generator function of the response variables in the simulation.
#'Its only imput is 'dose' and output is in the same format of the
#'response variable used to create object 'step_zero'.
#'
#'@export
trial_simulation <- function(step0, n_sim, sample_size, alpha_strategy = "fixed",
                             alpha_rate = NULL, response_sim){
  UseMethod("trial_simulation")
}

#'@export
next_dose <- function(data) {
  UseMethod("next_dose")
}

#'@export
ewoc_jags <- function(data, n_adapt, burn_in, n_mcmc, n_thin, n_chains) {
  UseMethod("ewoc_jags")
}

#'@export
feasibility <- function(current_alpha, strategy, dlt, resolution, rate){
  if (strategy == "constant")
    next_alpha <- current_alpha
  if (strategy == "increasing")
    next_alpha <- ifelse(current_alpha > 0.49, 0.5, current_alpha + rate)
  if (strategy == "conditional")
    next_alpha <- ifelse(current_alpha > 0.49, 0.5,
                         current_alpha + sum(resultion)*rate)
  out <- next_alpha
  return(out)
}



