#'Generating a stop rule function
#
#'@param step an object from the classes 'ewoc_d1classic', 'ewoc_d1extended',
#''ewoc_d1ph' generated during the simulation.
#'
#'@details The stop rule function is evaluated at each step of the trial.
#'It can defined based on any information contained in the object 'step' that
#'is the output from one of the functions 'ewoc_d1classic',
#''ewoc_d1extended', 'ewoc_d1ph'.
#'
#'@return a logical character indicating if the trial should be stopped or not.
#'
#'@examples
#'\dontshow{
#'DLT <- 0
#'dose <- 30
#'step_zero <- ewoc_d1classic(DLT ~ dose, type = 'discrete',
#'                            theta = 0.33, alpha = 0.25,
#'                            min_dose = 0, max_dose = 100,
#'                            dose_set = seq(0, 100, 20),
#'                            rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                            mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                            rounding = "nearest")
#'stop_rule_sim(step_zero)
#'response_sim <- response_d1classic(rho = 0.05, mtd = 20, theta = 0.33,
#'                                   min_dose = 10, max_dose = 50)
#'sim <- ewoc_simulation(step_zero = step_zero,
#'                        n_sim = 1, sample_size = 2,
#'                        alpha_strategy = "increasing",
#'                        response_sim = response_sim,
#'                        stop_rule_sim = stop_rule_sim,
#'                        ncores = 2)
#'}
#'
#'\dontrun{
#'DLT <- 0
#'dose <- 30
#'step_zero <- ewoc_d1classic(DLT ~ dose, type = 'discrete',
#'                            theta = 0.33, alpha = 0.25,
#'                            min_dose = 0, max_dose = 100,
#'                            dose_set = seq(0, 100, 20),
#'                            rho_prior = matrix(1, ncol = 2, nrow = 1),
#'                            mtd_prior = matrix(1, ncol = 2, nrow = 1),
#'                            rounding = "nearest")
#'stop_rule_sim(step_zero)
#'response_sim <- response_d1classic(rho = 0.05, mtd = 20, theta = 0.33,
#'                                   min_dose = 10, max_dose = 50)
#'sim <- ewoc_simulation(step_zero = step_zero,
#'                        n_sim = 1, sample_size = 30,
#'                        alpha_strategy = "increasing",
#'                        response_sim = response_sim,
#'                        stop_rule_sim = stop_rule_sim,
#'                        ncores = 2)
#'}
#'
#'@export
stop_rule_sim <- function(step){
  out <- FALSE
}
