#'Generating a stop rule function
#
#'@param step an object from the classes 'ewoc_d1classic', 'ewoc_d1extended',
#''ewoc_d1ph' generated during the simulation.
#'
#'@details The stop rule function is evaluated at each step of the trial.
#'It can defined based on any information contained in the object 'step'.
#'The object 'step' is the output from one of the functions 'ewoc_d1classic',
#''ewoc_d1extended', 'ewoc_d1ph'. In this way, a dummy step can be created with
#'a specific function in order to explore the information available.
#'
#'@return a logical character indicating if the trial should be stopped or not.
#'
#'@export
stop_rule_sim <- function(step){
  out <- FALSE
}
