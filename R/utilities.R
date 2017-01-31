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

#'@export
trial_simulation <- function(object, n_sim, sample_size, alpha_strategy,
                             response_sim, pdlt_sim, rcovariate){
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
feasibility <- function(current_alpha, strategy, dlt, resolution){

  next_alpha <- current_alpha

  index <- max(which(!is.na(resolution)))

  if (strategy == "increasing")
    next_alpha <- ifelse(next_alpha > 0.49, 0.5, next_alpha + 0.05)
  if (strategy == "conditional" & is.finite(index))
    next_alpha <- ifelse(next_alpha > 0.49, 0.5,
                         ifelse(dlt[index] == 0,
                                next_alpha + 0.05, next_alpha))
  out <- next_alpha
  return(out)
}


#'@export
plot.ewoc <- function(object){

  mtd <- as.numeric(object$mtd)
  data_plot <- data.frame(x = mtd)
  dens <- density(mtd)
  shade <- with(dens, data.frame(x, y))

  label <- paste("Next dose:", round(object$next_dose, 2))

  gp <- ggplot(data_plot, aes(x = mtd)) + geom_density() +
    geom_vline(xintercept = as.numeric(object$hpd_dose), linetype = 2, size = 1.2) +
    geom_ribbon(data =
                subset(shade, x > object$hpd_dose[1] & x < object$next_dose),
                aes(ymax = y, x = x), ymin = 0, fill="red", alpha=0.3) +
    labs(y = "Density", x = "MTD") +
    annotate("text", x = object$next_dose + 0.10*object$next_dose,
             y = max(shade$y)/2,
             label = label) +
    theme_bw()

  plot(gp)
  return(gp)
}

#'@export
plot.ewoc_covariate <- function(object){

  mtd <- as.numeric(object$mtd)
  data_plot <- data.frame(x = mtd)
  dens <- density(mtd)
  shade <- with(dens, data.frame(x, y))

  label <- paste("Next dose:", round(object$next_dose, 2))

  gp <- ggplot(data_plot, aes(x = mtd)) + geom_density() +
    geom_vline(xintercept = as.numeric(object$hpd_dose), linetype = 2, size = 1.2) +
    geom_ribbon(data =
                  subset(shade, x > object$hpd_dose[1] & x < object$next_dose),
                aes(ymax = y, x = x), ymin = 0, fill="red", alpha=0.3) +
    labs(y = "Density", x = "MTD") +
    annotate("text", x = object$next_dose + 0.10*object$next_dose,
             y = max(shade$y)/2,
             label = label) +
    theme_bw()

  plot(gp)
  return(gp)
}


