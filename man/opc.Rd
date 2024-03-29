% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/operating_characteristics.R
\name{opc}
\alias{opc}
\title{Operating characteristics for EWOC simulations}
\usage{
opc(sim_list, pdlt_list, mtd_list, toxicity_margin = NULL, mtd_margin = NULL)
}
\arguments{
\item{sim_list}{a list of 'ewoc_simulation' objects for different scenarios
created using the \code{\link[ewoc]{ewoc_simulation}} function.}

\item{pdlt_list}{a list of functions to calculate the probability of toxicity with a numeric vector
of doses as input and a numeric vector of probabilities as output.}

\item{mtd_list}{a list of numerical values indicating the true MTD for each scenario.}

\item{toxicity_margin}{a numerical value of the acceptable margin of distance from the
\code{target_rate}.}

\item{mtd_margin}{a numerical value of the acceptable margin of distance from the
\code{mtd_list}.}
}
\value{
\code{dlt_rate} See \code{\link[ewoc]{dlt_rate}}.

\code{dose_toxicity} See \code{\link[ewoc]{optimal_toxicity}}.

\code{mtd_toxicity} See \code{\link[ewoc]{optimal_toxicity}}.

\code{bias_mse} See \code{\link[ewoc]{mtd_bias}} and \code{\link[ewoc]{mtd_mse}}.

\code{dose_efficiency} See \code{\link[ewoc]{optimal_mtd}}.

\code{mtd_efficiency} See \code{\link[ewoc]{optimal_mtd}}.

\code{stop} See \code{\link[ewoc]{stop_rule}}.
}
\description{
Generic operating characteristics for one or more scenarios in EWOC simulations.
}
\examples{
\dontshow{
### Only one simulation
DLT <- 0
dose <- 20
step_zero <- ewoc_d1classical(DLT ~ dose, type = 'discrete',
                          theta = 0.33, alpha = 0.25,
                          min_dose = 20, max_dose = 100,
                          dose_set = seq(20, 100, 20),
                          rho_prior = matrix(1, ncol = 2, nrow = 1),
                          mtd_prior = matrix(1, ncol = 2, nrow = 1),
                          rounding = "nearest")
response_sim <- response_d1classical(rho = 0.05, mtd = 60, theta = 0.33,
                                 min_dose = 20, max_dose = 100)
sim <- ewoc_simulation(step_zero = step_zero,
                     n_sim = 1, sample_size = 2, n_cohort = 1,
                     alpha_strategy = "conditional",
                     response_sim = response_sim,
                     fixed_first_cohort =  TRUE,
                     ncores = 1)

pdlt <- pdlt_d1classical(rho = 0.05, mtd = 60, theta = 0.33,
                      min_dose = 20, max_dose = 100)

opc(sim_list = list(sim), pdlt_list = list(pdlt),
   mtd_list = list(60), toxicity_margin = 0.05, mtd_margin = 0.10)

### Two or more simulations

sim_list <- list()
mtd_list <- list()
pdlt_list <- list()

DLT <- 0
dose <- 20
step_zero <- ewoc_d1classical(DLT ~ dose, type = 'discrete',
                          theta = 0.33, alpha = 0.25,
                          min_dose = 20, max_dose = 100,
                          dose_set = seq(20, 100, 20),
                          rho_prior = matrix(1, ncol = 2, nrow = 1),
                          mtd_prior = matrix(1, ncol = 2, nrow = 1),
                          rounding = "nearest")
mtd_list[[1]] <- 60
response_sim <- response_d1classical(rho = 0.05, mtd = mtd_list[[1]],
                                 theta = 0.33,
                                 min_dose = 20, max_dose = 100)
sim_list[[1]] <- ewoc_simulation(step_zero = step_zero,
                     n_sim = 1, sample_size = 2, n_cohort = 1,
                     alpha_strategy = "conditional",
                     response_sim = response_sim,
                     fixed_first_cohort =  TRUE,
                     ncores = 1)
pdlt_list[[1]] <- pdlt_d1classical(rho = 0.05, mtd = mtd_list[[1]],
                               theta = 0.33,
                               min_dose = 20, max_dose = 100)

mtd_list[[2]] <- 40
response_sim <- response_d1classical(rho = 0.05, mtd = mtd_list[[2]],
                                 theta = 0.33,
                                 min_dose = 20, max_dose = 100)
sim_list[[2]] <- ewoc_simulation(step_zero = step_zero,
                     n_sim = 1, sample_size = 2, n_cohort = 1,
                     alpha_strategy = "conditional",
                     response_sim = response_sim,
                     fixed_first_cohort =  TRUE,
                     ncores = 1)

pdlt_list[[2]] <- pdlt_d1classical(rho = 0.05, mtd = mtd_list[[2]],
                               theta = 0.33,
                               min_dose = 20, max_dose = 100)

opc(sim_list = sim_list, pdlt_list = pdlt_list,
   mtd_list = mtd_list, toxicity_margin = 0.05, mtd_margin = 0.10)
}


\dontrun{
### Only one simulation
DLT <- 0
dose <- 20
step_zero <- ewoc_d1classical(DLT ~ dose, type = 'discrete',
                          theta = 0.33, alpha = 0.25,
                          min_dose = 20, max_dose = 100,
                          dose_set = seq(20, 100, 20),
                          rho_prior = matrix(1, ncol = 2, nrow = 1),
                          mtd_prior = matrix(1, ncol = 2, nrow = 1),
                          rounding = "nearest")
response_sim <- response_d1classical(rho = 0.05, mtd = 60, theta = 0.33,
                                 min_dose = 20, max_dose = 100)
sim <- ewoc_simulation(step_zero = step_zero,
                     n_sim = 1, sample_size = 30, n_cohort = 1,
                     alpha_strategy = "conditional",
                     response_sim = response_sim,
                     fixed_first_cohort =  TRUE,
                     ncores = 1)

pdlt <- pdlt_d1classical(rho = 0.05, mtd = 60, theta = 0.33,
                      min_dose = 20, max_dose = 100)

opc(sim_list = list(sim), pdlt_list = list(pdlt),
   mtd_list = list(60), toxicity_margin = 0.05, mtd_margin = 6)

### Two or more simulations

sim_list <- list()
mtd_list <- list()
pdlt_list <- list()

DLT <- 0
dose <- 20
step_zero <- ewoc_d1classical(DLT ~ dose, type = 'discrete',
                          theta = 0.33, alpha = 0.25,
                          min_dose = 20, max_dose = 100,
                          dose_set = seq(20, 100, 20),
                          rho_prior = matrix(1, ncol = 2, nrow = 1),
                          mtd_prior = matrix(1, ncol = 2, nrow = 1),
                          rounding = "nearest")
mtd_list[[1]] <- 60
response_sim <- response_d1classical(rho = 0.05, mtd = mtd_list[[1]],
                                 theta = 0.33,
                                 min_dose = 20, max_dose = 100)
sim_list[[1]] <- ewoc_simulation(step_zero = step_zero,
                     n_sim = 1, sample_size = 30, n_cohort = 1,
                     alpha_strategy = "conditional",
                     response_sim = response_sim,
                     fixed_first_cohort =  TRUE,
                     ncores = 1)
pdlt_list[[1]] <- pdlt_d1classical(rho = 0.05, mtd = mtd_list[[1]],
                               theta = 0.33,
                               min_dose = 20, max_dose = 100)

mtd_list[[2]] <- 40
response_sim <- response_d1classical(rho = 0.05, mtd = mtd_list[[2]],
                                 theta = 0.33,
                                 min_dose = 20, max_dose = 100)
sim_list[[2]] <- ewoc_simulation(step_zero = step_zero,
                     n_sim = 1, sample_size = 30, n_cohort = 1,
                     alpha_strategy = "conditional",
                     response_sim = response_sim,
                     fixed_first_cohort =  TRUE,
                     ncores = 1)

pdlt_list[[2]] <- pdlt_d1classical(rho = 0.05, mtd = mtd_list[[2]],
                               theta = 0.33,
                               min_dose = 20, max_dose = 100)

opc(sim_list = sim_list, pdlt_list = pdlt_list,
   mtd_list = mtd_list, toxicity_margin = 0.05, mtd_margin = 6)
}

}
\references{
Diniz, M. A., Tighiouart, M., & Rogatko, A. (2019). Comparison between continuous and discrete doses for model based designs in cancer dose finding. PloS one, 14(1).
}
