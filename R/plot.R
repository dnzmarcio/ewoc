#'@export
plot.ewoc_d1basic <- function(x, ...){

  object <- x
  sm <- summary(object, print = FALSE)

  mtd <- as.numeric(object$mtd)
  data_plot <- data.frame(mtd)

  dens <- density(mtd)
  shade <- with(dens, data.frame(x, y))

  label <- paste("Next dose:", round(sm$next_dose, 2))

  gp <- ggplot(data_plot, aes(x = mtd)) + geom_density() +
    geom_vline(xintercept = as.numeric(sm$next_dose),
               linetype = 2, size = 1.2) +
    geom_ribbon(data =
                subset(shade, x > sm$hpd_dose[1] & x < sm$hpd_dose[2]),
                aes(ymax = y, x = x), ymin = 0, fill="red", alpha=0.3) +
    labs(y = "Density", x = "MTD") +
    annotate("text", x = sm$next_dose + 0.10*sm$next_dose,
             y = max(shade$y)/2,
             label = label) +
    theme_bw()

  plot(gp)
  return(gp)
}

#'@export
plot.ewoc_d1extended <- function(x, ...){

  object <- x
  sm <- summary(object, print = FALSE)

  mtd <- as.numeric(object$mtd)
  data_plot <- data.frame(x = mtd)

  dens <- density(mtd, n = 2^15)
  shade <- with(dens, data.frame(x, y))

  label <- paste("Next dose:", round(sm$next_dose, 2))

  gp <- ggplot(data_plot, aes(x = mtd)) + geom_density() +
    geom_vline(xintercept = as.numeric(sm$next_dose),
               linetype = 2, size = 1.2) +
    geom_ribbon(data =
                subset(shade, x > max(sm$hpd_dose[1], object$trial$min_dose()) &
                         x < min(sm$hpd_dose[2], object$trial$max_dose())),
                aes(ymax = y, x = x), ymin = 0, fill = "red", alpha = 0.3) +
    labs(y = "Density", x = "MTD") +
    annotate("text", x = sm$next_dose + 0.30*sm$next_dose,
             y = max(shade$y)/2,
             label = label) +
    theme_bw()

  plot(gp)
  return(gp)
}

#'@export
plot.ewoc_d1ph <- function(x, ...){

  object <- x
  sm <- summary(object, print = FALSE)

  mtd <- as.numeric(object$mtd)
  data_plot <- data.frame(mtd)

  dens <- density(mtd)
  shade <- with(dens, data.frame(x, y))

  label <- paste("Next dose:", round(sm$next_dose, 2))

  gp <- ggplot(data_plot, aes(x = mtd)) + geom_density() +
    geom_vline(xintercept = as.numeric(next_dose),
               linetype = 2, size = 1.2) +
    geom_ribbon(data =
                  subset(shade, x > sm$hpd_dose[1] & x < sm$hpd_dose[2]),
                aes(ymax = y, x = x), ymin = 0, fill="red", alpha=0.3) +
    labs(y = "Density", x = "MTD") +
    annotate("text", x = sm$next_dose + 0.10*sm$next_dose,
             y = max(shade$y)/2,
             label = label) +
    theme_bw()

  plot(gp)
  return(gp)
}

#'@export
plot.ewoc_d1multinomial <- function(x, ..., next_covariable = NULL){

  object <- x

  if (is.null(next_covariable))
    next_covariable <- object$trial$levels_cov
  index <- which(next_covariable == object$trial$levels_cov)

  sm <- summary(object, print = FALSE, next_covariable = next_covariable)

  gp <- list()

  for(i in index){
    mtd <- as.numeric(object$mtd[, i])
    data_plot <- data.frame(mtd)

    dens <- density(mtd)
    shade <- with(dens, data.frame(x, y))

    label <- paste("Next dose:", round(sm$next_dose[i], 2))

    gp[[i]] <- ggplot(data_plot, aes(x = mtd)) + geom_density() +
      geom_vline(xintercept = as.numeric(sm$next_dose[i]),
                 linetype = 2, size = 1.2) +
      geom_ribbon(data = subset(shade, x > sm$hpd_dose[i, 1] &
                                  x < sm$hpd_dose[i, 2]),
                  aes(ymax = y, x = x), ymin = 0, fill = "red", alpha = 0.3) +
      labs(y = "Density", x = "MTD",
           title = paste("Group:", next_covariable[i])) +
      annotate("text", x = sm$next_dose[i] + 0.10*sm$next_dose[i],
               y = max(shade$y)/2,
               label = label) +
      theme_bw()

    plot(gp[[i]])
    invisible(readline(prompt="Press [enter] to continue"))
  }

  return(gp)
}

#'@export
plot.ewoc_d1ordinal <- function(x, ..., next_covariable = NULL){

  object <- x

  if (is.null(next_covariable))
    next_covariable <- object$trial$levels_cov
  index <- which(next_covariable == object$trial$levels_cov)

  sm <- summary(object, print = FALSE, next_covariable = next_covariable)

  gp <- list()

  for(i in index){
    mtd <- as.numeric(object$mtd[, i])
    data_plot <- data.frame(mtd)

    dens <- density(mtd)
    shade <- with(dens, data.frame(x, y))

    label <- paste("Next dose:", round(sm$next_dose[i], 2))

    gp[[i]] <- ggplot(data_plot, aes(x = mtd)) + geom_density() +
      geom_vline(xintercept = as.numeric(sm$next_dose[i]),
                 linetype = 2, size = 1.2) +
      geom_ribbon(data = subset(shade, x > sm$hpd_dose[i, 1] &
                                  x < sm$hpd_dose[i, 2]),
                  aes(ymax = y, x = x), ymin = 0, fill = "red", alpha = 0.3) +
      labs(y = "Density", x = "MTD",
           title = paste("Group:", next_covariable[i])) +
      annotate("text", x = sm$next_dose[i] + 0.10*sm$next_dose[i],
               y = max(shade$y)/2,
               label = label) +
      theme_bw()

    plot(gp[[i]])
    invisible(readline(prompt="Press [enter] to continue"))
  }

  return(gp)
}
