#'@export
plot.ewoc_d1basic <- function(x, ...){

  object <- x
  sm <- summary(object, print = FALSE)

  mtd <- as.numeric(object$mtd)
  data_plot <- data.frame(x = mtd)
  dens <- density(mtd)
  shade <- with(dens, data.frame(x, y))

  label <- paste("Next dose:", round(sm$next_dose, 2))

  gp <- ggplot(data_plot, aes(x = mtd)) + geom_density() +
    geom_vline(xintercept = as.numeric(sm$hpd_dose),
               linetype = 2, size = 1.2) +
    geom_ribbon(data =
                subset(shade, x > sm$hpd_dose[1] & x < sm$next_dose),
                aes(ymax = y, x = x), ymin = 0, fill="red", alpha=0.3) +
    labs(y = "Density", x = "MTD") +
    annotate("text", x = sm$next_dose + 0.10*sm$next_dose,
             y = max(shade$y)/2,
             label = label) +
    theme_bw()

  plot(gp)
  return(gp)
}
