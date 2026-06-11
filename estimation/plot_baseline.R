plot_baseline_CHF <- function(xi_mean, sp, data, config,
                              ngrid = 300,
                              main = "Baseline Cumulative Hazard",
                              file = NULL,
                              width = 6,
                              height = 5,
                              res = 600) {

  if (!is.null(file)) {
    png(filename = file, width = width, height = height,
        units = "in", res = res)
    on.exit(dev.off())
  }

  t_min <- min(data$Stop, na.rm = TRUE)
  t_max <- max(data$Stop, na.rm = TRUE)
  tgrid <- seq(t_min, t_max, length.out = ngrid)

  H_true <- config$H0_true(tgrid)
  H_hat  <- M0(tgrid, xi_mean, sp)

  plot(tgrid, H_true, type = "l", lwd = 2, col = "black",
       xlab = "t", ylab = expression(Lambda[0](t)), main = main)
  lines(tgrid, H_hat, lwd = 2, col = "red", lty = 2)

  legend("topleft",
         legend = c("True baseline CHF", "M/I-spline estimate"),
         col = c("black", "red"), lwd = 2, lty = c(1, 2), bty = "n")
}
