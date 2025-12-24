plot_baseline_CHF <- function(b_mean, knots, data, config,
                              ngrid = 300,
                              main = "Baseline Cumulative Hazard",
                              file = NULL,
                              width = 6,
                              height = 5,
                              res = 600) {

  # 如果指定文件名，则打开图形设备
  if (!is.null(file)) {
    png(filename = file, width = width, height = height,
        units = "in", res = res)
    on.exit(dev.off())
  }
  
  # 时间范围
  t_min <- min(data$Stop, na.rm = TRUE)
  t_max <- max(data$Stop, na.rm = TRUE)
  #tgrid <- seq(t_min, t_max, length.out = ngrid)
  tgrid <- sort(unique(c(seq(knots[1], knots[length(knots)], length.out = 300), knots)))  #取每个 knot 点之间插值
  # tgrid <- seq(knots[1], knots[K+1], length.out = 300) # 只画knots区间的

  print(b_mean)
  # 真实 vs 估计
  H_true <- config$H0_true(tgrid)
  H_hat  <- M0(tgrid, b_mean, knots)

  plot(
    tgrid, H_true,
    type = "l",
    lwd  = 2,
    col  = "black",
    xlab = "t",
    ylab = expression(Lambda[0](t)),
    main = main
  )

  lines(tgrid, H_hat, lwd = 2, col = "red", lty = 2)

  abline(v = knots, col = "grey70", lty = 3)

  legend(
    "topleft",
    legend = c("True baseline CHF", "Piecewise linear estimate"),
    col    = c("black", "red"),
    lwd    = 2,
    lty    = c(1, 2),
    bty    = "n"
  )
}
