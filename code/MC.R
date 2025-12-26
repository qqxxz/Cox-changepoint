get_true_baseline <- function(config) {

  if (config$Distribution == "Exp") {

    S0 <- function(t) exp(-config$Lambda * t)
    H0 <- function(t) config$Lambda * t

  } else if (config$Distribution %in% c("WD", "WI")) {

    # Weibull (both increasing/decreasing share same H0 form)
    S0 <- function(t) exp(-config$Lambda * t^config$V)
    H0 <- function(t) config$Lambda * t^config$V

  } else if (config$Distribution == "Gtz") {

    # Gompertz: H0(t) = (Lambda/Alpha) * (exp(Alpha * t) - 1)
    Alpha <- if(!is.null(config$Alpha)) config$Alpha else 0.1
    Lambda <- if(!is.null(config$Lambda)) config$Lambda else 0.3
    H0 <- function(t) (Lambda / Alpha) * (exp(Alpha * t) - 1)
    S0 <- function(t) exp(-H0(t))

  } else if (config$Distribution == "Quadratic") {

    a <- config$Coeff$a
    H0 <- function(t) a * t^2
    S0 <- function(t) exp(-H0(t))

  } else if (config$Distribution == "QuadraticLinear") {

    a <- config$Coeff$a
    b <- config$Coeff$b
    H0 <- function(t) a * t^2 + b * t
    S0 <- function(t) exp(-H0(t))

  } else if (config$Distribution == "PiecewiseWeibull") {

    # 使用 config 中的分段参数
    alpha1 <- config$Coeff$alpha1
    lambda1 <- config$Coeff$lambda1
    alpha2 <- config$Coeff$alpha2
    lambda2 <- config$Coeff$lambda2
    t_star <- config$Coeff$t_star

    H0 <- function(t) {
      sapply(t, function(tt) {
        if (tt <= t_star) {
          lambda1 * tt^alpha1
        } else {
          lambda1 * t_star^alpha1 + lambda2 * (tt^alpha2 - t_star^alpha2)
        }
      })
    }
    S0 <- function(t) exp(-H0(t))

  } else {
    stop("Unknown distribution")
  }

  list(S0_true = S0, H0_true = H0)
}

##################### Monte Carlo 仿真 #####################
run_simulation <- function(config) {
  true_base <- get_true_baseline(config) # 生成真实基准基线函数
  config$S0_true <- true_base$S0_true
  config$H0_true <- true_base$H0_true

  B <- config$B
  par_mat <- NULL # 初始化  
  AIE_SF_vec  <- numeric(B)
  AIE_CHF_vec <- numeric(B)
  b_mat <- NULL

  for (b in 1:B) {
    cat("Replication:", b, "\n")

    dat <- TimeindepLTRC_gnrt_ChangepointPH(
      N = config$n,
      Distribution = config$Distribution,
      eta = config$eta,
      Beta = config$Beta,
      Gamma = config$Gamma,
      truncation.percent = config$truncation,
      censor.percent = config$censor,
      x1.mean = config$x1.mean,
      x1.sd   = config$x1.sd,
      x2.mean = config$x2.mean,
      x2.sd   = config$x2.sd,
      adjust.censor = TRUE
    )$Data # 生成数据

    fit <- fit_piecewise(dat, p = config$p) # 模型拟合
    theta_hat <- fit$par

    par_mat <- rbind(par_mat, theta_hat)  # 储存参数估计
    b_mat <- rbind(b_mat, fit$b)

    tgrid <- seq(min(dat$Stop), max(dat$Stop), length.out = 200) # 构造时间网格

    AIE_SF_vec[b]  <- compute_AIE_SF(fit, dat, config)
    AIE_CHF_vec[b] <- compute_AIE_CHF(fit, dat, config)
  }

  b_mean <- colMeans(b_mat, na.rm = TRUE)
  plot_baseline_CHF(
    b_mean  = b_mean,
    knots   = fit$knots,
    data    = dat,
    config  = config,
    file    = "baseline_CHF.png"
  )

  list(
    par_mat = par_mat,
    AIE_SF  = AIE_SF_vec,
    AIE_CHF = AIE_CHF_vec,
    AIE_SF_MC  = mean(AIE_SF_vec,  na.rm = TRUE),
    AIE_CHF_MC = mean(AIE_CHF_vec, na.rm = TRUE)
  )
}

##################### 仿真统计量 #####################
summary_MC <- function(par_mat, true_par) {
    stopifnot(
    ncol(par_mat) == length(true_par)
  ) # 确保维度一致

  est_mean <- colMeans(par_mat)
  bias <- est_mean - true_par
  sse <- colMeans((par_mat - matrix(true_par,
                                    nrow(par_mat),
                                    length(true_par),
                                    byrow = TRUE))^2) # 均方误差
 
  see <- apply(par_mat, 2, sd) # 标准误

  # 95% / 99% CP 覆盖率
  CP95 <- CP99 <- rep(NA, ncol(par_mat))

  for (j in 1:ncol(par_mat)) {

    CP95[j] <- mean(
      abs(par_mat[, j] - true_par[j]) <= qnorm(0.975) * see[j],
      na.rm = TRUE
    )

    CP99[j] <- mean(
      abs(par_mat[, j] - true_par[j]) <= qnorm(0.995) * see[j],
      na.rm = TRUE
    )
  }



  return(data.frame(
    Estimate = est_mean,
    Bias = bias,
    SSE = sse,
    SEE = see,
    CP95 = CP95,
    CP99 = CP99
  ))
}

