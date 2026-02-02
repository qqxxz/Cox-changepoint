library(parallel)
ncore <- detectCores() - 1  # 留一个核给系统

################################无变点模型（γ=0）的拟合######################
fit_null_model <- function(data, knots, p) {

  K <- length(knots) - 1

  neg_loglik_null <- function(par, data_input, knots, p) {

    beta  <- par[1:p]
    theta_b <- par[(p+1):(p+K)]
    b <- exp(theta_b)

    X <- as.matrix(data_input[, paste0("X", 1:p)])
    Z  <- data_input$Stop
    T0 <- data_input$Start
    delta <- data_input$Event

    psi <- as.vector(X %*% beta)
    psi <- pmin(pmax(psi, -20), 20)

    mZ <- m0(Z, b, knots)
    MZ <- M0(Z, b, knots)
    MT <- M0(T0, b, knots)

    ll <- sum(delta * (log(mZ) + psi)) -
          sum((MZ - MT) * exp(psi))

    return(-ll)
  }

  init <- c(rep(1, p), rep(0.1, K))

  fit <- optim(
    par = init,
    fn  = neg_loglik_null,
    data_input = data,
    knots = knots,
    p = p,
    method = "L-BFGS-B",
    lower = c(rep(-5, p), rep(log(1e-2), K)),
    upper = c(rep(5, p), rep(log(10), K))
  )

  beta_hat <- fit$par[1:p]
  b_hat    <- exp(fit$par[(p+1):(p+K)])

  list(
    beta = beta_hat,
    b    = b_hat,
    logLik = -fit$value
  )
}


##########################SUP 统计量######################
compute_SUP_stat_score <- function(data, fit0, knots, eta_grid, p) {

  beta_hat <- fit0$beta
  b_hat    <- fit0$b

  X  <- as.matrix(data[, paste0("X", 1:p)])
  Z  <- data$Stop
  T0 <- data$Start
  delta <- data$Event
  x1 <- data$X1

  MZ <- M0(Z,  b_hat, knots)
  MT <- M0(T0, b_hat, knots)

  Xt <- if (p > 1) {
    cbind(data[, paste0("X", 2:p)], 1)
  } else {
    matrix(1, nrow = nrow(data), ncol = 1)
  }

  psi <- as.vector(X %*% beta_hat)
  mu  <- as.numeric((MZ - MT) * exp(psi))
  r   <- as.numeric(delta - mu)

  SUP_vals <- numeric(length(eta_grid))

  for (k in seq_along(eta_grid)) {

    eta <- eta_grid[k]
    Ieta <- as.numeric(x1 > eta)

    U <- colSums(sweep(Xt, 1, r * Ieta, "*"))

    Sigma <- matrix(0, ncol(Xt), ncol(Xt))
    idx <- which(Ieta == 1)
    for (i in idx) {
    xi <- matrix(Xt[i, ], ncol = 1)
    Sigma <- Sigma + (r[i]^2) * (xi %*% t(xi))
    }
    Sigma <- Sigma + diag(1e-6, ncol(Sigma))

    SUP_vals[k] <- as.numeric(t(U) %*% solve(Sigma, U))
  }

  list(
    SUP = max(SUP_vals, na.rm = TRUE),
    SUP_path = SUP_vals
  )
}


#####################置换检验##########################3
single_cp_test_score <- function(data, config) {

  p      <- config$p
  B_perm <- config$B_perm
  k      <- config$k
  trim   <- config$eta.trim

  # 对一份样本数据做一次变点存在性检验
  knots <- select_knots(data, p)

  x1 <- data$X1
  x1_min <- min(x1, na.rm = TRUE)
  x1_max <- max(x1, na.rm = TRUE)

  # 为避免极端点数值不稳定，可在两端稍作截断
  lower <- x1_min + trim * (x1_max - x1_min)
  upper <- x1_max - trim * (x1_max - x1_min)

  eta_grid <- seq(lower, upper, length.out = k)

  fit0 <- fit_null_model(data, knots, p)

  sup_obs <- compute_SUP_stat_score(
    data, fit0, knots, eta_grid, p
  )$SUP

  sup_perm <- numeric(B_perm)

  # 构造置换分布
  for (b in 1:B_perm) {
    data_perm <- data
    data_perm$X1 <- sample(data$X1) # 打乱变点协变量

    # 重新拟合零假设模型
    fit0_perm <- fit_null_model(data_perm, knots, p)

    sup_perm[b] <- compute_SUP_stat_score(
      data_perm, fit0_perm, knots, eta_grid, p
    )$SUP
  }

  crit <- quantile(sup_perm, probs = 1 - config$alpha)

  list(
    SUP_obs = sup_obs, # 一个数
    sup_perm = sup_perm, # 一列数（分布）
    crit = crit # 一个分位数
  )
}


#################  MC 主函数#####################
MC_changepoint_test <- function(config) {

  B_mc <- config$B
  alpha_vec <- config$alpha

  one_run <- function(b) {

    # Type I / Power 不互相污染
    gamma_tag <- sum(abs(config$Gamma))  # 0 for null, >0 for power

    # ===== 每次 MC replication 单独设随机种子 =====
    mc_seed <- config$base_seed +  # 全局起点
              config$seed_multiplier * b +  # mc循环不同随机种子
              987 * config$test_id +  # type and power 不同随机种子
              654 * config$exp_id # 不同实验不同随机种子

    set.seed(mc_seed)

    cat("MC replication:", b, " seed =", mc_seed, "\n")

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
    )$Data

    test_res <- single_cp_test_score(dat, config)

    SUP_obs <- test_res$SUP_obs
    pval    <- mean(test_res$sup_perm >= SUP_obs)
    reject  <- as.numeric(SUP_obs > test_res$crit)

    list(
      SUP = SUP_obs,
      pval = pval,
      reject = reject
    )
  }

  res_list <- mclapply(1:B_mc, one_run, mc.cores = ncore)

  SUP_obs_vec <- sapply(res_list, `[[`, "SUP")
  pval_vec    <- sapply(res_list, `[[`, "pval")
  reject_mat  <- do.call(rbind, lapply(res_list, `[[`, "reject"))

  Rejection_Rate <- colMeans(reject_mat)
  names(Rejection_Rate) <- paste0("alpha=", alpha_vec)

  list(
    SUP = SUP_obs_vec,
    p_value = pval_vec,
    Rejection_Rate = Rejection_Rate
  )
}

############# MC模拟函数（封装Type I & Power）############################
run_MC_test <- function(config, gamma_vec, out_dir="results") {
  
  config$Gamma <- gamma_vec
  MC_res <- MC_changepoint_test(config)
  Rejection <- MC_res$Rejection_Rate
  
  # -------- Summary --------
  df_summary <- data.frame(
    alpha = names(Rejection),
    Rejection_Rate = as.numeric(Rejection)
  )
  
  # -------- 每次MC的SUP和p值 --------
  df_pval <- data.frame(
    Iteration = seq_along(MC_res$p_value),
    SUP = MC_res$SUP,
    p_value = MC_res$p_value
  )
  
  # -------- 实验配置 --------
    df_config <- data.frame(
    Parameter = names(config),
    Value = sapply(config, function(x) {
        if (length(x) > 1) {
        paste(x, collapse = ", ")
        } else {
        as.character(x)
        }
    }),
    stringsAsFactors = FALSE
    )
  
  # -------- Excel文件名 --------
  fname_xlsx <- paste0(
    out_dir, "/MC_Dist=", config$Distribution,
    "_n=", config$n,
    "_LT=", config$truncation,
    "_C=", config$censor,
    "_Gamma=", paste(round(gamma_vec,2), collapse="_"),
    "_B=", config$B, ".xlsx"
  )
  
  # -------- 写入Excel（三个sheet）--------
  wb <- createWorkbook()
  
  addWorksheet(wb, "Summary")
  writeData(wb, "Summary", df_summary)
  
  addWorksheet(wb, "pval_SUP")
  writeData(wb, "pval_SUP", df_pval)
  
  addWorksheet(wb, "Config")
  writeData(wb, "Config", df_config)
  
  saveWorkbook(wb, fname_xlsx, overwrite = TRUE)
  cat("Results saved to Excel:", fname_xlsx, "\n")
  
  return(list(
    MC_res = MC_res,
    file_XLSX = fname_xlsx
  ))
}
