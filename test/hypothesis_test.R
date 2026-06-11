library(parallel)
ncore <- detectCores() - 1  # 留一个核给系统

################################无变点模型（γ=0）的拟合######################
fit_null_model <- function(data, sp, p) {
  sp <- spline_setup(data$Stop, sp$k, sp$od, knotc = sp$knotc,
                     times_T0 = data$Start)

  neg_loglik_null <- function(par, data_input, sp, p) {
    beta     <- par[1:p]
    theta_xi <- par[(p + 1):(p + sp$p)]
    xi       <- theta_xi^2

    X <- as.matrix(data_input[, paste0("X", 1:p)])
    psi <- pmin(pmax(as.vector(X %*% beta), -20), 20)

    mZ <- pmax(as.vector(xi %*% sp$Msp), 1e-10)
    MZ <- as.vector(xi %*% sp$Isp)
    MT <- as.vector(xi %*% sp$Isp_T0)

    delta <- data_input$Event
    event_ll <- sum(ifelse(delta > 0, log(mZ) + psi, 0))
    risk_ll  <- sum((MZ - MT) * exp(psi))
    nll <- -(event_ll - risk_ll)
    if (!is.finite(nll)) return(1e10)
    nll
  }

  init <- c(rep(1, p), rep(0.3, sp$p))
  fit <- optim(
    par = init, fn = neg_loglik_null,
    data_input = data, sp = sp, p = p,
    method = "L-BFGS-B",
    lower = c(rep(-5, p), rep(0, sp$p)),
    upper = c(rep(5, p), rep(5, sp$p))
  )

  xi_hat <- fit$par[(p + 1):(p + sp$p)]^2
  list(beta = fit$par[1:p], xi = xi_hat, b = xi_hat, logLik = -fit$value)
}


##########################SUP 统计量######################
compute_SUP_stat_score <- function(data, fit0, sp, eta_grid, p) {

  beta_hat <- fit0$beta
  xi_hat   <- fit0$xi

  X  <- data.matrix(data[, paste0("X", 1:p), drop = FALSE])
  delta <- data$Event
  u_cp <- as.numeric(data$U)

  MZ <- as.vector(xi_hat %*% sp$Isp)
  MT <- as.vector(xi_hat %*% sp$Isp_T0)

  Xt <- cbind(1, X)  # X̃ = (1, X1,...,Xp)

  psi <- as.vector(X %*% beta_hat)
  mu  <- as.numeric((MZ - MT) * exp(psi))
  r   <- as.numeric(delta - mu)

  SUP_vals <- numeric(length(eta_grid))

  for (k in seq_along(eta_grid)) {

    eta <- eta_grid[k]
    Ieta <- as.numeric(u_cp > eta)

    U_score <- colSums(sweep(Xt, 1, r * Ieta, "*"))

    Sigma <- matrix(0, ncol(Xt), ncol(Xt))
    idx <- which(Ieta == 1)
    for (i in idx) {
      xi <- matrix(as.numeric(Xt[i, ]), ncol = 1)
      Sigma <- Sigma + (r[i]^2) * (xi %*% t(xi))
    }
    Sigma <- Sigma + diag(1e-6, ncol(Sigma))

    SUP_vals[k] <- as.numeric(t(U_score) %*% solve(Sigma, U_score))
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
  sp <- select_knots(data, p)$sp

  U <- data$U
  u_min <- min(U, na.rm = TRUE)
  u_max <- max(U, na.rm = TRUE)

  # 为避免极端点数值不稳定，可在两端稍作截断
  lower <- u_min + trim * (u_max - u_min)
  upper <- u_max - trim * (u_max - u_min)

  eta_grid <- seq(lower, upper, length.out = k)

  fit0 <- fit_null_model(data, sp, p)

  sup_obs <- compute_SUP_stat_score(
    data, fit0, sp, eta_grid, p
  )$SUP

  sup_perm <- numeric(B_perm)

  # 构造置换分布
  for (b in 1:B_perm) {
    data_perm <- data
    data_perm$U <- sample(data$U) # 打乱变点协变量 U

    # 重新拟合零假设模型
    fit0_perm <- fit_null_model(data_perm, sp, p)

    sup_perm[b] <- compute_SUP_stat_score(
      data_perm, fit0_perm, sp, eta_grid, p
    )$SUP
  }

  crit <- quantile(sup_perm, probs = 1 - config$alpha)

  list(
    SUP_obs = sup_obs, # 一个数
    sup_perm = sup_perm, # 一列数（分布）
    crit = crit # 一个分位数
  )
}


# 多次循环：每轮随机抽取整数作 set.seed，再对同一数据做 SUP 置换检验；保存 seed 与 p 值，可选写 Excel。
# generator_seed 为整数时先 set.seed(generator_seed)，再 sample 各轮检验用种子（整表可复现）；为 NA 时不固定。
replicate_cp_test_random_seeds <- function(
    data,
    config,
    n_rep = 100L,
    outfile = NULL,
    generator_seed = NA,
    max_seed = 2147483647L
) {
  n_rep <- as.integer(n_rep)
  if (n_rep < 1L) {
    stop("n_rep 须为正整数")
  }
  max_seed <- as.integer(max_seed)
  if (!is.na(generator_seed)) {
    set.seed(as.integer(generator_seed))
  }
  seeds <- sample.int(max_seed, size = n_rep, replace = FALSE)

  df <- data.frame(
    replication = seq_len(n_rep),
    seed = NA_integer_,
    p_value = NA_real_,
    SUP_obs = NA_real_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(n_rep)) {
    s <- seeds[i]
    set.seed(s)
    res <- single_cp_test_score(data, config)
    df$seed[i] <- s
    df$p_value[i] <- mean(res$sup_perm >= res$SUP_obs)
    df$SUP_obs[i] <- res$SUP_obs
  }

  if (!is.null(outfile)) {
    if (!requireNamespace("openxlsx", quietly = TRUE)) {
      stop("写入 Excel 需要安装 openxlsx 包")
    }
    openxlsx::write.xlsx(df, outfile, rowNames = FALSE)
    cat("replicate_cp_test_random_seeds: 已保存 ", outfile, "\n", sep = "")
  }

  invisible(list(df = df, outfile = outfile))
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
      u.mean = config$u.mean,
      u.sd   = config$u.sd,
      x1.mean = config$x1.mean,
      x1.sd   = config$x1.sd,
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

  ### ===== 跨平台并行 =====
  if (.Platform$OS.type == "windows") {

    cl <- parallel::makeCluster(ncore)

    project_dir <- config$project_dir

    parallel::clusterExport(
      cl,
      varlist = c("config", "project_dir"),
      envir = environment()
    )

    parallel::clusterEvalQ(cl, {
      source(file.path(project_dir, "data/TimeindepLTRC_gnrt_ChangepointPH.R"))
      source(file.path(project_dir, "estimation/estimate.R"))
      source(file.path(project_dir, "main/MC.R"))
      source(file.path(project_dir, "main/config.R"))
      source(file.path(project_dir, "main/save.R"))
      setwd(project_dir)
      source(file.path(project_dir, "estimation/plot_baseline.R"))
      source(file.path(project_dir, "test/hypothesis_test.R"))
    })

    res_list <- parallel::parLapply(cl, 1:B_mc, one_run)

    parallel::stopCluster(cl)

  } else {

    res_list <- parallel::mclapply(1:B_mc, one_run, mc.cores = ncore)

  }

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