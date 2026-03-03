source("E:/BNU/BA4/毕业论文/LTRC-changepoint/data/TimeindepLTRC_gnrt_ChangepointPH.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/estimation/estimate.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/code/MC.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/code/config.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/code/save.R")
setwd("E:/BNU/BA4/毕业论文/LTRC-changepoint")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/estimation/plot_baseline.R")


# 参数组合
n_list <- c(300, 500)
trunc_list <- c(0.1, 0.3)
censor_list <- c(0.2, 0.4)

exp_grid <- expand.grid(
  n = n_list,
  truncation = trunc_list,
  censor = censor_list
)

results_all <- list()  # 保存每次实验结果

for(i in 1:nrow(exp_grid)) {

  cat("====================================\n")
  cat("Experiment", i, "of", nrow(exp_grid), "\n")
  
  SIM_CONFIG$exp_id <- i   # 控制不同实验的随机数

  # 更新 SIM_CONFIG
  SIM_CONFIG$n <- exp_grid$n[i]
  SIM_CONFIG$truncation <- exp_grid$truncation[i]
  SIM_CONFIG$censor <- exp_grid$censor[i]

  cat("Running experiment", i, ": n =", SIM_CONFIG$n,
      ", truncation =", SIM_CONFIG$truncation,
      ", censor =", SIM_CONFIG$censor, "\n")

  # 运行仿真
  res <- run_simulation(SIM_CONFIG)

  # 计算感兴趣参数索引
  p <- SIM_CONFIG$p
  K <- ncol(res$par_mat) - (2*p + 1)
  idx_interest <- c(1:p, (p+1):(2*p), (2*p + K + 1))

  par_mat_interest <- res$par_mat[, idx_interest]
  # se_mat_interest  <- res$se_mat[, idx_interest]
  true_par <- c(SIM_CONFIG$Beta, SIM_CONFIG$Gamma, SIM_CONFIG$eta)

  # 统计总结
  summary_res <- summary_MC(par_mat_interest,true_par)

  # 保存 Excel
  save_MC_to_excel(
    res        = res,
    config     = SIM_CONFIG,
    true_par   = true_par,
    out_dir    = "E:/BNU/BA4/毕业论文/LTRC-changepoint/results",
    summary_df = summary_res
  )

  # 保存到列表
  results_all[[i]] <- list(
    config = SIM_CONFIG,
    result = res,
    summary = summary_res
  )
}
