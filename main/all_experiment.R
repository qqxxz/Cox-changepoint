source("E:/BNU/BA4/Cox-changepoint/code/data/TimeindepLTRC_gnrt_ChangepointPH.R")
source("E:/BNU/BA4/Cox-changepoint/code/estimation/estimate.R")
source("E:/BNU/BA4/Cox-changepoint/code/main/MC.R")
source("E:/BNU/BA4/Cox-changepoint/code/main/config.R")
source("E:/BNU/BA4/Cox-changepoint/code/main/save.R")
setwd("E:/BNU/BA4/Cox-changepoint/code")
source("E:/BNU/BA4/Cox-changepoint/code/estimation/plot_baseline.R")

# 参数组合
n_list <- c(300,500)
trunc_list <- c(0.1,0.3)
censor_list <- c(0.2,0.4)

exp_grid <- expand.grid(
  n = n_list,
  truncation = trunc_list,
  censor = censor_list
)

results_all <- list()  # 保存每次实验结果 

for(i in 1:nrow(exp_grid)) {

  cat("====================================\n")
  cat("Experiment", i, "of", nrow(exp_grid), "\n")
  
  # 同一截断/删失设定共用 exp_id，便于 n=300 与 n=500 在相同随机流下比较
  SIM_CONFIG$exp_id <- as.integer(10 * exp_grid$truncation[i] +
                                   100 * exp_grid$censor[i])

  # 更新 SIM_CONFIG
  SIM_CONFIG$n <- exp_grid$n[i]
  SIM_CONFIG$truncation <- exp_grid$truncation[i]
  SIM_CONFIG$censor <- exp_grid$censor[i]

  cat("Running experiment", i, ": n =", SIM_CONFIG$n,
      ", truncation =", SIM_CONFIG$truncation,
      ", censor =", SIM_CONFIG$censor, "\n")

  # 运行仿真
  res <- run_simulation(SIM_CONFIG)

  # par_mat 列顺序: beta, gamma, eta（每次 rep 内 BIC 选结点，维度固定）
  true_par <- c(SIM_CONFIG$Beta, SIM_CONFIG$Gamma, SIM_CONFIG$eta)
  summary_res <- summary_MC(res$par_mat, true_par, res$se_mat)

  # 保存 Excel
  save_MC_to_excel(
    res        = res,
    config     = SIM_CONFIG,
    true_par   = true_par,
    out_dir    = "E:/BNU/BA4/Cox-changepoint/code/results_est/results6",
    summary_df = summary_res
  )

  # 保存到列表
  results_all[[i]] <- list(
    config = SIM_CONFIG,
    result = res,
    summary = summary_res
  )
}
