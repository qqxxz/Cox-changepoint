source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/data/TimeindepLTRC_gnrt_ChangepointPH.R")
source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/estimation/estimate.R")
source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/code/MC.R")
source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/code/config.R")
source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/code/save.R")
setwd("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main")
source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/estimation/plot_baseline.R")

set.seed(123)

data_pilot <- TimeindepLTRC_gnrt_ChangepointPH(
  N = 5000,
  Distribution = SIM_CONFIG$Distribution,
  eta = SIM_CONFIG$eta,
  Beta = SIM_CONFIG$Beta,
  Gamma = SIM_CONFIG$Gamma,
  truncation.percent = SIM_CONFIG$truncation,
  censor.percent = SIM_CONFIG$censor,
  x1.mean = SIM_CONFIG$x1.mean,
  x1.sd   = SIM_CONFIG$x1.sd,
  x2.mean = SIM_CONFIG$x2.mean,
  x2.sd   = SIM_CONFIG$x2.sd,
  adjust.censor = TRUE
)$Data

knots_fixed <- select_knots(data_pilot, p = SIM_CONFIG$p)

cat("Fixed knots:", paste(knots_fixed, collapse = ", "), "\n")


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

  # 更新 SIM_CONFIG
  SIM_CONFIG$n <- exp_grid$n[i]
  SIM_CONFIG$truncation <- exp_grid$truncation[i]
  SIM_CONFIG$censor <- exp_grid$censor[i]

  cat("Running experiment", i, ": n =", SIM_CONFIG$n,
      ", truncation =", SIM_CONFIG$truncation,
      ", censor =", SIM_CONFIG$censor, "\n")

  # 运行仿真
  res <- run_simulation(SIM_CONFIG, knots_fixed)

  # 计算感兴趣参数索引
  p <- SIM_CONFIG$p
  K <- ncol(res$par_mat) - (2*p + 1)
  idx_interest <- c(1:p, (p+1):(2*p), (2*p + K + 1))

  par_mat_interest <- res$par_mat[, idx_interest]
  se_mat_interest  <- res$se_mat[, idx_interest]
  true_par <- c(SIM_CONFIG$Beta, SIM_CONFIG$Gamma, SIM_CONFIG$eta)

  # 统计总结
  summary_res <- summary_MC(par_mat_interest, se_mat_interest, true_par)

  # 保存 Excel
  save_MC_to_excel(
    res        = res,
    config     = SIM_CONFIG,
    true_par   = true_par,
    out_dir    = "/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/results",
    summary_df = summary_res
  )

  # 保存到列表
  results_all[[i]] <- list(
    config = SIM_CONFIG,
    result = res,
    summary = summary_res
  )
}
