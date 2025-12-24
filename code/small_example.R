source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/data/TimeindepLTRC_gnrt_ChangepointPH.R")
source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/estimation/estimate.R")
source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/code/MC.R")
source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/code/config.R")
source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/code/save.R")
setwd("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main")
source("/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/estimation/plot_baseline.R")
## ================== 运行实验 ==================
res <- run_simulation(SIM_CONFIG)

p <- SIM_CONFIG$p
K <- ncol(res$par_mat) - (2*p + 1)

idx_interest <- c(
  1:p,  # beta
  (p+1):(2*p),  # gamma
  (2*p + K + 1)  # eta
)

par_mat_interest <- res$par_mat[, idx_interest]
se_mat_interest  <- res$se_mat[, idx_interest]

true_par <- c(
  SIM_CONFIG$Beta,
  SIM_CONFIG$Gamma,
  SIM_CONFIG$eta
)

summary_res <- summary_MC(
  par_mat = par_mat_interest,
  se_mat  = se_mat_interest,
  true_par = true_par
)

## ================== 保存结果 ==================
save_MC_to_excel(
  res      = res,
  config   = SIM_CONFIG,
  true_par = true_par,
  out_dir  = "/Users/xiaodanqi/Desktop/code/【danqi】ltrc-changepoint-main/results",
  summary_df = summary_res  
)
