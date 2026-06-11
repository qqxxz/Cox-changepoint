PROJECT_DIR <- "E:/BNU/BA4/Cox-changepoint/code"

source(file.path(PROJECT_DIR, "data/TimeindepLTRC_gnrt_ChangepointPH.R"))
source(file.path(PROJECT_DIR, "estimation/estimate.R"))
source(file.path(PROJECT_DIR, "main/MC.R"))
source(file.path(PROJECT_DIR, "main/config.R"))
source(file.path(PROJECT_DIR, "main/save.R"))
setwd(PROJECT_DIR)
source(file.path(PROJECT_DIR, "estimation/plot_baseline.R"))

## ================== 运行实验 ==================
res <- run_simulation(SIM_CONFIG)

true_par <- c(
  SIM_CONFIG$Beta,
  SIM_CONFIG$Gamma,
  SIM_CONFIG$eta
)

summary_res <- summary_MC(res$par_mat, true_par)

## ================== 保存结果 ==================
save_MC_to_excel(
  res      = res,
  config   = SIM_CONFIG,
  true_par = true_par,
  out_dir  = file.path(PROJECT_DIR, "results"),
  summary_df = summary_res
)
