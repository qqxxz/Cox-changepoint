source("/Users/xiaodanqi/Desktop/LTRC-changepoint/data/TimeindepLTRC_gnrt_ChangepointPH.R")
source("/Users/xiaodanqi/Desktop/LTRC-changepoint/estimation/estimate.R")
source("/Users/xiaodanqi/Desktop/LTRC-changepoint/code/MC.R")
source("/Users/xiaodanqi/Desktop/LTRC-changepoint/code/config.R")
source("/Users/xiaodanqi/Desktop/LTRC-changepoint/code/save.R")
setwd("/Users/xiaodanqi/Desktop/LTRC-changepoint/")
source("/Users/xiaodanqi/Desktop/LTRC-changepoint/estimation/plot_baseline.R")
source("/Users/xiaodanqi/Desktop/LTRC-changepoint/test/hypothesis_test.R")
source("/Users/xiaodanqi/Desktop/LTRC-changepoint/code/config.R")
set.seed(66)
library(openxlsx)

param_grid <- expand.grid(
  n = c(300,500),
  truncation = c(0.1, 0.3),
  censor = c(0.2, 0.4)
)

all_results <- list()
out_dir <- "results_test"
dir.create(out_dir, showWarnings = FALSE)

for(i in 1:nrow(param_grid)) {
  SIM_CONFIG$n <- param_grid$n[i]
  SIM_CONFIG$truncation <- param_grid$truncation[i]
  SIM_CONFIG$censor <- param_grid$censor[i]
  
  cat("=== Running parameter set", i, "of", nrow(param_grid), "===\n")
  cat("n =", SIM_CONFIG$n, "truncation =", SIM_CONFIG$truncation, "censor =", SIM_CONFIG$censor, "\n")
  
  # Type I error (Gamma = 0)
  res_typeI <- run_MC_test(SIM_CONFIG, gamma_vec=c(0,0), out_dir=out_dir)

  # Power (Gamma = 非零)
  res_power <- run_MC_test(SIM_CONFIG, gamma_vec=c(0.9,0.5), out_dir=out_dir)
  
  # 保存汇总
  all_results[[i]] <- list(
    config = SIM_CONFIG,
    TypeI_Error = res_typeI$MC_res$Rejection_Rate,
    Power = res_power$MC_res$Rejection_Rate,
    pval_null = res_typeI$MC_res$p_value,
    pval_alt  = res_power$MC_res$p_value,
    Excel_files = c(res_typeI$file_XLSX, res_power$file_XLSX)
  )
}
