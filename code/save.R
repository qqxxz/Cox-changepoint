library(openxlsx)

make_xlsx_name <- function(config, prefix = "MC") {
  paste0(
    prefix,
    "_Dist=", config$Distribution,
    "_n=", config$n,
    "_LT=", config$truncation,
    "_C=",  config$censor,
    "_B=",  config$B,
    ".xlsx"
  )
}

save_MC_to_excel <- function(res, config, true_par, out_dir, summary_df = NULL) {
    
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  wb <- createWorkbook()

  ## ===== summary sheet（论文主表）=====
  if (is.null(summary_df)) {
    p <- SIM_CONFIG$p
    K <- ncol(res$par_mat) - (2*p + 1)
    idx_interest <- c(
      1:p,  # beta
      (p+1):(2*p),  # gamma
      (2*p + K + 1)  # eta
    )
    par_mat_interest <- res$par_mat[, idx_interest, drop = FALSE]
    se_mat_interest  <- res$se_mat[, idx_interest, drop = FALSE]
    summary_df <- summary_MC(par_mat_interest, se_mat_interest, true_par)
  }

  addWorksheet(wb, "summary")
  writeData(wb, "summary", summary_df, rowNames = TRUE)

  ## ===== config sheet =====
  cfg_df <- data.frame(
    Parameter = names(config),
    Value = sapply(config, function(x) {
      if (is.function(x)) "function"
      else paste(x, collapse = ",")
    })
  )
  addWorksheet(wb, "config")
  writeData(wb, "config", cfg_df)

  ## ===== par_mat sheet =====
  addWorksheet(wb, "par_mat")
  writeData(wb, "par_mat", res$par_mat)

  ## ===== AIE sheet =====
  AIE_df <- data.frame(
    AIE_SF  = res$AIE_SF,
    AIE_CHF = res$AIE_CHF
  )
  addWorksheet(wb, "AIE")
  writeData(wb, "AIE", AIE_df)
  
  writeData(
  wb, "AIE",
  data.frame(
    AIE_SF_MC  = res$AIE_SF_MC,
    AIE_CHF_MC = res$AIE_CHF_MC
  ),
  startRow = nrow(AIE_df) + 3,
  colNames = TRUE
  ) 

  ## ===== info sheet（可选）=====
  addWorksheet(wb, "info")
  writeData(
    wb, "info",
    data.frame(
      RunTime = Sys.time(),
      R_version = R.version.string
    )
  )

  ## ===== 保存 =====
  fname <- make_xlsx_name(config)
  saveWorkbook(
    wb,
    file = file.path(out_dir, fname),
    overwrite = TRUE
  )

  message("Saved: ", fname)
}
