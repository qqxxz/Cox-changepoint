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
    summary_df <- summary_MC(res$par_mat, true_par, res$se_mat)
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

  ## ===== par_mat sheet（beta, gamma, eta）=====
  addWorksheet(wb, "par_mat")
  writeData(wb, "par_mat", res$par_mat)

  if (!is.null(res$se_mat)) {
    addWorksheet(wb, "se_mat")
    writeData(wb, "se_mat", res$se_mat)
  }

  if (!is.null(res$k_vec)) {
    addWorksheet(wb, "k_selected")
    writeData(wb, "k_selected", data.frame(k = res$k_vec))
  }

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
