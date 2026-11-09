library(parallel)
ncore <- detectCores() - 1  # 留一个核给系统


get_true_baseline <- function(config) {

  if (config$Distribution == "Exp") {

    S0 <- function(t) exp(-config$Lambda * t)
    H0 <- function(t) config$Lambda * t

  } else if (config$Distribution %in% c("WD", "WI")) {

    # Weibull (both increasing/decreasing share same H0 form)
    S0 <- function(t) exp(-config$Lambda * t^config$V)
    H0 <- function(t) config$Lambda * t^config$V

  } else if (config$Distribution == "Gtz") {

    Alpha <- if (!is.null(config$Alpha)) config$Alpha else 0.1
    Lambda <- if (!is.null(config$Lambda)) config$Lambda else 0.3
    H0 <- function(t) (Lambda / Alpha) * (exp(Alpha * t) - 1)
    S0 <- function(t) exp(-H0(t))

  } else if (config$Distribution == "Quadratic") {

    a <- config$Coeff$a
    H0 <- function(t) a * t^2
    S0 <- function(t) exp(-H0(t))

  } else if (config$Distribution == "QuadraticLinear") {

    a <- config$Coeff$a
    b <- config$Coeff$b
    H0 <- function(t) a * t^2 + b * t
    S0 <- function(t) exp(-H0(t))

  } else if (config$Distribution == "PiecewiseWeibull") {

    alpha1 <- config$Coeff$alpha1
    lambda1 <- config$Coeff$lambda1
    alpha2 <- config$Coeff$alpha2
    lambda2 <- config$Coeff$lambda2
    t_star <- config$Coeff$t_star

    H0 <- function(t) {
      sapply(t, function(tt) {
        if (tt <= t_star) {
          lambda1 * tt^alpha1
        } else {
          lambda1 * t_star^alpha1 + lambda2 * (tt^alpha2 - t_star^alpha2)
        }
      })
    }
    S0 <- function(t) exp(-H0(t))

  } else {
    stop("Unknown distribution")
  }

  list(S0_true = S0, H0_true = H0)
}

##################### Monte Carlo 仿真 #####################
run_simulation <- function(config) {
  true_base <- get_true_baseline(config)
  config$S0_true <- true_base$S0_true
  config$H0_true <- true_base$H0_true

  B <- config$B

  one_rep <- function(b) {
    mc_seed <- config$base_seed +
              config$seed_multiplier * b +
              987 * config$exp_id

    set.seed(mc_seed)
    cat("Replication:", b, " Seed =", mc_seed, "\n")

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

    fit <- tryCatch(
      select_knots(
        dat, config$p,
        k_candidates = config$k_candidates,
        od = config$spline_order,
        k_fixed = config$k_fixed
      ),
      error = function(e) {
        structure(list(error = conditionMessage(e)), class = "mc_failed")
      }
    )

    if (inherits(fit, "mc_failed")) return(fit)
    if (is.null(fit)) return(structure(list(error = "unknown NULL"), class = "mc_failed"))

    list(
      par = c(fit$beta, fit$gamma, fit$eta),
      se_par = fit$se_par,
      par_full = fit$par,
      k = fit$k,
      xi  = fit$xi,
      AIE_SF  = compute_AIE_SF(fit, dat, config),
      AIE_CHF = compute_AIE_CHF(fit, dat, config),
      sp = fit$sp,
      data  = dat
    )
  }

  # 主进程先跑第 1 次 replication，便于在并行前暴露错误
  cat("Testing replication 1...\n")
  test1 <- one_rep(1L)
  if (inherits(test1, "mc_failed")) {
    stop("Replication 1 failed: ", test1$error)
  }
  cat("Replication 1 OK (k =", test1$k, ")\n")

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
    })

    parallel::clusterExport(cl, varlist = c("one_rep"), envir = environment())

    res_list <- if (B >= 2L) {
      parallel::parLapply(cl, 2:B, one_rep)
    } else {
      list()
    }
    parallel::stopCluster(cl)
    res_list <- c(list(test1), res_list)

  } else {

    res_list <- if (B >= 2L) {
      c(list(test1), parallel::mclapply(2:B, one_rep, mc.cores = ncore))
    } else {
      list(test1)
    }

  }

  failed <- vapply(res_list, inherits, logical(1), "mc_failed")
  if (any(failed)) {
    first_err <- res_list[[which(failed)[1L]]]$error
    stop("MC replication failed: ", first_err)
  }

  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0L) {
    stop("All MC replications failed.")
  }

  par_mat <- do.call(rbind, lapply(res_list, `[[`, "par"))
  se_mat  <- do.call(rbind, lapply(res_list, `[[`, "se_par"))
  k_vec   <- vapply(res_list, `[[`, numeric(1), "k")
  AIE_SF_vec  <- sapply(res_list, `[[`, "AIE_SF")
  AIE_CHF_vec <- sapply(res_list, `[[`, "AIE_CHF")

  rep_plot <- res_list[[1L]]
  plot_baseline_CHF(
    xi_mean = rep_plot$xi,
    sp      = rep_plot$sp,
    data    = rep_plot$data,
    config  = config,
    file = paste0("baseline_CHF_n", config$n,
              "_tr", config$truncation,
              "_c", config$censor, ".png")
  )

  list(
    par_mat = par_mat,
    se_mat = se_mat,
    k_vec = k_vec,
    AIE_SF  = AIE_SF_vec,
    AIE_CHF = AIE_CHF_vec,
    AIE_SF_MC  = mean(AIE_SF_vec,  na.rm = TRUE),
    AIE_CHF_MC = mean(AIE_CHF_vec, na.rm = TRUE)
  )
}

##################### 仿真统计量 #####################
summary_MC <- function(par_mat, true_par, se_mat = NULL) {
  stopifnot(ncol(par_mat) == length(true_par))

  est_mean <- colMeans(par_mat)
  bias <- est_mean - true_par
  sse <- colMeans((par_mat - matrix(true_par,
                                    nrow(par_mat),
                                    length(true_par),
                                    byrow = TRUE))^2)
  see <- apply(par_mat, 2, sd)

  if (is.null(se_mat)) {
    se_mat <- matrix(see, nrow = nrow(par_mat), ncol = ncol(par_mat), byrow = TRUE)
  }
  ase <- colMeans(se_mat, na.rm = TRUE)

  CP95 <- CP99 <- rep(NA, ncol(par_mat))
  for (j in seq_len(ncol(par_mat))) {
    CP95[j] <- mean(
      abs(par_mat[, j] - true_par[j]) <= qnorm(0.975) * se_mat[, j],
      na.rm = TRUE
    )
    CP99[j] <- mean(
      abs(par_mat[, j] - true_par[j]) <= qnorm(0.995) * se_mat[, j],
      na.rm = TRUE
    )
  }

  data.frame(
    Estimate = est_mean,
    Bias = bias,
    SSE = sse,
    SEE = see,
    ASE = ase,
    CP95 = CP95,
    CP99 = CP99
  )
}
