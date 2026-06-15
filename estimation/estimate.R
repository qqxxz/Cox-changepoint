################# M/I 样条基线：h₀(t)=ξ'M(t)，Λ₀(t)=ξ'I(t) #################
library(orthogonalsplinebasis)
library(stats)

# 构造 M/I 样条基矩阵（在 times 上），并在观测时刻上计算 M/I 矩阵
# sp$p = k + od（样条基个数），与回归维数 p 不同
# knotc：已选定的结点；为 NULL 时用 Y_for_knots 的分位数新建结点
spline_setup <- function(times, k, od = 3L, knotc = NULL,
                         Y_for_knots = NULL, times_T0 = NULL) {
  if (is.null(knotc)) {
    Y <- if (!is.null(Y_for_knots)) Y_for_knots else times
    Y <- Y[is.finite(Y) & Y >= 0]
    if (length(Y) < 2L) stop("构造结点的时间点不足")
    knotc <- expand.knots(
      quantile(Y, probs = seq(0, 1, 1 / (k + 1)), na.rm = TRUE),
      order = od
    )
  }
  n_basis <- k + od
  basc  <- SplineBasis(knotc, order = od) # B1​(t),…,Bp​(t) p=k+od
  dbasc <- deriv(basc) # B1′​(t),…,Bp′​(t)

  valc  <- evaluate(basc, times) # dim(valc)=n*p
  valdc <- evaluate(dbasc, times) # Bi′​(Yj​)
  if (is.null(dim(valc))) valc <- matrix(valc, nrow = 1L)
  valc[!is.finite(valc)]  <- 0
  valdc[!is.finite(valdc)] <- 0
  n <- length(times)
  Isp <- matrix(0, n_basis, n)
  Msp <- matrix(0, n_basis, n)
  for (i in seq_len(n_basis)) {
    for (j in seq_len(n)) {
      Isp[i, j] <- sum(valc[j, i:n_basis]) # Ii​(Yj​)=∑^p_{r=i} ​Br​(Yj​)
      Msp[i, j] <- sum(valdc[j, i:n_basis]) # Msp(ij)​=∑^p_{r=i} ​Mr′​(Yj​)
    }
  }
  cat("样条基矩阵维度: Isp =", dim(Isp), ", Msp =", dim(Msp), "\n")
  out <- list(k = k, od = od, p = n_basis, knotc = knotc, basc = basc, dbasc = dbasc,
              Isp = Isp, Msp = Msp)

  if (!is.null(times_T0)) {
    valc_T0 <- evaluate(basc, times_T0)
    if (is.null(dim(valc_T0))) valc_T0 <- matrix(valc_T0, nrow = 1L)
    valc_T0[!is.finite(valc_T0)] <- 0
    n0 <- length(times_T0)
    Isp_T0 <- matrix(0, n_basis, n0)
    for (i in seq_len(n_basis)) {
      for (j in seq_len(n0)) {
        Isp_T0[i, j] <- sum(valc_T0[j, i:n_basis])
      }
    }
    out$Isp_T0 <- Isp_T0
  }
  out
}

# 基线风险 / 累积风险（用于作图与 AIE）
m0 <- function(t, xi, sp) {
  valdc <- evaluate(sp$dbasc, t)
  if (is.null(dim(valdc))) valdc <- matrix(valdc, nrow = 1L)
  p <- sp$p
  out <- vapply(seq_len(nrow(valdc)), function(j) {
    acd <- valdc[j, ]
    sum(xi * vapply(seq_len(p), function(i) sum(acd[i:p]), numeric(1)))
  }, numeric(1))
  pmax(out, 1e-10)
}

M0 <- function(t, xi, sp) {
  valc <- evaluate(sp$basc, t)
  if (is.null(dim(valc))) valc <- matrix(valc, nrow = 1L)
  p <- sp$p
  out <- vapply(seq_len(nrow(valc)), function(j) {
    ac <- valc[j, ]
    sum(xi * vapply(seq_len(p), function(i) sum(ac[i:p]), numeric(1)))
  }, numeric(1))
  pmax(out, 1e-10)
}

################### 似然、梯度 #######################
# 计算偏似然核心量；block 指定返回哪部分负梯度，供分块 optim 用
lik_grad <- function(par, data, sp, p, eta, block = c("fn", "all", "beta", "gamma", "theta")) {
  block <- match.arg(block)
  if (!"U" %in% names(data)) stop("data 须包含变点协变量 U")

  beta  <- par[1:p]
  gamma <- par[(p + 1):(2 * p + 1)]
  theta <- par[(2 * p + 2):(2 * p + 1 + sp$p)]

  X  <- as.matrix(data[, paste0("X", 1:p), drop = FALSE]) # 设计矩阵 X = (X1,...,Xp)
  Xt <- cbind(1, X) # X̃ = (1, X1,...,Xp)
  delta <- data$Event
  ind   <- as.numeric(data$U > eta)
  #psi <- as.vector(X %*% beta + (Xt %*% gamma) * ind)
  psi   <- pmin(pmax(as.vector(X %*% beta + (Xt %*% gamma) * ind), -30), 30)
  exp_psi <- exp(psi)

  xi <- theta^2
  mZ <- pmax(as.vector(xi %*% sp$Msp), 1e-10)
  MZ <- as.vector(xi %*% sp$Isp)
  MT <- as.vector(xi %*% sp$Isp_T0)
  risk <- (MZ - MT) * exp_psi
  r    <- delta - risk

  if (any(!is.finite(c(psi, mZ, MZ, MT)))) return(1e10)
  if (block == "fn") {
    nll <- -(sum(ifelse(delta > 0, log(mZ) + psi, 0)) - sum(risk))
    if (!is.finite(nll)) return(1e10)
    return(nll)
  }

  g_beta  <- -colSums(X * r)
  g_gamma <- -colSums((Xt * ind) * r)
  g_theta <- -2 * theta * colSums(
    t(sp$Msp) * (delta / mZ) - t(sp$Isp - sp$Isp_T0) * exp_psi
  )

  switch(block,
         beta  = g_beta,
         gamma = g_gamma,
         theta = g_theta,
         all   = c(g_beta, g_gamma, g_theta))
}

neg_loglik_profile <- function(par, data_input, sp, p, eta) {
  lik_grad(par, data_input, sp, p, eta, "fn")
}
######################### 三明治渐进方差 #########################
asymptotic_se <- function(par, data, sp, p, eta) {
  n <- nrow(data)
  d_reg <- 2 * p + 1
  d <- d_reg + sp$p

  beta  <- par[1:p]
  gamma <- par[(p + 1):(2 * p + 1)]
  theta <- par[(2 * p + 2):(2 * p + 1 + sp$p)]
  X  <- as.matrix(data[, paste0("X", 1:p), drop = FALSE])
  Xt <- cbind(1, X)
  delta <- data$Event
  ind   <- as.numeric(data$U > eta)
  psi   <- pmin(pmax(as.vector(X %*% beta + (Xt %*% gamma) * ind), -30), 30)
  exp_psi <- exp(psi)
  xi <- theta^2
  mZ <- pmax(as.vector(xi %*% sp$Msp), 1e-10)
  MZ <- as.vector(xi %*% sp$Isp)
  MT <- as.vector(xi %*% sp$Isp_T0)
  r  <- delta - (MZ - MT) * exp_psi

  S <- matrix(0, n, d)
  S[, 1:p] <- X * r
  S[, (p + 1):d_reg] <- (Xt * ind) * r
  for (k in seq_len(sp$p)) {
    S[, d_reg + k] <-
      (delta / mZ) * 2 * theta[k] * sp$Msp[k, ] -
      exp_psi * 2 * theta[k] * (sp$Isp[k, ] - sp$Isp_T0[k, ])
  }
  B <- crossprod(S) / n

  g0 <- lik_grad(par, data, sp, p, eta, "all")  # 总梯度
  eps <- 1e-5
  H <- matrix(0, d, d)
  for (j in seq_len(d)) {
    par_p <- par
    par_p[j] <- par_p[j] + eps
    H[, j] <- (lik_grad(par_p, data, sp, p, eta, "all") - g0) / eps # 有限差分
  }
  H <- (H + t(H)) / 2 # 对称化
  Hinv <- tryCatch(solve(H), error = function(e) {
    if (requireNamespace("MASS", quietly = TRUE)) MASS::ginv(H) else NULL
  }) # 求逆

  se_reg <- rep(NA_real_, d_reg)
  vcov <- NULL
  if (!is.null(Hinv)) {
    vcov <- Hinv %*% B %*% Hinv * n # 三明治公式
    se_reg <- sqrt(pmax(diag(vcov)[seq_len(d_reg)], 0))
  }

  # η：profile 负对数似然（每点重新优化 β,γ,ξ）的二阶导 → 三明治型 profile 方差
  prof_nll <- function(e) {
    fit <- fit_given_eta(data, sp, p, e, max_iter = 12L, par0 = par)
    if (is.null(fit)) 1e10 else fit$value
  }
  u_rng <- diff(range(data$U, na.rm = TRUE))
  h <- max(0.02 * u_rng, 0.05)
  vals <- vapply(c(eta - h, eta, eta + h), prof_nll, numeric(1))
  d2 <- (vals[3] - 2 * vals[2] + vals[1]) / h^2
  se_eta <- sqrt(1 / max(d2, 1e-8))
  # 网格搜索带来的额外不确定性（profile_eta 在 U 分位数上网格化）
  eta_grid <- unique(stats::quantile(data$U, probs = seq(0.15, 0.85, by = 0.05),
                                   na.rm = TRUE))
  if (length(eta_grid) >= 2L) {
    grid_step <- min(diff(sort(eta_grid)))
    se_eta <- sqrt(se_eta^2 + (grid_step / 2)^2)
  }

  list(
    se_beta  = se_reg[1:p],
    se_gamma = se_reg[(p + 1):d_reg],
    se_eta   = se_eta,
    se_par   = c(se_reg, se_eta),
    vcov     = vcov
  )
}

###################### 给定 η：分块 optim + 联合 polish #####################
fit_given_eta <- function(data, sp, p, eta, max_iter = 30L, tol = 1e-5,
                          par0 = NULL, polish = TRUE) {
  sp <- spline_setup(data$Stop, sp$k, sp$od, knotc = sp$knotc,
                     times_T0 = data$Start)
  if (!is.null(par0)) {
    beta  <- par0[1:p]
    gamma <- par0[(p + 1):(2 * p + 1)]
    theta <- par0[(2 * p + 2):(2 * p + 1 + sp$p)]
  } else {
    h0 <- sqrt(sum(data$Event) / sum(pmax(data$Stop - data$Start, 1e-6)))
    beta  <- rep(0, p)
    gamma <- rep(0, p + 1)
    theta <- rep(max(h0, 0.1), sp$p)
  }
  opt_ctrl <- list(maxit = 500, factr = 1e7, pgtol = 1e-8)
  lo <- c(rep(-5, p), rep(-5, p + 1), rep(0, sp$p))
  hi <- c(rep(5, p), rep(5, p + 1), rep(5, sp$p))

  for (iter in seq_len(max_iter)) {
    old <- c(beta, gamma, theta)

    fit_b <- tryCatch(optim(
      beta,
      fn = function(b) lik_grad(c(b, gamma, theta), data, sp, p, eta, "fn"),
      gr = function(b) lik_grad(c(b, gamma, theta), data, sp, p, eta, "beta"),
      method = "L-BFGS-B", lower = lo[1:p], upper = hi[1:p], control = opt_ctrl
    ), error = function(e) NULL)
    if (!is.null(fit_b)) beta <- fit_b$par

    fit_g <- tryCatch(optim(
      gamma,
      fn = function(g) lik_grad(c(beta, g, theta), data, sp, p, eta, "fn"),
      gr = function(g) lik_grad(c(beta, g, theta), data, sp, p, eta, "gamma"),
      method = "L-BFGS-B", lower = lo[(p + 1):(2 * p + 1)], upper = hi[(p + 1):(2 * p + 1)],
      control = opt_ctrl
    ), error = function(e) NULL)
    if (!is.null(fit_g)) gamma <- fit_g$par

    fit_t <- tryCatch(optim(
      theta,
      fn = function(t) lik_grad(c(beta, gamma, t), data, sp, p, eta, "fn"),
      gr = function(t) lik_grad(c(beta, gamma, t), data, sp, p, eta, "theta"),
      method = "L-BFGS-B", lower = lo[(2 * p + 2):length(lo)], upper = hi[(2 * p + 2):length(hi)],
      control = opt_ctrl
    ), error = function(e) NULL)
    if (!is.null(fit_t)) theta <- fit_t$par

    if (max(abs(c(beta, gamma, theta) - old)) < tol) break
  }

  par <- c(beta, gamma, theta)
  if (polish) {
    fit_j <- tryCatch(optim(
      par,
      fn = function(x) lik_grad(x, data, sp, p, eta, "fn"),
      gr = function(x) lik_grad(x, data, sp, p, eta, "all"),
      method = "L-BFGS-B", lower = lo, upper = hi, control = opt_ctrl
    ), error = function(e) NULL)
    if (!is.null(fit_j) && is.finite(fit_j$value)) par <- fit_j$par
  }

  value <- lik_grad(par, data, sp, p, eta, "fn")
  if (!is.finite(value) || value >= 1e9) return(NULL)
  list(par = par, value = value, convergence = 0L, sp = sp)
}

###################### η 网格 + 局部 refine #####################
profile_eta <- function(data, sp, p,
                        probs = seq(0.15, 0.85, by = 0.05)) {
  eta_grid <- unique(stats::quantile(data$U, probs = probs, na.rm = TRUE))
  best_val <- Inf
  best_fit <- NULL
  best_eta <- NA
  par_cur  <- NULL

  for (eta in eta_grid) {
    fit <- fit_given_eta(data, sp, p, eta, max_iter = 15L,
                         par0 = par_cur, polish = FALSE)
    if (is.null(fit)) next
    par_cur <- fit$par
    if (fit$value < best_val) {
      best_val <- fit$value
      best_fit <- fit
      best_eta <- eta
    }
  }
  if (is.null(best_fit)) stop("所有 eta 的优化均失败")

  h <- max(0.02 * diff(range(data$U, na.rm = TRUE)), 0.05)
  fine_grid <- unique(c(seq(best_eta - 2 * h, best_eta + 2 * h, length.out = 11),
                        best_eta))
  fine_grid <- fine_grid[fine_grid >= min(data$U) & fine_grid <= max(data$U)]

  for (eta in fine_grid) {
    fit <- fit_given_eta(data, sp, p, eta, max_iter = 15L,
                         par0 = best_fit$par, polish = FALSE)
    if (is.null(fit)) next
    if (fit$value < best_val) {
      best_val <- fit$value
      best_fit <- fit
      best_eta <- eta
    }
  }

  best_fit <- fit_given_eta(data, sp, p, best_eta, max_iter = 30L,
                            par0 = best_fit$par, polish = TRUE)
  if (is.null(best_fit)) stop("最终 eta 拟合失败")

  list(eta = best_eta, par = best_fit$par, value = best_fit$value,
       convergence = best_fit$convergence)
}

##################### BIC 选结点 + 拟合 #######################
select_knots <- function(data, p,
                         k_candidates = c(3L, 5L, 7L, 9L),
                         od = 3L,
                         k_fixed = NULL) {
  Y <- c(data$Start, data$Stop)
  Y <- Y[is.finite(Y) & !is.na(Y)]
  n <- nrow(data)
  best_sp   <- NULL
  best_prof <- NULL

  if (!is.null(k_fixed)) {
    best_sp <- spline_setup(data$Stop, k_fixed, od, Y_for_knots = Y,
                            times_T0 = data$Start)
    best_prof <- profile_eta(data, best_sp, p)
  } else {
    best_BIC <- Inf
    for (k in k_candidates) {
      sp <- spline_setup(data$Stop, k, od, Y_for_knots = Y, times_T0 = data$Start)
      prof <- tryCatch(profile_eta(data, sp, p), error = function(e) NULL)
      if (is.null(prof)) next
      df <- 2 * p + 1 + sp$p + 1L
      BIC <- prof$value + log(n) * df
      if (BIC < best_BIC) {
        best_BIC  <- BIC
        best_sp   <- sp
        best_prof <- prof
      }
    }
    if (is.null(best_sp)) {
      warning("BIC 选择失败，使用默认 k = 5")
      best_sp <- spline_setup(data$Stop, 5L, od, Y_for_knots = Y,
                              times_T0 = data$Start)
      best_prof <- profile_eta(data, best_sp, p)
    }
  }

  prof <- best_prof
  sp <- spline_setup(data$Stop, best_sp$k, best_sp$od, knotc = best_sp$knotc,
                     Y_for_knots = Y, times_T0 = data$Start)
  infer <- asymptotic_se(prof$par, data, sp, p, prof$eta)

  beta  <- prof$par[1:p]
  gamma <- prof$par[(p + 1):(2 * p + 1)]
  xi    <- prof$par[(2 * p + 2):(2 * p + 1 + sp$p)]^2
  BIC_val <- prof$value + log(n) * (2 * p + 1 + sp$p + 1L)

  message("最佳 BIC: ", round(BIC_val, 3))
  message("最佳内结点数 k = ", sp$k,
          "，样条基个数 = ", sp$p, " (= k + od = ", sp$k, " + ", sp$od, ")")

  list(
    par = c(prof$par, prof$eta),
    beta = beta, gamma = gamma, xi = xi, b = xi,
    eta = prof$eta, sp = sp, knots = sp$knotc,
    k = sp$k, K = sp$p,
    logLik = -prof$value, BIC = BIC_val,
    se_beta = infer$se_beta, se_gamma = infer$se_gamma,
    se_eta = infer$se_eta, se_par = infer$se_par,
    vcov = infer$vcov,
    optim = prof
  )
}

##################### AIE #####################
compute_AIE_SF <- function(fit_result, data, config) {
  xi <- fit_result$xi
  sp <- fit_result$sp
  tgrid <- seq(min(data$Stop), max(data$Stop), length.out = 500)
  S_hat <- exp(-M0(tgrid, xi, sp))
  S_true <- exp(-config$H0_true(tgrid))
  dt <- diff(tgrid)[1]
  range_len <- diff(range(data$Stop))
  sum(abs(S_hat - S_true)) * dt / range_len
}

compute_AIE_CHF <- function(fit_result, data, config) {
  xi <- fit_result$xi
  sp <- fit_result$sp
  tgrid <- seq(min(data$Stop), max(data$Stop), length.out = 200)
  H_hat <- M0(tgrid, xi, sp)
  H_true <- config$H0_true(tgrid)
  dt <- diff(tgrid)[1]
  eps <- 1e-6
  range_len <- diff(range(data$Stop))
  sum(abs(H_hat - H_true) )* dt / range_len # / pmax(H_true, eps))
}
