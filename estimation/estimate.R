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

################### 负对数似然 #######################
neg_loglik_profile <- function(par, data_input, sp, p, eta) {
  beta     <- par[1:p]
  gamma    <- par[(p + 1):(2 * p + 1)]
  theta_xi <- par[(2 * p + 2):(2 * p + 1 + sp$p)]
  xi       <- theta_xi^2

  X  <- as.matrix(data_input[, paste0("X", 1:p)]) # 设计矩阵 X = (X1,...,Xp)
  Xt <- cbind(1, X) # X̃ = (1, X1,...,Xp)
  if (!"U" %in% names(data_input)) stop("data 须包含变点协变量 U")

  Z     <- data_input$Stop # 观察到的生存时间
  T0    <- data_input$Start # 左截断时间
  delta <- data_input$Event # 事件指示变量
  n     <- nrow(data_input)
  U     <- data_input$U
  ind   <- as.numeric(U > eta)  

  #psi <- as.vector(X %*% beta + (Xt %*% gamma) * ind)
  psi <- pmin(pmax(as.vector(X %*% beta + (Xt %*% gamma) * ind), -30), 30)
  mZ <- pmax(as.vector(xi %*% sp$Msp), 1e-10)
  MZ <- as.vector(xi %*% sp$Isp)
  MT <- as.vector(xi %*% sp$Isp_T0)

  if (any(!is.finite(psi)) || any(!is.finite(mZ)) ||
      any(!is.finite(MZ)) || any(!is.finite(MT))) {
    cat("数值问题: psi, mZ, MZ, MT 中存在非有限值\n")
    return(1e10)
  }

  nll <- -(sum(ifelse(delta > 0, log(mZ) + psi, 0)) - sum((MZ - MT) * exp(psi)))
  if (!is.finite(nll)){
    cat("数值问题: nll 非有限值\n")
    return(1e10)
  } 
  nll
}

###################### 给定 η 拟合 #####################
fit_given_eta <- function(data, sp, p, eta) {
  sp <- spline_setup(data$Stop, sp$k, sp$od, knotc = sp$knotc,
                     times_T0 = data$Start)
  dtime <- pmax(data$Stop - data$Start, 1e-6)
  h0 <- sqrt(sum(data$Event) / sum(dtime))

  #init <- c(rep(0.6, 2 * p + 1), rep(max(h0, 0.1), sp$p))
  beta_init <- rep(0, p)  # 假设 β 的初值为 0
  gamma_init <- rep(0, p + 1)  # 假设 γ 的初值为 0
  xi_init <- rep(1, sp$p)  # 假设 ξ 的初值为 1
  init <- c(beta_init, gamma_init, xi_init)

  fit <- tryCatch(
    optim(
      par = init,
      fn  = neg_loglik_profile,
      data_input = data,
      sp = sp,
      p = p,
      eta = eta,
      method = "L-BFGS-B",
      lower = c(rep(-5, 2 * p + 1), rep(0, sp$p)),
      upper = c(rep(5, 2 * p + 1), rep(5, sp$p)),
      control = list(maxit = 5000, factr = 1e5, pgtol = 1e-10),
      hessian = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(fit) || !is.finite(fit$value)) return(NULL)

  list(par = fit$par, value = fit$value, convergence = fit$convergence)
}

######################网格搜索#####################
profile_eta <- function(data, sp, p,
                        probs = seq(0.15, 0.85, by = 0.05)) {
 # 构造搜索范围，在 U 的 15%～85% 分位数之间，每 5% 取一个候选 η
  eta_grid <- unique(
    stats::quantile(data$U,
                    probs = probs,
                    na.rm = TRUE)
  )

  best_val <- Inf
  best_fit <- NULL
  best_eta <- NA

  for (eta in eta_grid) {
    fit <- fit_given_eta(data, sp, p, eta)
    if (is.null(fit)) next
    if (fit$value < best_val) {
      best_val <- fit$value
      best_fit <- fit
      best_eta <- eta
    }
  }

  if (is.null(best_fit)) stop("所有 eta 的优化均失败")

  list(eta = best_eta, par = best_fit$par, value = best_val,
       convergence = best_fit$convergence)
}

##################### BIC 选结点 + 拟合 #######################
# 返回完整拟合结果（含 sp）；hessian=TRUE 时额外计算 Hessian
select_knots <- function(data, p,
                         k_candidates = c(3L, 5L, 7L, 9L),
                         od = 3L,
                         hessian = TRUE) {
  Y <- data$Stop[!is.na(data$Stop)]
  n <- nrow(data)
  best_BIC  <- Inf
  best_sp   <- NULL
  best_prof <- NULL

  for (k in k_candidates) {
    sp <- spline_setup(data$Stop, k, od, Y_for_knots = Y, times_T0 = data$Start)

    prof <- tryCatch(profile_eta(data, sp, p), error = function(e) NULL)
    if (is.null(prof)) next

    # BIC = -logL + log(n) *df；df 为 β,γ,ξ,η（同示例思路，ξ 维数随 k 变化）
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
  }
  if (is.null(best_prof)) {
    best_prof <- profile_eta(data, best_sp, p)
  }

  prof <- best_prof
  sp <- best_sp

  if (hessian) {
    hess_fit <- tryCatch(
      optim(
        par = prof$par,
        fn  = neg_loglik_profile,
        data_input = data,
        sp = sp,
        p = p,
        eta = prof$eta,
        method = "L-BFGS-B",
        lower = c(rep(-5, 2 * p + 1), rep(0, sp$p)),
        upper = c(rep(5, 2 * p + 1), rep(5, sp$p)),
        control = list(maxit = 500, factr = 1e7, pgtol = 1e-8),
        hessian = TRUE
      ),
      error = function(e) NULL
    )
    if (!is.null(hess_fit) && is.finite(hess_fit$value)) {
      prof$hessian <- hess_fit$hessian
    }
  }

  par_inner <- prof$par
  beta     <- par_inner[1:p]
  gamma    <- par_inner[(p + 1):(2 * p + 1)]
  theta_xi <- par_inner[(2 * p + 2):(2 * p + 1 + sp$p)]
  xi       <- theta_xi^2
  logLik_val <- -prof$value
  df <- 2 * p + 1 + sp$p + 1L
  BIC_val <- prof$value + log(n) * df

  message("最佳 BIC: ", round(BIC_val, 3))
  message("最佳内结点数 k = ", sp$k,
          "，样条基个数 = ", sp$p, " (= k + od = ", sp$k, " + ", sp$od, ")")

  list(
    par = c(par_inner, prof$eta),
    beta = beta,
    gamma = gamma,
    xi = xi,
    b = xi,
    eta = prof$eta,
    sp = sp,
    knots = sp$knotc,
    k = sp$k,
    K = sp$p,
    logLik = logLik_val,
    BIC = BIC_val,
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
