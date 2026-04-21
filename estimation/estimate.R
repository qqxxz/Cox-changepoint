#################工具函数：计算 a_k、m₀(t)、M₀(t)#############3
compute_ak <- function(b, knots) {
    K <- length(b)
    a <- numeric(K)
    a[1] <-  -b[1] * knots[1]
    for (k in 2:K) {
        # 前面累积项 ∑_{i=1}^{k-1} b_i (t_i - t_{i-1})
        sum_part <- sum(b[1:(k-1)] * diff(knots[1:k]))

        # 减去 b_k * t_{k-1}
        a[k] <- sum_part - b[k] * knots[k-1]
    }
    return(a)
}

m0 <- function(t, b, knots) {
    K <- length(b)
    out <- rep(b[1], length(t))
    for (k in 1:K) {
        idx <- (t >= knots[k]) & (t < knots[k+1])
        out[idx] <- b[k]
    }
    # 最后一个区间右侧
    out[t >= knots[K+1]] <- b[K]
    # 第一个区间左侧
    out[t < knots[1]] <- b[1]
    pmax(out, 1e-10)  # 避免数值问题
}

M0 <- function(t, b, knots) {
    a <- compute_ak(b, knots)
    K <- length(b)
    out <- numeric(length(t))
    out[t < knots[1]] <- a[1] + b[1] * t[t < knots[1]]
    for (k in 1:K) {
        idx <- (t >= knots[k]) & (t < knots[k+1])
        out[idx] <- a[k] + b[k] * t[idx]
    }
    # 右侧外推（最后区间延伸）
    out[t >= knots[K+1]] <- a[K] + b[K] * t[t >= knots[K+1]]
    # 左侧外推
    out[t < knots[1]] <- a[1] + b[1] * t[t < knots[1]]
    out <- pmax(out, 1e-10)
}

###################负对数似然函数#######################
neg_loglik_profile <- function(par, data_input, knots, p, eta) {
  K <- length(knots) - 1  # 基准风险被分成的区间数

  beta  <- par[1:p]
  gamma <- par[(p+1):(2*p)]
  theta_b <- par[(2*p+1):(2*p+K)]
  b <- exp(theta_b)   # 保证 b_k > 0

  X <- as.matrix(data_input[, paste0("X", 1:p)])  # 设计矩阵
  if (!is.numeric(X)) stop("X 包含非数值数据")
  if (!is.numeric(beta)) stop("beta 不是数值向量")
  if (p > 1) {
      Xt <- cbind(
        as.matrix(data_input[, paste0("X", 2:p), drop = FALSE]), # 剔除 X1 的协变量
        1  # 添加常数项
  )
  } else {
    Xt <- matrix(1, nrow = nrow(data_input), ncol = 1)
  }

  Z  <- data_input$Stop # 观察到的生存时间
  T0 <- data_input$Start # 左删失时间
  delta <- data_input$Event # 事件指示变量

  n <- nrow(data_input)
  h <- 1.5 * sd(data_input$X1) * n^(-1/4)
  ind <- plogis((data_input$X1 - eta) / h) # 平滑 eta
  # ind <- as.numeric(data_input$X1 > eta) # 变点指示变量

  psi <- as.vector(X %*% beta + (Xt %*% gamma) * ind) 
  psi <- pmin(pmax(psi, -20), 20)  # 限制范围
  mZ <- m0(Z, b, knots)
  if (any(!is.finite(psi)) || any(!is.finite(mZ))) {
    return(1e10)
  }

  MZ <- M0(Z, b, knots)
  MT <- M0(T0, b, knots)

  ll <- sum(delta * (log(mZ) + psi)) -
        sum((MZ - MT) * exp(psi))

  return(-ll)
}

###################### 给定 η 的 BFGS 优化#####################
fit_given_eta <- function(data, knots, p, eta) {

  K <- length(knots) - 1

  init <- c(
    rep(1, 2*p),        # beta, gamma 0 0.8➡️1
    rep(log(mean(data$Event) / mean(data$Stop)), K)    # log(b) 0.1
  )

  fit <- optim(
    par = init,
    fn  = neg_loglik_profile,
    data_input = data,
    knots = knots,
    p = p,
    eta = eta,
    method = "L-BFGS-B",
    lower = c(rep(-5, 2*p), rep(log(1e-2), K)),
    upper = c(rep(5, 2*p), rep(log(10), K)),
    control = list(maxit = 1000, factr = 1e4, pgtol = 1e-10),  # 修改控制参数
    hessian = TRUE
  )

  if (fit$convergence != 0 || !is.finite(fit$value)) {
    return(NULL)
  }

  list(
    par = fit$par,
    value = fit$value,
    hessian = fit$hessian
  )
}

######################网格搜索#####################
get_eta_grid <- function(data, probs = seq(0.1, 0.9, by = 0.02)) {
  unique(quantile(data$X1, probs = probs, na.rm = TRUE))
} # 构造搜索范围，在 X1 的 15%～85% 分位数之间，每 5% 取一个候选 η

profile_eta <- function(data, knots, p) {

  eta_grid <- get_eta_grid(data)

  best_val <- Inf
  best_fit <- NULL
  best_eta <- NA

  for (eta in eta_grid) {

    fit <- fit_given_eta(data, knots, p, eta)

    if (is.null(fit)) next

    if (fit$value < best_val) {
      best_val <- fit$value
      best_fit <- fit
      best_eta <- eta
    }
  } # 对每一个候选值进行模型拟合

  if (is.null(best_fit)) {
    stop("所有 eta 的优化均失败")
  }

  list(
    eta = best_eta,
    par = best_fit$par,
    value = best_val,
    hessian = best_fit$hessian
  )
}

#####################计算 AIC 并选择最佳分割点#######################
select_knots <- function(data, p, quantiles = c(0.05,0.1,0.2,0.3,0.5,0.7,0.9)) {
    Z <- data$Stop[!is.na(data$Stop)]
    if (length(Z) == 0) {
        stop("没有有效的生存时间数据")
    }
    Zq <- quantile(data$Stop, quantiles, na.rm = TRUE) # 计算分位数
    Zq <- unique(Zq)  # 去掉重复或几乎重复的点
    all_sets <- unlist(
        lapply(1:length(Zq), function(k) {
        if (length(Zq) >= k) combn(Zq, k, simplify = FALSE)
        }),
        recursive = FALSE
    ) # 枚举所有子集

    # t0 = 0.9 * min(Z)
    best_AIC <- Inf
    best_knots <- NULL

    for (subset in all_sets) {
        # knots <- c(0.9*min(data$Stop,na.rm = TRUE), subset, max(data$Stop,na.rm = TRUE))
        t_min <- quantile(data$Stop, 0.05, na.rm = TRUE)
        t_max <- quantile(data$Stop, 0.95, na.rm = TRUE)
        knots <- sort(unique(c(t_min, subset, t_max)))  # 去重并排序
        if(length(knots) <= 1) next   # 跳过无效组合
        if(any(diff(knots) < 1e-3)) next  # 跳过间距过小的组合
        # 每段至少 3 个事件
        valid <- TRUE
        for (k in 1:(length(knots)-1)) {
          if (sum(data$Event == 1 &
                  data$Stop >= knots[k] &
                  data$Stop < knots[k+1]) < 2) {
            valid <- FALSE
            break
          }
        }
        if (!valid) next

        K <- length(knots) - 1
        if (K <= 0) next

        init <- c(
            rep(1, 2*p),  # 0.8➡️1
            rep(log(mean(data$Event) / mean(data$Stop)), K)
        ) # 初值

        eta_tmp <- median(data$X1, na.rm = TRUE) # 固定 eta

        fit <- try(
          optim(
            par = init,
            fn  = neg_loglik_profile,
            data_input = data,
            knots = knots,
            p = p,
            eta = eta_tmp,
            method = "L-BFGS-B",
            lower = c(rep(-5, 2*p), rep(log(1e-2), K)),
            upper = c(rep(5, 2*p), rep(log(10), K)),
            control = list(maxit = 500, reltol = 1e-8)
          ),
          silent = TRUE
        )

        if (inherits(fit, "try-error") ||
            fit$convergence != 0 ||
            !is.finite(fit$value)) next

        prof_tmp <- profile_eta(data, knots, p)
        fit_val  <- prof_tmp$value
        AIC_now  <- 2 * (length(prof_tmp$par) + 1) + 2 * fit_val

        if (AIC_now < best_AIC) {
            best_AIC <- AIC_now
            best_knots <- knots
        } # 更新分割点
        
    }

    # 兜底：如果所有优化都失败，则使用最简单的分割
    if (is.null(best_knots)) {
        warning("所有优化失败，使用默认分割点")
        best_knots <- c(0.9*min(data$Stop, na.rm=TRUE),median(data$Stop, na.rm=TRUE), max(data$Stop, na.rm=TRUE))
    }
    print(paste("最佳AIC:", best_AIC))
    print(paste("最佳分割点:", paste(best_knots, collapse = ", ")))
    
    return(best_knots)
}

#####################主函数：拟合分段线性基准风险模型#######################
fit_piecewise <- function(data_input, p, knots) {

  ## 1. 确定分割点
  K <- length(knots) - 1

  if (K <= 0) {
    stop("Invalid knots")
  }

  ## 2. profile likelihood over eta
  prof <- profile_eta(data_input, knots, p)

  best_eta <- prof$eta
  par_inner <- prof$par
  best_val <- prof$value

  ## 3. 参数整理
  beta  <- par_inner[1:p]
  gamma <- par_inner[(p + 1):(2 * p)]

  theta_b <- par_inner[(2 * p + 1):(2 * p + K)]
  b <- exp(theta_b)

  eta <- best_eta

  ## 4. 计算 a_k
  if (K > 0) {
    a <- tryCatch(
      compute_ak(b, knots),
      error = function(e) {
        warning("计算 a_k 失败: ", e$message)
        rep(NA, K)
      }
    )
  } else {
    a <- numeric(0)
  }

  ## 5. logLik 与 AIC
  logLik_val <- -best_val
  AIC_val <- 2 * (length(par_inner) + 1) - 2 * logLik_val
  # +1 是 eta
 
  ## 6. 返回
  return(list(
    par   = c(par_inner, eta),
    beta  = beta,
    gamma = gamma,
    b     = b,
    a     = a,
    eta   = best_eta,
    knots = knots,
    K     = K,
    logLik = logLik_val,
    AIC    = AIC_val,
    optim  = prof
  ))
}


##################### AIE 计算函数 #####################
compute_AIE_SF <- function(fit_result, data, config) {
    # 提取参数
    b <- fit_result$b
    knots <- fit_result$knots
    
    # 样本观测范围
    Z_min <- min(data$Stop, na.rm = TRUE)
    Z_max <- max(data$Stop, na.rm = TRUE)
    
    # 构造时间网格（在观测范围内）
    tgrid <- seq(Z_min, Z_max, length.out = 500)
    
    # 估计的基线生存函数（协变量全为0）
    M0_vals <- M0(tgrid, b, knots)
    S_hat <- exp(-M0_vals)  # 当所有协变量为0时，exp(psi) = exp(0) = 1
    
    # 真实的基线生存函数
    H0_true <- config$H0_true(tgrid)
    S_true <- exp(-H0_true)  # 同样，协变量为0
    
    # 数值积分计算绝对误差
    dt <- diff(tgrid)[1]
    eps <- 1e-6
    integral <- sum(abs(S_hat - S_true) ) * dt

    # 除以观测范围长度
    range_len <- Z_max - Z_min
    if (range_len < 1e-10) return(NA)
    
    integral / range_len
}


compute_AIE_CHF <- function(fit_result, data, config) {
    # 提取参数
    b <- fit_result$b
    knots <- fit_result$knots
    
    # 样本观测范围
    Z_min <- min(data$Stop, na.rm = TRUE)
    Z_max <- max(data$Stop, na.rm = TRUE)
    
    # 构造时间网格
    tgrid <- seq(Z_min, Z_max, length.out = 200)
    
    # 估计的基线累积风险函数
    H_hat <- M0(tgrid, b, knots)  # 协变量为0
    
    # 真实的基线累积风险函数
    H_true <- config$H0_true(tgrid)  # 协变量为0
    
    # 数值积分
    dt <- diff(tgrid)[1]
    eps <- 1e-6
    integral <- sum(abs(H_hat - H_true) / pmax(H_true, eps)) * dt

    # 除以观测范围长度
    range_len <- Z_max - Z_min
    if (range_len < 1e-10) return(NA)
    
    integral / range_len
}