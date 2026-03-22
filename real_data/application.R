source("E:/BNU/BA4/毕业论文/LTRC-changepoint/data/TimeindepLTRC_gnrt_ChangepointPH.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/estimation/estimate.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/code/MC.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/code/config.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/code/save.R")
setwd("E:/BNU/BA4/毕业论文/LTRC-changepoint")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/estimation/plot_baseline.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/test/hypothesis_test.R")


## ----------------1. 读取 Stanford Heart Transplant 数据-----------
library(survival)

data("heart")

str(heart)
summary(heart)
head(heart)

# 描述性统计
heart_clean <- heart
heart_clean$transplant <- as.numeric(as.character(heart$transplant))
vars <- c("start","stop","event","age","year","surgery","transplant")

desc_table <- data.frame(
  Parameter = c("start","stop","event","age","year","surgery","transplant"),
  Mean = sapply(heart_clean[vars], mean, na.rm = TRUE),
  SD   = sapply(heart_clean[vars], sd, na.rm = TRUE),
  N    = sapply(heart_clean[vars], function(x) sum(!is.na(x))),
  Min  = sapply(heart_clean[vars], min, na.rm = TRUE),
  Max  = sapply(heart_clean[vars], max, na.rm = TRUE)
)

desc_table$SE <- desc_table$SD / sqrt(desc_table$N)

# 保留三位小数
desc_table[,2:6] <- round(desc_table[,2:6], 3)

print(desc_table)
library(openxlsx)
write.xlsx(
  desc_table,
  file = "E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/desc_table.xlsx",
  rowNames = FALSE
)
## -------------------2. 数据整理---------------------

dat_raw <- heart

# 数据字段
# start : 左截断时间
# stop  : 观测结束时间
# event : 是否死亡 (1=死亡)
# age   : 年龄
# transplant : 是否接受移植

# 构造 Event
dat_raw$Event <- dat_raw$event

# 构造 Start / Stop
dat_raw$Start <- dat_raw$start
dat_raw$Stop  <- dat_raw$stop

# 变点变量 X1 = 年龄（标准化）
dat_raw$X1 <- as.numeric(scale(dat_raw$age))

# 协变量 X2 = 是否移植
dat_raw$X2 <- as.numeric(as.character(dat_raw$transplant))

# 最终分析数据
data_input <- dat_raw[, c("Start", "Stop", "Event", "X1", "X2")]

# 清理异常值
data_input <- subset(
  data_input,
  is.finite(Start) & is.finite(Stop) &
    is.finite(Event) & is.finite(X1) &
    is.finite(X2) &
    Stop > Start
)

cat("样本量:", nrow(data_input), "\n")
cat("事件数:", sum(data_input$Event), "\n")

summary(data_input)

## ---------------3. 标准 Cox 模型（左截断 + 右删失）作对照------------
library(survival)

cox_fit <- coxph(Surv(Start, Stop, Event) ~ X1 + X2, data = data_input)
summary(cox_fit)
logLik(cox_fit)
AIC(cox_fit)

## -------------------4. 拟合我的模型（两个变量）-----------------
p <- 2

# 选择最优分割点
knots <- select_knots(data_input, p = 2)

# 拟合模型
fit_empirical <- fit_piecewise(data_input, p = p, knots = knots)

# p值
library(MASS)
hessian <- fit_empirical$optim$hessian
hessian <- (hessian + t(hessian)) / 2   # 强制对称
cov_mat <- ginv(hessian)
se <- sqrt(diag(cov_mat))
est <- fit_empirical$par[1:length(se)]
z <- est / se
p_value <- 2 * pnorm(-abs(z))
CI_lower <- est - 1.96 * se
CI_upper <- est + 1.96 * se

cat("========== 模型拟合结果 ==========\n")
cat("eta =", fit_empirical$eta, "\n")
cat("knots =", paste(round(fit_empirical$knots, 4), collapse = ", "), "\n")
cat("beta =", paste(round(fit_empirical$beta, 4), collapse = ", "), "\n")
cat("gamma =", paste(round(fit_empirical$gamma, 4), collapse = ", "), "\n")
cat("b =", paste(round(fit_empirical$b, 4), collapse = ", "), "\n")
cat("a =", paste(round(fit_empirical$a, 4), collapse = ", "), "\n")
cat("logLik =", round(fit_empirical$logLik, 4), "\n")
cat("AIC =", round(fit_empirical$AIC, 4), "\n")

param_names <- c(
  paste0("beta",1:p),
  paste0("gamma",1:p),
  paste0("logb",1:(length(knots)-1))
)

result_table <- data.frame(
  Parameter = param_names,
  Estimate = est,
  SE = se,
  Z = z,
  P_value = p_value,
  CI_lower = CI_lower,
  CI_upper = CI_upper
)

print(result_table)
write.xlsx(
  result_table,
  file = "E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/result_table_2para.xlsx",
  rowNames = FALSE
)

## -------------------4b. 变点存在性检验（两变量：SUP 置换，H0: γ=0）-----------------
# 与 test/hypothesis_test.R 中 single_cp_test_score 一致；k、B_perm、eta.trim 与 simulation 默认可比
cp_test_config_2 <- list( # 设置检验参数
  p = 2L,
  B_perm = 1000L,
  k = 5L,
  eta.trim = 0.01,
  alpha = 0.05
)
set.seed(6)
cp_test_2 <- single_cp_test_score(data_input, cp_test_config_2) # 进行检验
cp_pval_2 <- mean(cp_test_2$sup_perm >= cp_test_2$SUP_obs) # 计算p值
cp_crit_2 <- as.numeric(cp_test_2$crit) # 临界值
cp_reject_2 <- cp_test_2$SUP_obs > cp_crit_2

cat("\n========== 变点存在性检验（两变量模型，SUP 置换）==========\n")
cat("H0: 无变点（γ=0）；变点协变量 X1 在置换下可交换\n")
cat("SUP_obs =", round(cp_test_2$SUP_obs, 6), "\n")
cat("置换 p 值 Pr(T_perm >= SUP_obs) =", round(cp_pval_2, 6), "\n")
cat("临界值 alpha=0.05 =", round(cp_crit_2, 6), "\n")
cat("是否拒绝 H0 (alpha=0.05):", ifelse(cp_reject_2, "是", "否"), "\n")

cp_test_table_2 <- data.frame(
  Model = "2 covariates (age, transplant)",
  SUP_obs = cp_test_2$SUP_obs,
  p_value_perm = cp_pval_2,
  crit_0.05 = cp_crit_2,
  reject_H0 = cp_reject_2,
  B_perm = cp_test_config_2$B_perm,
  k_grid = cp_test_config_2$k,
  stringsAsFactors = FALSE
)
write.xlsx(
  cp_test_table_2,
  file = "E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/cp_test_table_2para.xlsx"
)
## --------------------5. 将 eta 还原到原始年龄尺度------------------
age_mean <- mean(dat_raw$age, na.rm = TRUE)
age_sd   <- sd(dat_raw$age, na.rm = TRUE)

eta_std <- fit_empirical$eta # 变点（标准化尺度）
eta_age <- eta_std * age_sd + age_mean # # 还原到 age 变量尺度

eta_real_age <- eta_age + 48 # 还原真实年龄

cat("变点（标准化尺度）eta =", round(eta_std, 4), "\n")
cat("变点（进入年龄原始尺度）eta_age =", round(eta_age, 4), "\n")
cat("变点（真实年龄）eta_real_age =", round(eta_real_age, 4), "\n")

## --------------6. 作图：估计的基准累积风险函数与基准生存函数--------------
plot_empirical_baseline <- function(fit_result, data_input) {
  b <- fit_result$b
  knots <- fit_result$knots
  
  t_min <- min(data_input$Start, na.rm = TRUE)
  t_max <- max(data_input$Stop, na.rm = TRUE)
  tgrid <- seq(t_min, t_max, length.out = 400)
  
  H_hat <- M0(tgrid, b, knots) # 计算基准累积风险函数 
  S_hat <- exp(-H_hat) # 计算基准生存函数
  
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2))
  
  plot(
    tgrid, H_hat, type = "l", lwd = 2,
    xlab = "Age", ylab = expression(hat(H)[0](t)),
    main = "Estimated Baseline Cumulative Hazard"
  )
  abline(v = knots, lty = 2)
  
  plot(
    tgrid, S_hat, type = "l", lwd = 2,
    xlab = "Age", ylab = expression(hat(S)[0](t)),
    main = "Estimated Baseline Survival"
  )
  abline(v = knots, lty = 2)
  
  par(oldpar)
}

plot_empirical_baseline(fit_empirical, data_input)

## --------------7. profile likelihood for eta-----------------
# 绘制 eta 和 negloglik 的关系图，显示 profile likelihood 曲线
profile_eta_curve <- function(data, knots, p) {
  eta_grid <- get_eta_grid(data)  # 获取候选的 eta 网格
  obj_vals <- rep(NA, length(eta_grid))  # 初始化存储负对数似然值的向量
  
  for (i in seq_along(eta_grid)) {
    fit_i <- fit_given_eta(data, knots, p, eta_grid[i])  # 给定 eta 拟合模型
    if (!is.null(fit_i)) {
      obj_vals[i] <- fit_i$value  # 存储负对数似然值
    }
  }
  
  data.frame(eta = eta_grid, negloglik = obj_vals)  # 返回 eta 和对应的负对数似然值
}

prof_df <- profile_eta_curve(data_input, fit_empirical$knots, p = 2)

plot(
  prof_df$eta, prof_df$negloglik,
  type = "b", pch = 19,
  xlab = expression(eta),
  ylab = "Profile negative log-likelihood",
  main = expression("Profile likelihood for " * eta)
)
abline(v = fit_empirical$eta, lty = 2, col = 2)  # 标记最优 eta

# 还原成正常年龄尺度再画
prof_df$eta_age <- prof_df$eta * age_sd + age_mean + 48 # 还原到年龄尺度
png("E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/profile_eta_age_2para.png", width = 2000, height = 1600, res = 300)

plot(
  prof_df$eta_age, prof_df$negloglik,
  type = "b", pch = 19,
  xlab = "Age",
  ylab = "Negative log-likelihood",
  main = "Profile likelihood for threshold age",
  
  cex.lab = 1.2,   # 坐标轴标题
  cex.axis = 1.0,  # 坐标刻度
  cex.main = 1.2   # 主标题
)
abline(v = eta_real_age, lty = 2, col = 2)  # 标记最优年龄阈值

dev.off()

## -------------------8.绘制基于模型的 Kaplan-Meier 曲线---------------
library(MASS)

# 还原真实年龄
data_input$age_centered <- data_input$X1 * age_sd + age_mean
data_input$age_real <- data_input$age_centered + 48

# 在给定 theta=(beta,gamma,log b) 与固定 eta 下，计算各组平均生存曲线（一列一组）
.mean_surv_by_groups <- function(theta, eta, knots, tgrid, group_data_list, vars, p) {
  K <- length(theta) - 2 * p
  beta <- theta[1:p]
  gamma <- theta[(p + 1):(2 * p)]
  theta_b <- theta[(2 * p + 1):(2 * p + K)]
  b <- exp(theta_b)
  n_g <- length(group_data_list)
  S_mat <- matrix(NA_real_, nrow = length(tgrid), ncol = n_g)
  for (k in seq_len(n_g)) {
    group_data <- group_data_list[[k]]
    if (nrow(group_data) == 0L) next
    X <- as.matrix(group_data[, vars, drop = FALSE])
    ind <- as.numeric(group_data$X1 > eta)
    psi <- as.vector(X %*% beta + (X %*% gamma) * ind)
    psi <- pmin(pmax(psi, -20), 20)
    S_mat[, k] <- colMeans(sapply(tgrid, function(t) {
      exp(-M0(t, b, knots) * exp(psi)) # 计算生存函数
    }))
  }
  S_mat
}

# 变点模型平均预测生存曲线：共用绘图逻辑
# 置信带：对 optim 中 Hessian 反演得协方差，对 (beta,gamma,log b) 做参数 bootstrap 分位数带；eta 固定为估计值
model_mean_survival <- function(
  fit_result,
  data_input,
  group_data_list,
  cols,
  ylab,
  main,
  legend_pos,
  outfile,
  width,
  height,
  res,
  cex.legend,
  ci = TRUE,
  ci_boot = 400,
  conf_level = 0.95,
  ribbon_alpha = 0.22
) {
  p <- length(fit_result$beta)
  eta <- fit_result$eta
  knots <- fit_result$knots
  vars <- paste0("X", seq_len(p))
  if (!"age_real" %in% names(data_input)) {
    stop("data_input 需包含列 age_real（请先按 age_centered、age_real 还原）")
  }
  if (!all(vars %in% names(data_input))) {
    stop("data_input 需包含与模型一致的协变量列 X1,...,Xp")
  }
  t_min <- min(data_input$Start, na.rm = TRUE)
  t_max <- max(data_input$Stop, na.rm = TRUE)
  tgrid <- seq(t_min, t_max, length.out = 400) # 构造时间轴
  n_g <- length(group_data_list)
  par_hat <- fit_result$optim$par
  S_mat <- .mean_surv_by_groups(par_hat, eta, knots, tgrid, group_data_list, vars, p) # 计算生存函数
  lo <- hi <- NULL
  if (isTRUE(ci) && !is.null(fit_result$optim$hessian)) {
    H <- fit_result$optim$hessian
    if (nrow(H) == length(par_hat) && ncol(H) == length(par_hat)) {
      H <- (H + t(H)) / 2
      cov_mat <- tryCatch(
        ginv(H), # hessian矩阵的逆为协方差矩阵
        error = function(e) NULL
      )
      if (!is.null(cov_mat)) {
        draws <- tryCatch(
          mvrnorm(ci_boot, par_hat, cov_mat), # 正态分布中采样，进行参数bootstrap
          error = function(e) NULL
        )
        if (!is.null(draws)) {
          qlo <- (1 - conf_level) / 2
          qhi <- 1 - qlo
          lo <- matrix(NA_real_, nrow = length(tgrid), ncol = n_g)
          hi <- matrix(NA_real_, nrow = length(tgrid), ncol = n_g)
          S_arr <- array(NA_real_, dim = c(length(tgrid), n_g, ci_boot))
          for (bb in seq_len(ci_boot)) { # 每组参数计算生存曲线
            S_arr[, , bb] <- .mean_surv_by_groups(
              draws[bb, ], eta, knots, tgrid, group_data_list, vars, p
            )
          }
          for (k in seq_len(n_g)) { # 对每个时间点t求上下界
            lo[, k] <- apply(S_arr[, k, ], 1, quantile, probs = qlo, na.rm = TRUE)
            hi[, k] <- apply(S_arr[, k, ], 1, quantile, probs = qhi, na.rm = TRUE)
          }
          lo <- pmax(pmin(lo, 1), 0)
          hi <- pmax(pmin(hi, 1), 0)
        }
      }
    }
  }
  ok <- !apply(S_mat, 2, function(x) all(is.na(x)))
  if (!is.null(outfile)) {
    png(outfile, width = width, height = height, res = res)
    on.exit(dev.off(), add = TRUE)
  }
  plot(
    NA, NA,
    xlim = c(t_min, t_max), ylim = c(0, 1),
    xlab = "Follow-up time", ylab = ylab, main = main
  )
  if (isTRUE(ci) && !is.null(lo) && !is.null(hi)) {
    for (k in seq_len(n_g)) {
      if (!ok[k]) next
      polygon( # 画置信区间阴影
        c(tgrid, rev(tgrid)),
        c(lo[, k], rev(hi[, k])),
        col = adjustcolor(cols[k], alpha.f = ribbon_alpha),
        border = NA
      )
    }
  }
  matlines( # 画曲线
    tgrid, S_mat[, ok, drop = FALSE],
    lty = 1, col = cols[ok], lwd = 2
  )
  legend( # 图例
    legend_pos, legend = names(group_data_list)[ok], col = cols[ok], lty = 1, lwd = 2,
    bty = "n", cex = cex.legend
  )
  invisible(list(tgrid = tgrid, S = S_mat, lower = lo, upper = hi))
}
# 按年龄分组
plot_group_km <- function( 
  fit_result,
  data_input,
  eta_real_age,
  main = "Model-based survival by age threshold",
  outfile = NULL,
  width = 2000,
  height = 1600,
  res = 300,
  ci = TRUE,
  ci_boot = 400,
  conf_level = 0.95,
  ribbon_alpha = 0.22
) {
  d <- data_input
  d$group <- ifelse(d$age_real > eta_real_age, "Age > threshold", "Age ≤ threshold")
  d$group <- factor(d$group, levels = c("Age ≤ threshold", "Age > threshold"))
  group_data_list <- list(
    "Age ≤ threshold" = subset(d, group == "Age ≤ threshold"),
    "Age > threshold" = subset(d, group == "Age > threshold")
  )
  cols <- c("#2E86C1", "#E74C3C")
  model_mean_survival(
    fit_result, data_input, group_data_list, cols,
    ylab = "Survival probability", main = main, legend_pos = "topright",
    outfile = outfile, width = width, height = height, res = res, cex.legend = 1,
    ci = ci, ci_boot = ci_boot, conf_level = conf_level, ribbon_alpha = ribbon_alpha
  )
}
# 分年龄+是否移植分组
plot_km_age_transplant <- function(
  fit_result,
  data_input,
  eta_real_age,
  main = "Model-based survival by age threshold and transplant",
  outfile = NULL,
  width = 2000,
  height = 1600,
  res = 300,
  ci = TRUE,
  ci_boot = 400,
  conf_level = 0.95,
  ribbon_alpha = 0.22
) {
  d <- data_input
  age_hi <- d$age_real > eta_real_age
  tx <- d$X2 == 1L
  group_data_list <- list(
    "Age ≤ threshold, No transplant" = d[!age_hi & !tx, , drop = FALSE],
    "Age ≤ threshold, Transplant" = d[!age_hi & tx, , drop = FALSE],
    "Age > threshold, No transplant" = d[age_hi & !tx, , drop = FALSE],
    "Age > threshold, Transplant" = d[age_hi & tx, , drop = FALSE]
  )
  cols <- c("#2E86C1", "#E74C3C", "#27AE60", "#8E44AD")
  model_mean_survival(
    fit_result, data_input, group_data_list, cols,
    ylab = "Survival probability", main = main, legend_pos = "topright",
    outfile = outfile, width = width, height = height, res = res, cex.legend = 0.85,
    ci = ci, ci_boot = ci_boot, conf_level = conf_level, ribbon_alpha = ribbon_alpha
  )
}

# 与 KM 图分层一致的非参数 log-rank 检验（LTRC：Surv(Start, Stop, Event)；rho=0 为标准 log-rank）
run_logrank_km <- function(data_input, eta_real_age, strat_type = c("age", "age_tx", "age_surgery")) {
  strat_type <- match.arg(strat_type)
  if (!"age_real" %in% names(data_input)) {
    stop("run_logrank_km: data_input 需含 age_real")
  }
  d <- data_input
  if (strat_type == "age") {
    d$strat <- ifelse(d$age_real > eta_real_age, "Age > threshold", "Age ≤ threshold")
    d$strat <- factor(d$strat, levels = c("Age ≤ threshold", "Age > threshold"))
  } else if (strat_type == "age_tx") {
    age_hi <- d$age_real > eta_real_age
    tx <- d$X2 == 1L
    d$strat <- interaction(
      ifelse(age_hi, "Age > threshold", "Age ≤ threshold"),
      ifelse(tx, "Transplant", "No transplant"),
      sep = ", "
    )
    d$strat <- factor(d$strat, levels = c(
      "Age ≤ threshold, No transplant",
      "Age ≤ threshold, Transplant",
      "Age > threshold, No transplant",
      "Age > threshold, Transplant"
    ))
  } else {
    if (!"X3" %in% names(d)) {
      stop("run_logrank_km: strat_type=age_surgery 需要列 X3")
    }
    age_hi <- d$age_real > eta_real_age
    sg <- d$X3 == 1L
    d$strat <- interaction(
      ifelse(age_hi, "Age > threshold", "Age ≤ threshold"),
      ifelse(sg, "Surgery", "No surgery"),
      sep = ", "
    )
    d$strat <- factor(d$strat, levels = c(
      "Age ≤ threshold, No surgery",
      "Age ≤ threshold, Surgery",
      "Age > threshold, No surgery",
      "Age > threshold, Surgery"
    ))
  }
  d$strat <- droplevels(d$strat) #去除空的分组
  if (nlevels(d$strat) < 2L) {
    warning("run_logrank_km: 有效分层少于 2 组，跳过")
    return(invisible(NULL))
  }
  survdiff(Surv(Start, Stop, Event) ~ strat, data = d, rho = 0) # 执行log-rank检验
}


# 调用
plot_group_km(
  fit_empirical,
  data_input,
  eta_real_age,
  main = "",
  outfile = "E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/km_age_2para_confi.png"
)

plot_km_age_transplant(
  fit_empirical,
  data_input,
  eta_real_age,
  main = "",
  outfile = "E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/km_age_transplant_2para_confi.png"
)

cat("\n========== Log-rank（两变量数据，分层与 KM 图一致）==========\n")
print(run_logrank_km(data_input, eta_real_age, "age"))
cat("\n")
print(run_logrank_km(data_input, eta_real_age, "age_tx"))


## ---------------2. 数据整理--------------

dat_raw <- heart

# 数据字段
# start : 左截断时间
# stop  : 观测结束时间
# event : 是否死亡 (1=死亡)
# age   : 年龄
# transplant : 是否接受移植

# 构造 Event
dat_raw$Event <- dat_raw$event

# 构造 Start / Stop
dat_raw$Start <- dat_raw$start
dat_raw$Stop  <- dat_raw$stop

# 变点变量 X1 = 年龄（标准化）
dat_raw$X1 <- as.numeric(scale(dat_raw$age))

# 协变量 X2 = 是否移植
dat_raw$X2 <- as.numeric(as.character(dat_raw$transplant))
# 协变量 X3 = 既往搭桥手术情况
dat_raw$X3 <- as.numeric(as.character(dat_raw$surgery))

# 最终分析数据
data_input <- dat_raw[, c("Start", "Stop", "Event", "X1", "X2","X3")]

# 清理异常值
data_input <- subset(
  data_input,
  is.finite(Start) & is.finite(Stop) &
    is.finite(Event) & is.finite(X1) &
    is.finite(X2)  &
    is.finite(X3) &
    Stop > Start
)

cat("样本量:", nrow(data_input), "\n")
cat("事件数:", sum(data_input$Event), "\n")

summary(data_input)

## -------------3. 标准 Cox 模型（左截断 + 右删失）作对照-----------
library(survival)

cox_fit <- coxph(Surv(Start, Stop, Event) ~ X1 + X2 + X3, data = data_input)
summary(cox_fit)

## -------------------4. 拟合我的模型（三个变量）------------
p <- 3

# 选择最优分割点
knots <- select_knots(data_input, p = 3)

# 拟合模型
fit_empirical <- fit_piecewise(data_input, p = p, knots = knots)

# p值
library(MASS)
hessian <- fit_empirical$optim$hessian
hessian <- (hessian + t(hessian)) / 2   # 强制对称
cov_mat <- ginv(hessian)
se <- sqrt(diag(cov_mat))
est <- fit_empirical$par[1:length(se)]
z <- est / se
p_value <- 2 * pnorm(-abs(z))
CI_lower <- est - 1.96 * se
CI_upper <- est + 1.96 * se

cat("========== 模型拟合结果 ==========\n")
cat("eta =", fit_empirical$eta, "\n")
cat("knots =", paste(round(fit_empirical$knots, 4), collapse = ", "), "\n")
cat("beta =", paste(round(fit_empirical$beta, 4), collapse = ", "), "\n")
cat("gamma =", paste(round(fit_empirical$gamma, 4), collapse = ", "), "\n")
cat("b =", paste(round(fit_empirical$b, 4), collapse = ", "), "\n")
cat("a =", paste(round(fit_empirical$a, 4), collapse = ", "), "\n")
cat("logLik =", round(fit_empirical$logLik, 4), "\n")
cat("AIC =", round(fit_empirical$AIC, 4), "\n")

param_names <- c(
  paste0("beta",1:p),
  paste0("gamma",1:p),
  paste0("logb",1:(length(knots)-1))
)

result_table <- data.frame(
  Parameter = param_names,
  Estimate = est,
  SE = se,
  Z = z,
  P_value = p_value,
  CI_lower = CI_lower,
  CI_upper = CI_upper
)

print(result_table)

write.xlsx(
  result_table,
  file = "E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/result_table_3para.xlsx",
  rowNames = FALSE
)

## -------------------4b. 变点存在性检验（三变量：SUP 置换，H0: γ=0）-----------------
cp_test_config_3 <- list(
  p = 3L,
  B_perm = 1000L,
  k = 20L,
  eta.trim = 0.01,
  alpha = 0.05
)
set.seed(20250323)
cp_test_3 <- single_cp_test_score(data_input, cp_test_config_3)
cp_pval_3 <- mean(cp_test_3$sup_perm >= cp_test_3$SUP_obs)
cp_crit_3 <- as.numeric(cp_test_3$crit)
cp_reject_3 <- cp_test_3$SUP_obs > cp_crit_3

cat("\n========== 变点存在性检验（三变量模型，SUP 置换）==========\n")
cat("H0: 无变点（γ=0）；变点协变量 X1 在置换下可交换\n")
cat("SUP_obs =", round(cp_test_3$SUP_obs, 6), "\n")
cat("置换 p 值 Pr(T_perm >= SUP_obs) =", round(cp_pval_3, 6), "\n")
cat("临界值 alpha=0.05 =", round(cp_crit_3, 6), "\n")
cat("是否拒绝 H0 (alpha=0.05):", ifelse(cp_reject_3, "是", "否"), "\n")

cp_test_table_3 <- data.frame(
  Model = "3 covariates (age, transplant, surgery)",
  SUP_obs = cp_test_3$SUP_obs,
  p_value_perm = cp_pval_3,
  crit_0.05 = cp_crit_3,
  reject_H0 = cp_reject_3,
  B_perm = cp_test_config_3$B_perm,
  k_grid = cp_test_config_3$k,
  stringsAsFactors = FALSE
)

write.xlsx(
  cp_test_table_3,
  file = "E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/cp_test_table_3para.xlsx"
)

## -------------5. 将 eta 还原到原始年龄尺度-----------
age_mean <- mean(dat_raw$age, na.rm = TRUE)
age_sd   <- sd(dat_raw$age, na.rm = TRUE)

eta_std <- fit_empirical$eta # 变点（标准化尺度）
eta_age <- eta_std * age_sd + age_mean # # 还原到 age 变量尺度

eta_real_age <- eta_age + 48 # 还原真实年龄

cat("变点（标准化尺度）eta =", round(eta_std, 4), "\n")
cat("变点（进入年龄原始尺度）eta_age =", round(eta_age, 4), "\n")
cat("变点（真实年龄）eta_real_age =", round(eta_real_age, 4), "\n")

# 还原真实年龄（与第 8 节一致，供 plot_group_km / plot_km_age_transplant）
data_input$age_centered <- data_input$X1 * age_sd + age_mean
data_input$age_real <- data_input$age_centered + 48

## ---------------6. 作图：估计的基准累积风险函数与基准生存函数-------------
plot_empirical_baseline <- function(fit_result, data_input) {
  b <- fit_result$b
  knots <- fit_result$knots
  
  t_min <- min(data_input$Start, na.rm = TRUE)
  t_max <- max(data_input$Stop, na.rm = TRUE)
  tgrid <- seq(t_min, t_max, length.out = 400)
  
  H_hat <- M0(tgrid, b, knots) # 计算基准累积风险函数 
  S_hat <- exp(-H_hat) # 计算基准生存函数
  
  oldpar <- par(no.readonly = TRUE)
  par(mfrow = c(1, 2))
  
  plot(
    tgrid, H_hat, type = "l", lwd = 2,
    xlab = "Age", ylab = expression(hat(H)[0](t)),
    main = "Estimated Baseline Cumulative Hazard"
  )
  abline(v = knots, lty = 2)
  
  plot(
    tgrid, S_hat, type = "l", lwd = 2,
    xlab = "Age", ylab = expression(hat(S)[0](t)),
    main = "Estimated Baseline Survival"
  )
  abline(v = knots, lty = 2)
  
  par(oldpar)
}

plot_empirical_baseline(fit_empirical, data_input)

## ----------------7. profile likelihood for eta---------------
# 绘制 eta 和 negloglik 的关系图，显示 profile likelihood 曲线
profile_eta_curve <- function(data, knots, p) {
  eta_grid <- get_eta_grid(data)  # 获取候选的 eta 网格
  obj_vals <- rep(NA, length(eta_grid))  # 初始化存储负对数似然值的向量
  
  for (i in seq_along(eta_grid)) {
    fit_i <- fit_given_eta(data, knots, p, eta_grid[i])  # 给定 eta 拟合模型
    if (!is.null(fit_i)) {
      obj_vals[i] <- fit_i$value  # 存储负对数似然值
    }
  }
  
  data.frame(eta = eta_grid, negloglik = obj_vals)  # 返回 eta 和对应的负对数似然值
}

prof_df <- profile_eta_curve(data_input, fit_empirical$knots, p = 3)

plot(
  prof_df$eta, prof_df$negloglik,
  type = "b", pch = 19,
  xlab = expression(eta),
  ylab = "Profile negative log-likelihood",
  main = expression("Profile likelihood for " * eta)
)
abline(v = fit_empirical$eta, lty = 2, col = 2)  # 标记最优 eta

# 还原成正常年龄尺度再画
prof_df$eta_age <- prof_df$eta * age_sd + age_mean + 48 # 还原到年龄尺度
png("E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/profile_eta_age_3para.png", width = 2000, height = 1600, res = 300)

plot(
  prof_df$eta_age, prof_df$negloglik,
  type = "b", pch = 19,
  xlab = "Age",
  ylab = "Negative log-likelihood",
  main = "Profile likelihood for threshold age",
  
  cex.lab = 1.2,   # 坐标轴标题
  cex.axis = 1.0,  # 坐标刻度
  cex.main = 1.2   # 主标题
)
abline(v = eta_real_age, lty = 2, col = 2)  # 标记最优年龄阈值

dev.off()

## -------------------8.绘制基于模型的 Kaplan-Meier 曲线---------------

# 分年龄+是否手术分组（需 data_input 含 X3，与三变量模型一致）
plot_km_age_surgery <- function(
    fit_result,
    data_input,
    eta_real_age,
    main = "Model-based survival by age threshold and surgery",
    outfile = NULL,
    width = 2000,
    height = 1600,
    res = 300,
    ci = TRUE,
    ci_boot = 400,
    conf_level = 0.95,
    ribbon_alpha = 0.22
) {
  d <- data_input
  if (!"X3" %in% names(d)) {
    stop("plot_km_age_surgery: data_input 需包含列 X3")
  }
  age_hi <- d$age_real > eta_real_age
  sg <- d$X3 == 1L
  group_data_list <- list(
    "Age ≤ threshold, No surgery" = d[!age_hi & !sg, , drop = FALSE],
    "Age ≤ threshold, Surgery" = d[!age_hi & sg, , drop = FALSE],
    "Age > threshold, No surgery" = d[age_hi & !sg, , drop = FALSE],
    "Age > threshold, Surgery" = d[age_hi & sg, , drop = FALSE]
  )
  cols <- c("#2E86C1", "#E74C3C", "#27AE60", "#8E44AD")
  model_mean_survival(
    fit_result, data_input, group_data_list, cols,
    ylab = "Survival probability", main = main, legend_pos = "topright",
    outfile = outfile, width = width, height = height, res = res, cex.legend = 0.85,
    ci = ci, ci_boot = ci_boot, conf_level = conf_level, ribbon_alpha = ribbon_alpha
  )
}

# 三变量模型：变点模型平均预测生存曲线 + 参数 bootstrap 置信带（与两变量脚本中用法一致）
plot_group_km(
  fit_empirical,
  data_input,
  eta_real_age,
  main = "",
  outfile = "E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/km_age_3para_confi.png"
)

plot_km_age_transplant(
  fit_empirical,
  data_input,
  eta_real_age,
  main = "",
  outfile = "E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/km_age_transplant_3para_confi.png"
)

plot_km_age_surgery(
  fit_empirical,
  data_input,
  eta_real_age,
  main = "",
  outfile = "E:/BNU/BA4/毕业论文/LTRC-changepoint/real_data/km_age_surgery_3para_confi.png"
)

cat("\n========== Log-rank（三变量数据，分层与 KM 图一致）==========\n")
print(run_logrank_km(data_input, eta_real_age, "age"))
cat("\n")
print(run_logrank_km(data_input, eta_real_age, "age_tx"))
cat("\n")
print(run_logrank_km(data_input, eta_real_age, "age_surgery"))
