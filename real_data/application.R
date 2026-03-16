source("E:/BNU/BA4/毕业论文/LTRC-changepoint/data/TimeindepLTRC_gnrt_ChangepointPH.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/estimation/estimate.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/code/MC.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/code/config.R")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/code/save.R")
setwd("E:/BNU/BA4/毕业论文/LTRC-changepoint")
source("E:/BNU/BA4/毕业论文/LTRC-changepoint/estimation/plot_baseline.R")


## 1. 读取 Channing House 数据
library(boot)
data("channing")

str(channing)
summary(channing)
head(channing)

## 2. 数据整理

dat_raw <- channing

# 检查字段名
names(dat_raw)

# channing 包含:
# sex   : 性别
# entry : 进入观察年龄（左截断时间）
# exit  : 退出观察年龄（右删失时间）
# time  : 生存时长 = exit - entry
# cens  : 删失指示（1=死亡, 0=删失）

# 构造 Event：1 = 事件发生, 0 = 右删失
dat_raw$Event <- dat_raw$cens

# 构造 Start / Stop
dat_raw$Start <- dat_raw$entry
dat_raw$Stop  <- dat_raw$exit

# 性别变量数值化
# 兼容 sex 为字符或因子
dat_raw$sex_num <- ifelse(as.character(dat_raw$sex) %in% c("Male", "male", "M", "m"), 1, 0)

# 令 X1 = 进入年龄（中心化/标准化）
# 令 X2 = 性别 1=男
dat_raw$X1 <- as.numeric(scale(dat_raw$entry))
dat_raw$X2 <- dat_raw$sex_num

# 最终分析数据
data_input <- dat_raw[, c("Start", "Stop", "Event", "X1", "X2")]

# 去掉异常值/缺失
data_input <- subset(
  data_input,
  is.finite(Start) & is.finite(Stop) & is.finite(Event) &
    is.finite(X1) & is.finite(X2) &
    Stop > Start
)

cat("样本量:", nrow(data_input), "\n")
cat("事件数:", sum(data_input$Event), "\n")
summary(data_input)

## 3. 标准 Cox 模型（左截断 + 右删失）作对照
library(survival)

cox_fit <- coxph(Surv(Start, Stop, Event) ~ X1 + X2, data = data_input)
summary(cox_fit)

## 4. 拟合我的模型
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
  P_value = p_value
)

print(result_table)

## 5. 将 eta 还原到原始年龄尺度
entry_mean <- mean(dat_raw$entry, na.rm = TRUE)
entry_sd   <- sd(dat_raw$entry, na.rm = TRUE)

eta_std <- fit_empirical$eta
eta_age <- eta_std * entry_sd + entry_mean

cat("变点（标准化尺度）eta =", round(eta_std, 4), "\n")
cat("变点（进入年龄原始尺度）eta_age =", round(eta_age, 4), "\n")

## 6. 作图：估计的基准累积风险函数与基准生存函数
plot_empirical_baseline <- function(fit_result, data_input) {
  b <- fit_result$b
  knots <- fit_result$knots
  
  t_min <- min(data_input$Start, na.rm = TRUE)
  t_max <- max(data_input$Stop, na.rm = TRUE)
  tgrid <- seq(t_min, t_max, length.out = 400)
  
  H_hat <- M0(tgrid, b, knots)
  S_hat <- exp(-H_hat)
  
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

## 7. profile likelihood for eta
profile_eta_curve <- function(data, knots, p) {
  eta_grid <- get_eta_grid(data)
  obj_vals <- rep(NA, length(eta_grid))
  
  for (i in seq_along(eta_grid)) {
    fit_i <- fit_given_eta(data, knots, p, eta_grid[i])
    if (!is.null(fit_i)) {
      obj_vals[i] <- fit_i$value
    }
  }
  
  data.frame(eta = eta_grid, negloglik = obj_vals)
}

prof_df <- profile_eta_curve(data_input, fit_empirical$knots, p = 2)

plot(
  prof_df$eta, prof_df$negloglik,
  type = "b", pch = 19,
  xlab = expression(eta),
  ylab = "Profile negative log-likelihood",
  main = expression("Profile likelihood for " * eta)
)
abline(v = fit_empirical$eta, lty = 2, col = 2)

# 还原成年龄尺度再画
prof_df$eta_age <- prof_df$eta * entry_sd + entry_mean

plot(
  prof_df$eta_age, prof_df$negloglik,
  type = "b", pch = 19,
  xlab = "Threshold age at entry",
  ylab = "Profile negative log-likelihood",
  main = "Profile likelihood for threshold age"
)
abline(v = eta_age, lty = 2, col = 2)

