################################ Description ########################################################
# This function generates dataset with time-fixed covariates for changepoint Cox PH model
# 生成时间不变协变量的变点Cox比例风险模型数据集
# 变点模型: h(t|W) = h₀(t)exp(β'X + γ'X̃·I(U>η))
# 其中 W=(X',U)'，X̃=(1,X')，U 为变点协变量，η 为变点参数

# 使用方法示例:
# result = TimeindepLTRC_gnrt_ChangepointPH(N=300, Distribution="WI", 
#                                            truncation.percent=0.1, censor.percent=0.2)
#####################################################################################################

# 根据协变量、分布类型和模型参数生成生存时间（Time）和风险水平（Xi）
# 变点Cox模型：h(t|W) = h₀(t)exp(β'X + γ'X̃·I(U>η))
Time_gnrt_ChangepointPH <- function(data, Distribution, Coeff, eta){
  u_cp = data$U        # 变点协变量 U
  x1 = data$X1         # 回归协变量 X1
  u_rand = data$U_rand # 均匀随机数

  Beta = Coeff$Beta
  Gamma = Coeff$Gamma
  Lambda = Coeff$Lambda
  V = Coeff$V
  Alpha = Coeff$Alpha

  x_vec = x1
  x_tilde = c(1, x1)

  if (u_cp > eta) {
    Param = exp(sum(Beta * x_vec) + sum(Gamma * x_tilde))
  } else {
    Param = exp(sum(Beta * x_vec))
  }

  if (Distribution == "Exp"){
    # 指数分布: S(t) = exp(-λt·Param)
    TT = -log(u_rand) / Lambda / Param
  } else if (Distribution == "WD"){
    # Weibull递减风险: S(t) = exp(-λt^V·Param)
    TT = (-log(u_rand) / Lambda / Param)^(1/V)
  } else if (Distribution == "WI"){
    # Weibull递增风险: S(t) = exp(-λt^V·Param)
    TT = (-log(u_rand) / Lambda / Param)^(1/V)
  } else if (Distribution == "Gtz"){
     # Gompertz分布
    TT = Alpha * (-log(u_rand)) / Lambda / Param
    TT = log(TT + 1) / Alpha
  } else if (Distribution == "Quadratic"){
    # 累积基线风险为时间的二次函数: Λ₀(t) = a·t²
    a = Coeff$a
    TT = sqrt(-log(u_rand) / a / Param)
  } else if (Distribution == "QuadraticLinear") {
    # Λ0(t) = a t^2 + b t
    a = Coeff$a
    b = Coeff$b
    H = -log(u_rand) / Param
    D = b^2 + 4 * a * H
    TT = (-b + sqrt(D)) / (2 * a)
    } else if (Distribution == "PiecewiseWeibull"){
    # 分段威布尔分布
    alpha1 = Coeff$alpha1
    lambda1 = Coeff$lambda1
    alpha2 = Coeff$alpha2
    lambda2 = Coeff$lambda2
    t_star = Coeff$t_star

    H = -log(u_rand) / Param

    if (H < lambda1 * t_star^alpha1) {
      TT = (H / lambda1)^(1/alpha1)
    } else {
      H1 = lambda1 * t_star^alpha1
      H2 = H - H1
      TT = ((H2 / lambda2) + t_star^alpha2)^(1/alpha2)
    }
  } else {
    stop("Wrong distribution is given.")
  }

  result = list(Time = TT, Xi = Param)
  return(result)
}

# 使用uniroot求解删失参数以达到目标删失率
adjust_censor_rate <- function(target_rate, Data, censor.func){
  # 使用uniroot求解删失参数
  # 目标：找到rate使得实际删失率 = 目标删失率
  
  # 检查输入数据有效性
  if(is.null(Data) || nrow(Data) == 0){
    stop("Data is empty or NULL")
  }
  if(any(is.na(Data$Y)) || any(is.na(Data$Start))){
    stop("Data contains NA values in Y or Start columns")
  }
  if(any(Data$Y <= 0) || any(Data$Start < 0)){
    stop("Data contains invalid values: Y must be > 0, Start must be >= 0")
  }
  N = nrow(Data)
  
  # 保存当前随机种子状态
  old_seed = NULL
  if(exists(".Random.seed", envir = .GlobalEnv)){
    old_seed = .Random.seed
  }
  
  # 定义目标函数：f(rate) = 实际删失率(rate) - 目标删失率
  # 需要找到使 f(rate) = 0 的 rate 值
  objective_func <- function(rate){
    # 检查rate是否有效
    if(is.na(rate) || is.infinite(rate) || rate <= 0){
      return(1e6)  # 返回一个很大的值，表示无效
    }
    
    # 保存当前随机种子状态
    old_seed_func = NULL
    if(exists(".Random.seed", envir = .GlobalEnv)){
      old_seed_func = .Random.seed
    }
    
    result = tryCatch({
      # 使用rate的哈希值作为种子，确保相同rate得到相同结果
      # 将rate转换为整数作为种子的一部分
      seed_value = as.integer(rate * 1e6) %% 1e6 + 1
      set.seed(seed_value)
      
      # 生成删失时间
      Censor.time = censor.func(N, rate)
      # 检查Censor.time是否有效
      if(any(is.na(Censor.time)) || any(is.infinite(Censor.time)) || any(Censor.time < 0)){
        return(1e6)  # 返回一个很大的值，表示无效
      }
      
      # 计算实际删失率
      # 使用na.rm=TRUE来处理可能的NA值
      actual_rate = mean(Data$Y > (Censor.time + Data$Start), na.rm = TRUE)
      # 检查actual_rate是否有效
      if(is.na(actual_rate) || is.infinite(actual_rate)){
        return(1e6)  # 返变点前的回归系数回一个很大的值，表示无效
      }
      
      # 返回差值
      actual_rate - target_rate
    }, error = function(e){
      # 如果出错，返回一个很大的值
      warning(paste("Error in objective_func:", e$message))
      return(1e6)
    })
    
    # 恢复随机种子状态（只在之前存在时才恢复）
    if(!is.null(old_seed_func)){
      .Random.seed <<- old_seed_func
    }
    
    return(result)
  }
  
  # 确定搜索区间
  # rate越小，删失时间越大，删失率越低
  # rate越大，删失时间越小，删失率越高
  
  rate_low = 1e-6
  rate_high = 100
  
  # 测试端点以确保函数值异号
  f_low = objective_func(rate_low)
  f_high = objective_func(rate_high)
  
  # 检查函数值是否有效
  if(is.na(f_low) || is.infinite(f_low)){
    stop("objective_func returned invalid value at lower bound")
  }
  if(is.na(f_high) || is.infinite(f_high)){
    stop("objective_func returned invalid value at upper bound")
  }
  
  # 如果端点值同号，调整区间
  if(!is.na(f_low) && !is.na(f_high) && f_low * f_high > 0){
    if(f_low > 0 && f_high > 0){
      # 都为正，说明删失率太高，需要更小的rate
      rate_high = rate_low
      rate_low = 1e-8
    } else if(f_low < 0 && f_high < 0){
      # 都为负，说明删失率太低，需要更大的rate
      rate_low = rate_high
      rate_high = 1000
    }
    
    # 再次测试
    f_low = objective_func(rate_low)
    f_high = objective_func(rate_high)
    
    # 再次检查
    if(is.na(f_low) || is.infinite(f_low) || is.na(f_high) || is.infinite(f_high)){
      stop("objective_func returned invalid value after interval adjustment")
    }
  }
  
  # 使用uniroot求解
  result = uniroot(objective_func, 
                   interval = c(rate_low, rate_high),
                   extendInt = "yes",
                   tol = 1e-4,
                   maxiter = 100)
  
  # 恢复随机种子状态
  if(!is.null(old_seed)){
    .Random.seed <<- old_seed
  }
  
  return(result$root)
}

# 生成时间不变协变量的变点Cox比例风险模型数据集
TimeindepLTRC_gnrt_ChangepointPH <- function(N = 300, 
                                             Distribution = c("Exp", "WI", "WD", "Gtz", "Quadratic","QuadraticLinear", "PiecewiseWeibull"), 
                                             eta = 2,  # 变点参数
                                             Beta = c(0.5),  # 变点前的回归系数 (β₁, β₂)
                                             Gamma = c(0, 0),  # 变点后的效应变化量 (γ₀, γ₁, γ₂)
                                             truncation.percent = 0.1,  # 左截断百分比 (10% 或 30%)
                                             censor.percent = 0.2,  # 右删失百分比 (20% 或 40%)
                                             u.mean = 2,   # 变点协变量 U 的均值
                                             u.sd = 1.5,   # 变点协变量 U 的标准差
                                             x1.mean = 0,  # 回归协变量 X1 的均值
                                             x1.sd = 1,    # 回归协变量 X1 的标准差
                                             weibull.shape = 0.5,  # Weibull形状参数 (0.5 或 3.0)
                                             adjust.censor = TRUE){  # 是否调整删失参数以达到目标删失率
  
  # 初始化数据框
  Data <- as.data.frame(matrix(NA, N * 3, 10))
  names(Data) <- c("I", "ID", "U", "X1", "Start", "Stop", "C", "Event", "Y", "Xi")
  Data$C <- 0
  Count = 0
  
  # 根据分布类型设置参数
  Coeff = NULL
  if (Distribution == "Exp"){
    Lambda = 0.3
    Alpha = 0
    V = 0
  } else if (Distribution == "WD"){
    # Weibull递减风险
    Lambda = 0.3
    V = 0.5 # 0.5 或 3.0
    Alpha = 0
  } else if (Distribution == "WI"){
    # Weibull递增风险
    Lambda = 0.3
    V = weibull.shape  # 0.5 或 3.0
    Alpha = 0
  } else if (Distribution == "Gtz"){
    Alpha = 0.1
    Lambda = 0.3
    V = 0
  } else if (Distribution == "Quadratic"){
    # 累积基线风险为时间的二次函数: Λ₀(t) = a·t²
    Coeff$a = 0.01
    Lambda = 0
    Alpha = 0
    V = 0
  } else if (Distribution == "QuadraticLinear"){
    # Λ0(t) = a t^2 + b t
    Coeff$a = 0.3
    Coeff$b = 3.0
    Lambda = 0
    Alpha = 0
    V = 0
  } else if (Distribution == "PiecewiseWeibull"){
    # 分段威布尔分布
    Coeff$alpha1 = 0.5
    Coeff$lambda1 = 0.3
    Coeff$alpha2 = 2.0
    Coeff$lambda2 = 0.3
    Coeff$t_star = 2.0
    Lambda = 0
    Alpha = 0
    V = 0
  } else {
    stop("Wrong distribution is given.")
  }
  
  Coeff$Alpha = Alpha
  Coeff$Beta = Beta
  Coeff$Gamma = Gamma
  Coeff$Lambda = Lambda
  Coeff$V = V
  
  # 预生成协变量：U 为变点变量，X1 为回归协变量
  max_samples = N * 5
  u_all = rnorm(max_samples, mean = u.mean, sd = u.sd)
  x1_all = rnorm(max_samples, mean = x1.mean, sd = x1.sd)
  
  # 第一步：生成所有候选样本（不考虑左截断）
  candidate_data = list()
  candidate_times = numeric(max_samples)
  
  for(i in 1:max_samples){
    data_temp = list(U = u_all[i], X1 = x1_all[i], U_rand = runif(1))
    ret_temp = Time_gnrt_ChangepointPH(data = data_temp, Distribution = Distribution, 
                                      Coeff = Coeff, eta = eta)
    candidate_times[i] = ret_temp$Time
    candidate_data[[i]] = list(U = u_all[i], X1 = x1_all[i],
                               Time = ret_temp$Time, Xi = ret_temp$Xi)
  }
  
  # 第二步：根据左截断百分比确定截断时间阈值
  truncation_threshold = quantile(candidate_times, truncation.percent)
  
  # 第三步：筛选满足左截断条件的样本
  valid_idx = which(candidate_times > truncation_threshold)
  if(length(valid_idx) < N){
    stop(paste("Not enough valid samples after left truncation. Got", 
               length(valid_idx), "but need", N))
  }
  
  # 随机选择N个有效样本
  selected_idx = sample(valid_idx, N)
  
  # 第四步：为选中的样本生成左截断时间和右删失时间
  for(i in 1:N){
    idx = selected_idx[i]
    data_item = candidate_data[[idx]]
    
    # 生成左截断时间（在[0, truncation_threshold]范围内）
    L = runif(1, 0, truncation_threshold)
    
    Data[i, "U"] = data_item$U
    Data[i, "X1"] = data_item$X1
    Data[i, "Start"] = L
    Data[i, "Y"] = data_item$Time
    Data[i, "Xi"] = data_item$Xi
  }
  
    # 第五步：生成右删失时间
  # 只使用前N行有效数据，并确保没有NA值
  Data_valid = Data[1:N, ]
  
  # 检查Data_valid中是否有NA值
  if(any(is.na(Data_valid$Y)) || any(is.na(Data_valid$Start))){
    # 如果还有NA值，说明填充有问题，只保留有效行
    valid_rows = !is.na(Data_valid$Y) & !is.na(Data_valid$Start) & 
                 Data_valid$Y > 0 & Data_valid$Start >= 0
    Data_valid = Data_valid[valid_rows, ]
    if(nrow(Data_valid) == 0){
      stop("No valid data rows after filtering NA values")
    }
  }
  
  # 使用指数分布生成删失时间
  mean_survival = mean(Data_valid$Y, na.rm = TRUE)
  if(is.na(mean_survival) || mean_survival <= 0){
    stop("Invalid mean survival time")
  }
  initial_censor_rate = 1 / (mean_survival * (1/censor.percent - 1))
  
  if(adjust.censor && censor.percent > 0 && censor.percent < 1){
    # 使用uniroot调整删失参数以达到目标删失率
    censor.func = function(n, rate) rexp(n, rate = rate)
    
    # 使用uniroot求解
    final_rate = adjust_censor_rate(target_rate = censor.percent, 
                                    Data = Data_valid, 
                                    censor.func = censor.func)
    
    # 生成最终的删失时间
    Censor.time = rexp(N, rate = final_rate)
  } else {
    # 不调整或删失率为0或1
    Censor.time = rexp(N, rate = initial_censor_rate)
  }
  
  # 第六步：确定观测时间和事件指示符
  for(i in 1:N){
    if(Data[i, "Y"] <= Censor.time[i] + Data[i, "Start"]){
      Data[i, "Stop"] = Data[i, "Y"]
      Data[i, "Event"] = 1
    } else {
      Data[i, "Stop"] = Censor.time[i] + Data[i, "Start"]
      Data[i, "Event"] = 0
      Data[i, "C"] = 1
    }
  }
  
  # 只保留前N行
  Data = Data[1:N, ]
  
  Data$I = seq(1, N) 
  Data$ID = seq(1, N) 
  
  # 计算实际截断和删失率
  actual_truncation_rate = 1 - length(valid_idx) / max_samples
  actual_censor_rate = mean(Data$Event == 0)
  
  RES = NULL
  RES$Data <- Data
  RES$Info = list(Set = "ChangepointPH", 
                  Coeff = Coeff, 
                  Dist = Distribution,
                  eta = eta,
                  truncation.percent = truncation.percent,
                  censor.percent = censor.percent,
                  actual_truncation_rate = actual_truncation_rate,
                  actual_censor_rate = actual_censor_rate)
  return(RES)
}