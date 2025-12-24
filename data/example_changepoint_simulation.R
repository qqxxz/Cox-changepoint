############################################################
# 变点Cox比例风险模型模拟实验示例
# 基于 TimeindepLTRC_gnrt_ChangepointPH.R
############################################################

source("E:/BNU/BA4/毕业论文/code/【danqi】ltrc-changepoint-main/data/TimeindepLTRC_gnrt_ChangepointPH.R")
# 参数设置
eta = 2  # 变点参数
Beta = c(-1, 1.5, 0.5)  # 变点前的回归系数 (β₀, β₁, β₂)
Gamma = c(0, 0, 0)  # 变点后的效应变化量（可根据需要调整）

# 协变量分布
x1.mean = 2
x1.sd = 1.5
x2.mean = 0
x2.sd = 1

# 实验设计
sample_sizes = c(300, 500)
truncation_percents = c(0.1, 0.3)  # 10%, 30%
censor_percents = c(0.2, 0.4)  # 20%, 40%

# 基准累积风险函数设置
# (a) 威布尔分布
distributions_weibull = list(
  list(name = "WI", shape = 0.5, label = "Weibull_Increasing_shape0.5"),
  list(name = "WI", shape = 3.0, label = "Weibull_Increasing_shape3.0"),
  list(name = "WD", shape = 0.5, label = "Weibull_Decreasing_shape0.5"),
  list(name = "WD", shape = 3.0, label = "Weibull_Decreasing_shape3.0")
)

# (b) 累积基线风险为时间的二次函数
distributions_quadratic = list(
  list(name = "Quadratic", label = "Quadratic_Baseline")
)

# (c) 分段威布尔分布
distributions_piecewise = list(
  list(name = "PiecewiseWeibull", label = "PiecewiseWeibull")
)

# 合并所有分布
all_distributions = c(distributions_weibull, distributions_quadratic, distributions_piecewise)

# ============================================================
# 单个模拟实验示例
# ============================================================

# 示例1：威布尔分布，n=300, 左截断10%, 右删失20%
result1 = TimeindepLTRC_gnrt_ChangepointPH(
  N = 300,
  Distribution = "WI",
  eta = eta,
  Beta = Beta,
  Gamma = Gamma,
  truncation.percent = 0.1,
  censor.percent = 0.2,
  x1.mean = x1.mean,
  x1.sd = x1.sd,
  x2.mean = x2.mean,
  x2.sd = x2.sd,
  weibull.shape = 0.5,
  adjust.censor = TRUE
)

# 查看结果
print("示例1结果:")
print(paste("实际左截断率:", round(result1$Info$actual_truncation_rate, 3)))
print(paste("实际右删失率:", round(result1$Info$actual_censor_rate, 3)))
print(head(result1$Data))

# ============================================================
# 批量模拟实验函数
# ============================================================

run_simulation_batch = function(n_sim = 100, 
                                N = 300,
                                Distribution = "WI",
                                weibull.shape = 0.5,
                                truncation.percent = 0.1,
                                censor.percent = 0.2,
                                eta = 2,
                                Beta = c(-1, 1.5, 0.5),
                                Gamma = c(0, 0, 0)){
  # 运行n_sim次模拟实验
  results = list()
  
  for(i in 1:n_sim){
    cat("运行第", i, "次模拟...\n")
    
    result = TimeindepLTRC_gnrt_ChangepointPH(
      N = N,
      Distribution = Distribution,
      eta = eta,
      Beta = Beta,
      Gamma = Gamma,
      truncation.percent = truncation.percent,
      censor.percent = censor.percent,
      weibull.shape = weibull.shape,
      adjust.censor = TRUE
    )
    
    results[[i]] = result
    
    # 每10次输出一次进度
    if(i %% 10 == 0){
      cat("已完成", i, "次模拟\n")
    }
  }
  
  return(results)
}

# ============================================================
# 完整实验设计示例
# ============================================================

run_full_experiment = function(n_sim = 100){
  all_results = list()
  exp_count = 0
  
  # 遍历所有实验设置
  for(N in sample_sizes){
    for(trunc_pct in truncation_percents){
      for(censor_pct in censor_percents){
        for(dist_info in all_distributions){
          exp_count = exp_count + 1
          
          # 确定分布参数
          if(dist_info$name %in% c("WI", "WD")){
            shape = dist_info$shape
            dist_name = dist_info$name
          } else {
            shape = 0.5  # 默认值，对于非Weibull分布不使用
            dist_name = dist_info$name
          }
          
          cat("\n========================================\n")
          cat("实验", exp_count, ":\n")
          cat("样本量:", N, "\n")
          cat("左截断率:", trunc_pct, "\n")
          cat("右删失率:", censor_pct, "\n")
          cat("分布:", dist_info$label, "\n")
          cat("========================================\n")
          
          # 运行模拟
          results = run_simulation_batch(
            n_sim = n_sim,
            N = N,
            Distribution = dist_name,
            weibull.shape = shape,
            truncation.percent = trunc_pct,
            censor.percent = censor_pct
          )
          
          # 保存结果
          all_results[[exp_count]] = list(
            N = N,
            truncation.percent = trunc_pct,
            censor.percent = censor_pct,
            distribution = dist_info$label,
            results = results
          )
        }
      }
    }
  }
  
  return(all_results)
}

# ============================================================
# 结果分析函数
# ============================================================

# 计算模拟结果的统计量
analyze_results = function(results_list){
  # results_list 是 run_simulation_batch 的返回值
  
  n_sim = length(results_list)
  
  # 提取关键统计量
  truncation_rates = sapply(results_list, function(x) x$Info$actual_truncation_rate)
  censor_rates = sapply(results_list, function(x) x$Info$actual_censor_rate)
  
  summary_stats = list(
    mean_truncation_rate = mean(truncation_rates),
    sd_truncation_rate = sd(truncation_rates),
    mean_censor_rate = mean(censor_rates),
    sd_censor_rate = sd(censor_rates)
  )
  
  return(summary_stats)
}

# ============================================================
# 使用示例
# ============================================================

# 小规模测试（建议先运行这个）
cat("运行小规模测试...\n")
test_results = run_simulation_batch(
  n_sim = 5,
  N = 300,
  Distribution = "WI",
  weibull.shape = 0.5,
  truncation.percent = 0.1,
  censor.percent = 0.2
)

# 分析结果
test_summary = analyze_results(test_results)
print(test_summary)

# 完整实验（需要较长时间，谨慎运行）
# full_results = run_full_experiment(n_sim = 100)

