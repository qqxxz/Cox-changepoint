SIM_CONFIG <- list(

  ## ========= Monte Carlo =========
  B = 100,
  n = 300,
  Distribution = "Exp",
  truncation = 0.1, # 10%, 30%
  censor = 0.2,  # 20%, 40%

  ## ========= 分布参数 =========
  Lambda = 0.3,
  V = 0.5,                 # Weibull shape
  Coeff = list(a =0.3, b = 3.0),

  ## ========= 真值参数 =========
  p = 2,
  eta = 2,  # 变点参数
  Beta = c(-1,-1), # 变点前的回归系数 (β₁, β₂)
  Gamma = c(1.5,1),  # 变点后的效应变化量  2 1 ➡️ 1.5 1 调低一点

  ## ========= 协变量分布 =========
  x1.mean = 2,
  x1.sd   = 1,  # 小一点 1.5 ➡️ 1
  x2.mean = 0,
  x2.sd   = 1
)
