SIM_CONFIG <- list(

  ## ========= Monte Carlo =========
  B = 500,
  n = 300,
  Distribution = "Exp",
  truncation = 0.1, # 10%, 30%
  censor = 0.2,  # 20%, 40%

  ## ========= 分布参数 =========
  Lambda = 0.3,
  V = 0.5,                 # Weibull shape
  Coeff = list(a =0.3, b = 3.0),

  ## ========= 真值参数 =========
  p = 1,
  eta = 2,  # 变点参数
  Beta = c(-1),  # 变点前的回归系数 β₁
  Gamma = c(0.9, 0.5),  # 变点后的效应变化量 (γ₀, γ₁)

  ## ========= 协变量分布 =========
  u.mean = 2,   # 变点协变量 U 的均值
  u.sd   = 0.5, # 变点协变量 U 的标准差
  x1.mean = 0,  # 回归协变量 X1 的均值
  x1.sd   = 1,  # 回归协变量 X1 的标准差

  ## ========= M/I 样条 =========
  spline_order = 3L,
  k_candidates = c(3L, 5L, 7L, 9L),

  ## ========= SUP 检验参数 =========
  k = 20,                  # 变点候选网格数 |H|
  eta.trim = 0.1,        # 去除两端极端点比例
  B_perm = 1000,           # 每次检验的置换次数

  ## ========= 显著性水平 =========
  alpha = c(0.10, 0.05, 0.01),

  ## ========= RNG 控制 =========
  base_seed = 123,
  seed_multiplier = 6699,
  exp_id = 0L,
  test_id = 0L,

  ## ========= 项目路径（Windows 并行 cluster 也会用到）=========
  project_dir = "E:/BNU/BA4/Cox-changepoint/code"
)

