# 变点Cox比例风险模型模拟数据生成

## 概述

本代码用于生成**变点Cox比例风险模型**的左截断右删失（LTRC）模拟数据。

## 模型设定

### 变点Cox比例风险模型

风险函数：
```
h(t|x) = h₀(t) · exp(β'x + γ'x̃·I(x₁>η))
```

其中：
- `x₁` 是变点协变量，服从 N(2, 1.5²)
- `x₂` 是回归协变量，服从 N(0, 1)
- `η = 2` 是变点参数
- `β = (-1, 1.5, 0.5)` 是变点前的回归系数向量
- `γ` 是变点后的效应变化量向量
- `x = (1, x₁, x₂)` 是完整设计向量
- `x̃ = (1, x₂)` 是剔除x₁后的设计向量

## 文件说明

### 1. `TimeindepLTRC_gnrt_ChangepointPH_v2.R`

主要的数据生成函数，包含：

- **`Time_gnrt_ChangepointPH()`**: 根据协变量和模型参数生成生存时间
- **`adjust_censor_rate()`**: 调整删失参数以达到目标删失率
- **`TimeindepLTRC_gnrt_ChangepointPH()`**: 主函数，生成完整的LTRC数据集

### 2. `example_changepoint_simulation.R`

使用示例和批量模拟实验框架。

## 主要修改点（相对于原文件）

### 1. 模型结构修改

**原文件** (`TimeindepLTRC_gnrt_PH.R`):
- 使用线性或交互项模型：`exp(β'x)` 或复杂的交互项

**修改后**:
- 实现变点模型：`exp(β'x + γ'x̃·I(x₁>η))`
- 支持变点协变量x₁和变点参数η

### 2. 协变量生成

**原文件**:
- 多个协变量（X1-X20），分布多样

**修改后**:
- 简化为两个协变量：
  - `x₁ ~ N(2, 1.5²)` (变点协变量)
  - `x₂ ~ N(0, 1)` (回归协变量)

### 3. 基准累积风险函数

**原文件**:
- 支持 Exp, WD, WI, Gtz 分布

**修改后**:
- 保留原有分布
- 新增：
  - **Quadratic**: 累积基线风险为时间的二次函数 `Λ₀(t) = a·t²`
  - **PiecewiseWeibull**: 分段威布尔分布

### 4. 左截断和右删失控制

**原文件**:
- 左截断：使用固定的 `truncation` 参数
- 右删失：使用 `censor.rate` (0, 1, 2) 控制

**修改后**:
- 左截断：通过分位数方法精确控制左截断百分比（10% 或 30%）
- 右删失：通过调整删失时间分布参数精确控制右删失百分比（20% 或 40%）
- 添加 `adjust_censor_rate()` 函数自动调整删失参数

### 5. 参数设置

**原文件**:
- 参数硬编码在函数内部

**修改后**:
- 参数可通过函数参数灵活设置
- 支持实验设计中的各种参数组合

## 使用方法

### 基本使用

```r
# 加载函数
source("TimeindepLTRC_gnrt_ChangepointPH_v2.R")

# 生成数据
result = TimeindepLTRC_gnrt_ChangepointPH(
  N = 300,                    # 样本量
  Distribution = "WI",        # 分布类型
  eta = 2,                    # 变点参数
  Beta = c(-1, 1.5, 0.5),    # 变点前回归系数
  Gamma = c(0, 0, 0),        # 变点后效应变化量
  truncation.percent = 0.1,  # 左截断百分比 (10%)
  censor.percent = 0.2,      # 右删失百分比 (20%)
  weibull.shape = 0.5        # Weibull形状参数
)

# 查看数据
head(result$Data)
print(result$Info)
```

### 支持的分布类型

1. **"Exp"**: 指数分布
2. **"WD"**: Weibull递减风险 (V < 1)
3. **"WI"**: Weibull递增风险 (V > 1)
4. **"Gtz"**: Gompertz分布
5. **"Quadratic"**: 累积基线风险为时间的二次函数
6. **"PiecewiseWeibull"**: 分段威布尔分布

### 实验设计参数

根据你的研究设计，可以设置：

```r
# 样本量
sample_sizes = c(300, 500)

# 左截断百分比
truncation_percents = c(0.1, 0.3)  # 10%, 30%

# 右删失百分比
censor_percents = c(0.2, 0.4)  # 20%, 40%

# Weibull形状参数
weibull_shapes = c(0.5, 3.0)
```

## 输出结构

函数返回一个列表，包含：

- **`Data`**: 数据框，包含以下列：
  - `I`, `ID`: 样本标识
  - `X1`, `X2`: 协变量
  - `Start`: 左截断时间（进入研究时间）
  - `Stop`: 观测结束时间
  - `T`: 真实生存时间
  - `Event`: 事件指示符（1=发生事件，0=删失）
  - `C`: 删失指示符
  - `Xi`: 风险水平

- **`Info`**: 信息列表，包含：
  - `Set`: 模型类型
  - `Coeff`: 模型系数
  - `Dist`: 分布类型
  - `eta`: 变点参数
  - `truncation.percent`: 目标左截断百分比
  - `censor.percent`: 目标右删失百分比
  - `actual_truncation_rate`: 实际左截断率
  - `actual_censor_rate`: 实际右删失率

## 注意事项

1. **左截断实现**：代码通过生成大量候选样本，然后根据分位数筛选满足左截断条件的样本。如果生成的样本数不足，会给出警告。

2. **右删失调整**：默认情况下 `adjust.censor = TRUE`，函数会自动调整删失参数以达到目标删失率。如果关闭此选项，实际删失率可能与目标值有偏差。

3. **计算时间**：由于需要生成候选样本并筛选，生成数据可能需要一些时间，特别是当左截断百分比较高时。

4. **参数验证**：确保 `Beta` 和 `Gamma` 的长度正确（都是3维向量）。

## 批量模拟实验

使用 `example_changepoint_simulation.R` 中的函数进行批量模拟：

```r
source("example_changepoint_simulation.R")

# 运行小规模测试
test_results = run_simulation_batch(
  n_sim = 5,
  N = 300,
  Distribution = "WI",
  weibull.shape = 0.5,
  truncation.percent = 0.1,
  censor.percent = 0.2
)
```

## 与原代码的对应关系

| 原代码功能 | 修改后功能 |
|----------|----------|
| `Time_gnrt_PH()` | `Time_gnrt_ChangepointPH()` |
| `TimeindepLTRC_gnrt_PH()` | `TimeindepLTRC_gnrt_ChangepointPH()` |
| 多个协变量模型 | 变点模型（2个协变量） |
| 固定截断参数 | 百分比控制截断 |
| 固定删失率 | 百分比控制删失 |

## 后续工作建议

1. **参数估计**：在生成数据后，需要实现变点Cox模型的参数估计方法
2. **分割点选择**：实现AIC准则下的分割点选择算法
3. **性能评估**：评估估计量的偏差、方差、覆盖率等
4. **敏感性分析**：分析不同参数设置对结果的影响

