# 加载必要的函数
# %%
source("cml_firendly/样本生成函数/三联体家庭结构生成函数.R") # 加载数据生成函数 library(dplyr) library(MASS) 里面加载了这些包
source("cml_firendly/FGWAS 优化版本/FGWAS 函数.R") # 加载 FGWAS 函数

# %%
# 生成数据
data <- generate_multiple_datasets_v3()
data_1_full <- data[[1]] # 完整的暴露数据
data_2_full <- data[[2]] # 完整的结局数据 (包含 snp_type 列)

# 进行FGWAS处理
hat_expose <- fgwas_for_data_optimized(data_1_full)
hat_outcome <- fgwas_for_data_optimized(data_2_full)

# 处理成可以进行 MR 的样子
hat_expose_trMR <- fgwas_to_mr(hat_expose)
hat_outcome_trMR <- fgwas_to_mr(hat_outcome)
