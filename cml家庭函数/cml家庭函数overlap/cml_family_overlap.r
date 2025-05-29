# %% 导入数据生成函数以及处理函数
source("样本生成函数/三联体家庭结构生成函数.R") # 加载数据生成函数 library(dplyr) library(MASS) 里面加载了这些包
source("FGWAS 优化版本/FGWAS 函数.R")
# %% 生成数据然后尝试构建overlap版
phase_two_data_full <- generate_multiple_datasets_v3(n = 5)
phase_one_data <- phase_two_data_full[[1]]
phase_two_data <- phase_two_data_full[[2]]

hat_expose <- fgwas_for_data_optimized(phase_one_data)
hat_outcome <- fgwas_for_data_optimized(phase_two_data)

cor_matrix <- matrix(0, ncol = 2, nrow = 2)
hat_expose$Sigma_inv
sigma_c_i <- matrix(0, ncol = 4, nrow = 4)
# 把新构建的矩阵存入一个list
