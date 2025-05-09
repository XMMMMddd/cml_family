# 加载必要的函数
# %%
source("cml家庭/样本生成函数/三联体家庭结构生成函数.R") # 加载数据生成函数 library(dplyr) library(MASS) 里面加载了这些包
source("cml家庭/FGWAS 优化版本/FGWAS 函数.R") # 加载 FGWAS 函数

# %%
# 生成数据
data <- generate_multiple_datasets_v3(n = 10, num_pleiotropic = 1, beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1, beta_OStoOE_exp = 0.3, mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0.05, mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0.05, mean_beta_OStoOE_out = 0.1, sd_beta_OStoOE_out = 0.05, prop_negative_pleiotropy = 0.5, compatibility_selection_prop = 0, compatibility_selection_geno = "independent", correlation_param = 0.5, compatibility_selection_factor_exp = 0, compatibility_selection_factor_out = 0, crowd_stratification_differences = 0, beta_exp_to_out = 0.4, beta_confounding_exp = 0.2, beta_confounding_out = 0.2, correlation = 0.2, seed = NULL)
data_1_full <- data[[1]] # 完整的暴露数据
data_2_full <- data[[2]] # 完整的结局数据 (包含 snp_type 列)

# 进行FGWAS处理
hat_expose <- fgwas_for_data_optimized(data_1_full)
hat_outcome <- fgwas_for_data_optimized(data_2_full)

# 处理成可以进行 MR 的样子
hat_expose_trMR <- fgwas_to_mr(hat_expose)
hat_outcome_trMR <- fgwas_to_mr(hat_outcome)



a <- cMl_para(
    beta_y_hat = hat_outcome$beta_hat,
    beta_x_hat = hat_expose$beta_hat,
    Sigma_inv_x = hat_expose$Sigma_inv,
    Sigma_inv_y = hat_outcome$Sigma_inv
)

b <- cMl_para_yuanshi_cpp(
    beta_y_hat = hat_outcome$beta_hat,
    beta_x_hat = hat_expose$beta_hat,
    Sigma_inv_x = hat_expose$Sigma_inv,
    Sigma_inv_y = hat_outcome$Sigma_inv
)

c <- cml_family_2_cn_cpp(
    beta_y_hat = hat_outcome$beta_hat,
    beta_x_hat = hat_expose$beta_hat,
    Sigma_inv_x = hat_expose$Sigma_inv,
    Sigma_inv_y = hat_outcome$Sigma_inv
)
a$theta
b$theta
c$theta_weight
