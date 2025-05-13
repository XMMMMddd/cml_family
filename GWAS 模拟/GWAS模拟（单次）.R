# %%

source("样本生成函数/三联体家庭结构生成函数.R") # 加载数据生成函数 library(dplyr) library(MASS) 里面加载了这些包
source("C:/users/Administrator/Desktop/cml_family/FGWAS 优化版本/FGWAS 函数.R") # 加载 FGWAS 函数

# %%

gwas_simulation_once <- function(
    n = 1, num_pleiotropic = 0,
    N_out = 1000,
    # 样本重叠
    overlap_prop = 0,
    # 工具变量的强度
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3,
    # 水平多效性的强度
    mean_beta_FStoOE_out = 0.1,
    sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0.1,
    sd_beta_MStoOE_out = 0.05,
    mean_beta_OStoOE_out = 0.1,
    sd_beta_OStoOE_out = 0.05,
    # 不满足 inside 假设的比例
    prop_negative_pleiotropy = 0.5,
    # 选型婚配基因相关
    compatibility_selection_prop = 0,
    compatibility_selection_geno = "independent",
    correlation_param = 0.5,
    # 选项混配暴露相关
    compatibility_selection_factor_exp = 0,
    compatibility_selection_factor_out = 0,
    crowd_stratification_differences = 0,
    # 因果效应
    beta_exp_to_out = 0,
    beta_confounding_exp = 0.2,
    beta_confounding_out = 0.2, correlation = 0.2, seed = NULL) {
    results <- list(
        beta_exp = NA,
        beta_se_exp = NA,
        beta_out = NA,
        beta_se_out = NA
    )
    data <- generate_multiple_datasets_v3(
        n = n, num_pleiotropic = num_pleiotropic,
        N_out = N_out,
        # 样本重叠
        overlap_prop = overlap_prop,
        # 工具变量的强度
        beta_FStoOE_exp = beta_FStoOE_exp, beta_MStoOE_exp = beta_MStoOE_exp,
        beta_OStoOE_exp = beta_OStoOE_exp,
        # 水平多效性的强度
        mean_beta_FStoOE_out = mean_beta_FStoOE_out,
        sd_beta_FStoOE_out = sd_beta_FStoOE_out,
        mean_beta_MStoOE_out = mean_beta_MStoOE_out,
        sd_beta_MStoOE_out = sd_beta_MStoOE_out,
        mean_beta_OStoOE_out = mean_beta_OStoOE_out,
        sd_beta_OStoOE_out = sd_beta_OStoOE_out,
        # 不满足 inside 假设的比例
        prop_negative_pleiotropy = prop_negative_pleiotropy,
        # 选型婚配基因相关
        compatibility_selection_prop = compatibility_selection_prop,
        compatibility_selection_geno = compatibility_selection_geno,
        correlation_param = correlation_param,
        # 选项混配暴露相关
        compatibility_selection_factor_exp = compatibility_selection_factor_exp,
        compatibility_selection_factor_out = compatibility_selection_factor_out,
        crowd_stratification_differences = crowd_stratification_differences,
        # 因果效应
        beta_exp_to_out = beta_exp_to_out,
        beta_confounding_exp = beta_confounding_exp,
        beta_confounding_out = beta_confounding_out, correlation = correlation, seed = seed
    )
    data_1_full <- data[[1]] # 完整的暴露数据
    data_2_full <- data[[2]] # 完整的结局数据 (包含 snp_type 列)

    # 进行FGWAS处理
    hat_expose <- fgwas_for_data_optimized(data_1_full)
    hat_outcome <- fgwas_for_data_optimized(data_2_full)
    hat_expose_trMR <- fgwas_to_mr(hat_expose)
    hat_outcome_trMR <- fgwas_to_mr(hat_outcome)

    results$beta_exp <- hat_expose_trMR$results_of_fgwas_beta
    results$beta_se_exp <- hat_expose_trMR$beta_se
    results$beta_out <- hat_outcome_trMR$results_of_fgwas_beta
    results$beta_se_out <- hat_outcome_trMR$beta_se
    return(results)
}

