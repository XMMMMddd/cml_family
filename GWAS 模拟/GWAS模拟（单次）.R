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
triplet_family_simulation_once_2 <- function(
    n_snps = 3, n_pleiotropy = 1,
    n_independent = 1000, p_trio = 0.5,
    p_exp_out = 0.5, p_overlap = 0,
    p_f = 0.3, p_m = 0.3, # p_m 当前未在SNP生成中使用
    # 暴露效应
    beta_fs_oe_exp = 0.3, beta_ms_oe_exp = 0.3,
    beta_os_oe_exp = 0.3,
    # 结局效应 (直接多效性 / 遗传叠加效应)
    beta_fs_oe_out = 0, beta_ms_oe_out = 0,
    beta_os_oe_out = 0, p_negative_pleiotropy = 0,
    # 因果效应
    beta_exp_to_out = 0,
    # 混杂效应
    var_confounding_exp = 0.2, var_confounding_out = 0.2,
    # 其他参数
    r_correlation = 0.2, n_seed = NULL,
    # 选型婚配强度(跨性状)
    assortative_mating_strength = 1000) {
    # 首先声明结果
    results <- list(
        a_theta_point = NA, a_theta_se = NA, a_z = NA,
        a_p_value = NA, a_duration = NA,
        b_theta_point = NA, b_theta_se = NA, b_z = NA,
        b_p_value = NA, b_duration = NA
    ) # 初始化结果列表

    # 生成数据
    phase_two_data_full <- generate_mr_trio_data_ultra_updata(
        n_snps = n_snps, n_pleiotropy = n_pleiotropy,
        n_independent = n_independent, p_trio = p_trio,
        p_exp_out = p_exp_out, p_overlap = p_overlap,
        p_f = p_f, p_m = p_m, # p_m 当前未在SNP生成中使用
        # 暴露效应
        beta_fs_oe_exp = beta_fs_oe_exp, beta_ms_oe_exp = beta_ms_oe_exp,
        beta_os_oe_exp = beta_os_oe_exp,
        # 结局效应 (直接多效性 / 遗传叠加效应)
        beta_fs_oe_out = beta_fs_oe_out, beta_ms_oe_out = beta_ms_oe_out,
        beta_os_oe_out = beta_os_oe_out,
        p_negative_pleiotropy = p_negative_pleiotropy,
        # 因果效应
        beta_exp_to_out = beta_exp_to_out,
        # 混杂效应
        var_confounding_exp = var_confounding_exp,
        var_confounding_out = var_confounding_out,
        # 其他参数
        r_correlation = r_correlation, n_seed = n_seed,
        # 选型婚配强度(跨性状)
        assortative_mating_strength = assortative_mating_strength
    )
    phase_two_data_analysis <- perform_fgwas_analysis_gls(phase_two_data_full,
        n_snps = n_snps
    )


    phase_two_data_analysis_exp <- list(
        beta_hat =
            phase_two_data_analysis$beta_hat_exp,
        sigma_inv =
            phase_two_data_analysis$beta_sigma_exp
    )
    phase_two_data_analysis_exp <- fgwas_to_mr(
        phase_two_data_analysis_exp
    )
    phase_two_data_analysis_exp$p_value <-
        p_value <- 2 * pnorm(
            abs(
                phase_two_data_analysis_out$results_of_fgwas_beta /
                    phase_two_data_analysis_out$beta_se
            ),
            lower.tail = FALSE
        )
    phase_two_data_analysis_out <- list(
        beta_hat =
            phase_two_data_analysis$beta_hat_out,
        sigma_inv =
            phase_two_data_analysis$beta_sigma_out
    )

    phase_two_data_analysis_out <- fgwas_to_mr(
        phase_two_data_analysis_out
    )
    phase_two_data_analysis_out$p_value <-
        p_value <- 2 * pnorm(
            abs(
                phase_two_data_analysis_out$results_of_fgwas_beta /
                    phase_two_data_analysis_out$beta_se
            ),
            lower.tail = FALSE
        )



    return(list(
        p_value_1 = a_1$coefficients[8], p_value_2 = a_2$coefficients[8],
        p_value_3 = a_3$coefficients[8],
        phase_two_data_analysis_exp, phase_two_data_analysis_out
    ))
}
# %%
triplet_family_simulation_once_2(
    n_snps = 3, n_pleiotropy = 0,
    n_independent = 2000, p_trio = 0.5,
    p_exp_out = 0.5, p_overlap = 0,
    p_f = 0.3, p_m = 0.3, # p_m 当前未在SNP生成中使用
    # 暴露效应
    beta_fs_oe_exp = 0, beta_ms_oe_exp = 0,
    beta_os_oe_exp = 0.9,
    # 结局效应 (直接多效性 / 遗传叠加效应)
    beta_fs_oe_out = 0, beta_ms_oe_out = 0,
    beta_os_oe_out = 0, p_negative_pleiotropy = 0,
    # 因果效应
    beta_exp_to_out = 0,
    # 混杂效应
    var_confounding_exp = 0.2, var_confounding_out = 0.2,
    # 其他参数
    r_correlation = 0.2, n_seed = NULL,
    # 选型婚配强度(跨性状)
    assortative_mating_strength = 1000
)
