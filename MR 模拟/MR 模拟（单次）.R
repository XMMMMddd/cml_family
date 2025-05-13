# 加载必要的函数
# %%
# install.packages("devtools") # 安装其他代码的包
# library(devtools)

# devtools::install_github("xue-hr/MRcML")
library(MRcML) # cML的实现

# devtools::install_github("rondolab/MR-PRESSO") # MR-PRSSO实现
library(MRPRESSO) # MR_PRESSO 实现
# install.packages("MendelianRandomization")
library(MendelianRandomization)

# devtools::install_github("gqi/MRMix")
library(MRMix)

# devtools::install_github("noahlorinczcomi/MRBEE")
library(MRBEE) # 暂时没用到
# %%
source("样本生成函数/三联体家庭结构生成函数.R") # 加载数据生成函数 library(dplyr) library(MASS) 里面加载了这些包
source("FGWAS 优化版本/FGWAS 函数.R") # 加载 FGWAS 函数1
source("cml家庭函数/cml_friamly.r") # library(RcppArmadillo) library(Rcpp)
source("cml家庭函数/cml_oracle.r")
# %%
## 计算se和估计值的函数
calculate_mr_stats <- function(estimate, se) {
    # 计算 Z 分数
    z_score <- estimate / se
    # 计算 P 值 (双侧检验)
    p_val <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
    # 返回包含 Z 分数和 P 值的列表
    return(list(z = z_score, p_value = p_val))
}

# %%
# ! 写一个数据生成的流程

triplet_family_simulation_once <- function(
    n = 10, num_pleiotropic = 0, N_out = 4000,
    # 样本重叠
    overlap_prop = 0,
    # 直接效应
    beta_FStoOE_exp = 0.3, beta_MStoOE_exp = 0.3, beta_OStoOE_exp = 0.3,
    # 水平多效性
    prop_negative_pleiotropy = 0, # 不满足inside假设
    mean_beta_FStoOE_out = 0, sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0, sd_beta_MStoOE_out = 0.05,
    mean_beta_OStoOE_out = 0, sd_beta_OStoOE_out = 0.05,
    # 选型婚配
    ## 选型婚配基因型
    compatibility_selection_prop = 0,
    compatibility_selection_geno = "independent", correlation_param = 0.5,
    ## 选型婚配暴露相关
    compatibility_selection_factor_exp = 0,
    compatibility_selection_factor_out = 0,
    # 人群分层（双人群差异）
    crowd_stratification_differences = 0,
    # 其他参数设置
    beta_exp_to_out = 0, # 因果效应
    beta_confounding_exp = 0.2,
    beta_confounding_out = 0.2,
    correlation = 0.2,
    seed = NULL) {
    # 首先声明结果
    results <- list(
        A_theta_point = NA, A_theta_se = NA, A_z = NA, A_p_value = NA, A_duration = NA,
        B_theta_point = NA, B_theta_se = NA, B_z = NA, B_p_value = NA, B_duration = NA,
        C_theta_point = NA, C_theta_se = NA, C_z = NA, C_p_value = NA, C_duration = NA,
        D_theta_point = NA, D_theta_se = NA, D_z = NA, D_p_value = NA, D_duration = NA,
        E_theta_point = NA, E_theta_se = NA, E_z = NA, E_p_value = NA, E_duration = NA,
        F_theta_point = NA, F_theta_se = NA, F_z = NA, F_p_value = NA, F_duration = NA,
        G_theta_point = NA, G_theta_se = NA, G_z = NA, G_p_value = NA, G_duration = NA,
        H_theta_point = NA, H_theta_se = NA, H_z = NA, H_p_value = NA, H_duration = NA,
        I_theta_point = NA, I_theta_se = NA, I_z = NA, I_p_value = NA, I_duration = NA
    ) # 初始化结果列表

    # 生成数据
    phase_two_data_full <- generate_multiple_datasets_v3(
        n = n, num_pleiotropic = num_pleiotropic, N_out = N_out,
        # 样本重叠
        overlap_prop = overlap_prop,
        # 直接效应
        beta_FStoOE_exp = beta_FStoOE_exp,
        beta_MStoOE_exp = beta_MStoOE_exp,
        beta_OStoOE_exp = beta_OStoOE_exp,
        # 水平多效性
        prop_negative_pleiotropy = prop_negative_pleiotropy, # 不满足inside假设
        mean_beta_FStoOE_out = mean_beta_FStoOE_out,
        sd_beta_FStoOE_out = sd_beta_FStoOE_out,
        mean_beta_MStoOE_out = mean_beta_MStoOE_out,
        sd_beta_MStoOE_out = sd_beta_MStoOE_out,
        mean_beta_OStoOE_out = mean_beta_OStoOE_out,
        sd_beta_OStoOE_out = sd_beta_OStoOE_out,
        # 选型婚配
        ## 选型婚配基因型
        compatibility_selection_prop = compatibility_selection_prop,
        compatibility_selection_geno = compatibility_selection_geno,
        correlation_param = correlation_param,
        ## 选型婚配暴露
        compatibility_selection_factor_exp = compatibility_selection_factor_exp,
        compatibility_selection_factor_out = compatibility_selection_factor_out,
        # 人群分层（双人群差异）
        crowd_stratification_differences = crowd_stratification_differences,
        # 其他参数设置
        beta_exp_to_out = beta_exp_to_out,
        beta_confounding_exp = beta_confounding_exp,
        beta_confounding_out = beta_confounding_out,
        correlation = correlation, seed = NULL
    )
    phase_one_data <- phase_two_data_full[[1]] # 完整的暴露数据
    phase_two_data <- phase_two_data_full[[2]] # 完整的结局数据 (包含 snp_type 列)
    phase_one_data_oracle <- phase_one_data %>%
        filter(snp_type != "pleiotropic_variable") # nolint
    phase_two_data_oracle <- phase_two_data %>% # nolint
        filter(snp_type != "pleiotropic_variable") # nolint
    # 进行FGWAS处理
    hat_expose <- fgwas_for_data_optimized(phase_one_data)
    hat_outcome <- fgwas_for_data_optimized(phase_two_data)
    # 只考虑好的 SNPs
    hat_expose_oracle <- fgwas_for_data_optimized(phase_one_data_oracle)
    hat_outcome_oracle <- fgwas_for_data_optimized(phase_two_data_oracle)
    # 处理成可以进行 MR 的样子
    hat_expose_trMR <- fgwas_to_mr(hat_expose)
    hat_outcome_trMR <- fgwas_to_mr(hat_outcome)

    # A: cml_overlap
    res_A_time_1 <- Sys.time()
    res_A <- mr_cML_Overlap(hat_expose_trMR$results_of_fgwas_beta,
        hat_outcome_trMR$results_of_fgwas_beta,
        hat_expose_trMR$beta_se,
        hat_outcome_trMR$beta_se,
        n = N_out, # 确认样本量参数
        rho = 0
    )
    res_A_time_2 <- Sys.time()
    res_A_estimtor <- calculate_mr_stats(res_A$MA_BIC_theta, res_A$MA_BIC_se)
    results$A_theta_point <- res_A$MA_BIC_theta
    results$A_theta_se <- res_A$MA_BIC_se
    results$A_z <- res_A_estimtor$z
    results$A_p_value <- res_A$MA_BIC_p
    results$A_duration <- as.numeric(res_A_time_2 - res_A_time_1)


    # B: cml_family 方法
    res_B_time_1 <- Sys.time()
    res_B <- cml_family_2_cn_cpp(
        beta_y_hat = hat_outcome$beta_hat,
        beta_x_hat = hat_expose$beta_hat,
        Sigma_inv_x = hat_expose$Sigma_inv,
        Sigma_inv_y = hat_outcome$Sigma_inv
    )
    res_B_time_2 <- Sys.time()
    results$B_theta_point <- res_B$theta_weight
    results$B_theta_se <- res_B$theta_se_weight
    res_B_estimtor <- calculate_mr_stats(res_B$theta_weight, res_B$theta_se_weight)
    results$B_z <- res_B_estimtor$z
    results$B_p_value <- res_B_estimtor$p_value
    results$B_duration <- as.numeric(res_B_time_2 - res_B_time_1)


    # C:  mr_presso 方法
    presso_method_input <- data.frame(
        beta_exp = hat_expose_trMR$results_of_fgwas_beta,
        se_exp = hat_expose_trMR$beta_se,
        beta_out = hat_outcome_trMR$results_of_fgwas_beta,
        se_out = hat_outcome_trMR$beta_se
    )
    res_C_time_1 <- Sys.time()
    res_C <- mr_presso(
        BetaOutcome = "beta_out",
        BetaExposure = "beta_exp", SdOutcome = "se_out",
        SdExposure = "se_exp",
        OUTLIERtest = TRUE,
        data = presso_method_input,
        DISTORTIONtest = TRUE,
        NbDistribution = 1000,
        SignifThreshold = 0.05
    )
    res_C_time_2 <- Sys.time()

    results$C_theta_point <- res_C$`Main MR results`[1, 3]
    results$C_theta_se <- res_C$`Main MR results`[1, 4]
    results$C_z <- res_C$`Main MR results`[1, 5]
    results$C_p_value <- res_C$`Main MR results`[1, 6]
    results$C_duration <- as.numeric(res_C_time_2 - res_C_time_1)

    # D: cml_oracle
    res_D_time_1 <- Sys.time()
    res_D <- orcacle_function(
        beta_y_hat = hat_outcome_oracle$beta_hat,
        beta_x_hat = hat_expose_oracle$beta_hat,
        Sigma_inv_x = hat_expose_oracle$Sigma_inv,
        Sigma_inv_y = hat_outcome_oracle$Sigma_inv
    )
    res_D_time_2 <- Sys.time()
    results$D_theta_point <- res_D[1]
    results$D_theta_se <- res_D[2]
    res_D_estimtor <- calculate_mr_stats(res_D[1], res_D[2])
    results$D_z <- res_D_estimtor$z
    results$D_p_value <- res_D_estimtor$p_value
    results$D_duration <- as.numeric(res_D_time_2 - res_D_time_1)

    # E: IVW
    ## 做成能用MendelianRandomization这个包导入的数据
    mr_input_obj <- MendelianRandomization::mr_input(
        bx = hat_expose_trMR$results_of_fgwas_beta,
        bxse = hat_expose_trMR$beta_se,
        by = hat_outcome_trMR$results_of_fgwas_beta,
        byse = hat_outcome_trMR$beta_se
    )
    res_E_time_1 <- Sys.time()
    res_E <- MendelianRandomization::mr_ivw(mr_input_obj)
    res_E_time_2 <- Sys.time()
    results$E_theta_point <- res_E$Estimate
    results$E_theta_se <- res_E$StdError
    res_E_estimtor <- calculate_mr_stats(res_E$Estimate, res_E$StdError)
    results$E_z <- res_E_estimtor$z
    results$E_p_value <- res_E$Pvalue
    results$E_duration <- as.numeric(res_E_time_2 - res_E_time_1)

    # F: MR-Egger
    res_F_time_1 <- Sys.time()
    res_F <- MendelianRandomization::mr_egger(mr_input_obj)
    res_F_time_2 <- Sys.time()
    results$F_theta_point <- res_F$Estimate
    results$F_theta_se <- res_F$StdError.Est
    res_F_estimtor <- calculate_mr_stats(res_F$Estimate, res_F$StdError.Est)
    results$F_z <- res_F_estimtor$z
    results$F_p_value <- res_F$Causal.pval
    results$F_duration <- as.numeric(res_F_time_2 - res_F_time_1)


    # G: MRMix
    res_G_time_1 <- Sys.time()
    res_G <- MRMix(
        betahat_x = hat_expose_trMR$results_of_fgwas_beta,
        betahat_y = hat_outcome_trMR$results_of_fgwas_beta,
        sx = hat_expose_trMR$beta_se,
        sy = hat_outcome_trMR$beta_se
    )
    res_G_time_2 <- Sys.time()
    results$G_theta_point <- res_G$theta
    results$G_theta_se <- res_G$SE_theta
    results$G_z <- res_G$zstat_theta
    results$G_p_value <- res_G$pvalue_theta
    results$G_duration <- as.numeric(res_G_time_2 - res_G_time_1)

    # H: MRBEE
    # debug(MRBEE.IMRP) # 开启调试模式
    # fit <- MRBEE.IMRP(
    #     by = matrix(hat_outcome_trMR$results_of_fgwas_beta),
    #     bX = matrix(hat_expose_trMR$results_of_fgwas_beta),
    #     byse = matrix(hat_outcome_trMR$beta_se),
    #     bXse = matrix(hat_expose_trMR$beta_se),
    #     Rxy = diag(1, 2)
    # )


    # I：MR-ConMix
    res_I_time_1 <- Sys.time()
    res_I <- MendelianRandomization::mr_conmix(mr_input_obj)
    res_I_time_2 <- Sys.time()
    results$I_theta_point <- res_I$Estimate
    # results$I_theta_se <- res_I$StdError.Est
    # res_I_estimtor <- calculate_mr_stats(res_I$Estimate, res_I$StdError.Est)
    # results$I_z <- res_I_estimtor$z
    results$I_p_value <- res_I$Pvalue
    results$I_duration <- as.numeric(res_I_time_2 - res_I_time_1)
    return(results)
}

