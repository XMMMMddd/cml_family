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
source("cml家庭函数/cml家庭函数一代/cml_friamly.r") # library(RcppArmadillo) library(Rcpp)
source("cml家庭函数/cml家庭函数一代/cml_oracle.r")
source("cml家庭函数/cml家庭函数二代/cml_family_version2.r")
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
    assortative_mating_prob = 0, # 选型婚配比例
    assortative_mating_strength = 1000, # 选型婚配强度
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
        assortative_mating_prob = assortative_mating_prob, # 选型婚配比例
        assortative_mating_strength = assortative_mating_strength, # 选型婚配强度
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


    # C:  cml_family_ver2 方法
    res_C_time_1 <- Sys.time()
    res_C <- cml_family_ver_2_cpp(
        beta_hat_out = hat_outcome$beta_hat,
        beta_hat_exp = hat_expose$beta_hat,
        beta_sigma_exp = hat_expose$Sigma_inv,
        beta_sigma_out = hat_outcome$Sigma_inv
    )
    res_C_time_2 <- Sys.time()
    results$C_theta_point <- res_C[[1]]
    results$C_theta_se <- res_C[[2]]
    res_C_estimtor <- calculate_mr_stats(results$C_theta_point, results$C_theta_se)
    results$C_z <- res_C_estimtor$z
    results$C_p_value <- res_C_estimtor$p_value
    results$C_duration <- as.numeric(res_C_time_2 - res_C_time_1)


    # D: cml_oracle
    res_D_time_1 <- Sys.time()
    res_D <- oracle_function(
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
    print("前4个成功")
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
    results$I_p_value <- res_I$Pvalue
    results$I_duration <- as.numeric(res_I_time_2 - res_I_time_1)


    return(results)
}

triplet_family_simulation_once_robust_mr <- function(
    n = 10, num_pleiotropic = 0, N_out = 4000,
    # 样本重叠
    overlap_prop = 0,
    # 直接效应
    beta_FStoOE_exp = 0.3, beta_MStoOE_exp = 0.3, beta_OStoOE_exp = 0.3,
    # 水平多效性
    prop_negative_pleiotropy = 0,
    mean_beta_FStoOE_out = 0, sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0, sd_beta_MStoOE_out = 0.05,
    mean_beta_OStoOE_out = 0, sd_beta_OStoOE_out = 0.05,
    # 选型婚配
    assortative_mating_prob = 0,
    assortative_mating_strength = 1000,
    # 人群分层（双人群差异）
    crowd_stratification_differences = 0,
    # 其他参数设置
    beta_exp_to_out = 0,
    beta_confounding_exp = 0.2,
    beta_confounding_out = 0.2,
    correlation = 0.2,
    seed = NULL) {

    # --- 初始化结果列表 (保持不变) ---
    results <- list(
        A_theta_point = NA, A_theta_se = NA, A_z = NA, A_p_value = NA, A_duration = NA,
        B_theta_point = NA, B_theta_se = NA, B_z = NA, B_p_value = NA, B_duration = NA,
        C_theta_point = NA, C_theta_se = NA, C_z = NA, C_p_value = NA, C_duration = NA,
        D_theta_point = NA, D_theta_se = NA, D_z = NA, D_p_value = NA, D_duration = NA,
        E_theta_point = NA, E_theta_se = NA, E_z = NA, E_p_value = NA, E_duration = NA,
        F_theta_point = NA, F_theta_se = NA, F_z = NA, F_p_value = NA, F_duration = NA,
        G_theta_point = NA, G_theta_se = NA, G_z = NA, G_p_value = NA, G_duration = NA,
        H_theta_point = NA, H_theta_se = NA, H_z = NA, H_p_value = NA, H_duration = NA, # 即使 H 未实现，也保留位置
        I_theta_point = NA, I_theta_se = NA, I_z = NA, I_p_value = NA, I_duration = NA,
        # 增加记录错误的字段
        A_error = NA, B_error = NA, C_error = NA, D_error = NA,
        E_error = NA, F_error = NA, G_error = NA, H_error = NA, I_error = NA
    )

    # --- 数据生成和预处理 (假设这部分没有 'wrong sign in by' 错误) ---
    # 但如果这里出错，整个函数会提前停止，除非也加上 tryCatch
    # 为了聚焦 MR 方法，暂时不在这里加 tryCatch
    phase_two_data_full <- generate_multiple_datasets_v3(
        n = n, num_pleiotropic = num_pleiotropic, N_out = N_out,
        overlap_prop = overlap_prop, beta_FStoOE_exp = beta_FStoOE_exp,
        beta_MStoOE_exp = beta_MStoOE_exp, beta_OStoOE_exp = beta_OStoOE_exp,
        prop_negative_pleiotropy = prop_negative_pleiotropy,
        mean_beta_FStoOE_out = mean_beta_FStoOE_out, sd_beta_FStoOE_out = sd_beta_FStoOE_out,
        mean_beta_MStoOE_out = mean_beta_MStoOE_out, sd_beta_MStoOE_out = sd_beta_MStoOE_out,
        mean_beta_OStoOE_out = mean_beta_OStoOE_out, sd_beta_OStoOE_out = sd_beta_OStoOE_out,
        assortative_mating_prob = assortative_mating_prob,
        assortative_mating_strength = assortative_mating_strength,
        crowd_stratification_differences = crowd_stratification_differences,
        beta_exp_to_out = beta_exp_to_out, beta_confounding_exp = beta_confounding_exp,
        beta_confounding_out = beta_confounding_out, correlation = correlation, seed = seed # 使用传入的 seed
    )
    phase_one_data <- phase_two_data_full[[1]]
    phase_two_data <- phase_two_data_full[[2]]
    phase_one_data_oracle <- phase_one_data %>%
        filter(snp_type != "pleiotropic_variable")
    phase_two_data_oracle <- phase_two_data %>%
        filter(snp_type != "pleiotropic_variable")

    hat_expose <- fgwas_for_data_optimized(phase_one_data)
    hat_outcome <- fgwas_for_data_optimized(phase_two_data)
    hat_expose_oracle <- fgwas_for_data_optimized(phase_one_data_oracle)
    hat_outcome_oracle <- fgwas_for_data_optimized(phase_two_data_oracle)

    hat_expose_trMR <- fgwas_to_mr(hat_expose)
    hat_outcome_trMR <- fgwas_to_mr(hat_outcome)

    # --- MR 方法计算 (增加 tryCatch) ---

    # A: cml_overlap
    res_A_time_1 <- Sys.time()
    tryCatch({
        res_A <- mr_cML_Overlap(hat_expose_trMR$results_of_fgwas_beta,
            hat_outcome_trMR$results_of_fgwas_beta,
            hat_expose_trMR$beta_se,
            hat_outcome_trMR$beta_se,
            n = N_out, rho = 0)
        res_A_estimtor <- calculate_mr_stats(res_A$MA_BIC_theta, res_A$MA_BIC_se)
        results$A_theta_point <- res_A$MA_BIC_theta
        results$A_theta_se <- res_A$MA_BIC_se
        results$A_z <- res_A_estimtor$z
        results$A_p_value <- res_A$MA_BIC_p
    }, error = function(e) {
        results$A_error <<- e$message # 使用 <<- 赋值给父环境的 results
        warning("Error in MR method A (cml_overlap): ", e$message, call. = FALSE)
    })
    res_A_time_2 <- Sys.time()
    results$A_duration <- as.numeric(difftime(res_A_time_2, res_A_time_1, units = "secs"))


    # B: cml_family 方法
    res_B_time_1 <- Sys.time()
    tryCatch({
        res_B <- cml_family_2_cn_cpp(
            beta_y_hat = hat_outcome$beta_hat, beta_x_hat = hat_expose$beta_hat,
            Sigma_inv_x = hat_expose$Sigma_inv, Sigma_inv_y = hat_outcome$Sigma_inv)
        res_B_estimtor <- calculate_mr_stats(res_B$theta_weight, res_B$theta_se_weight)
        results$B_theta_point <- res_B$theta_weight
        results$B_theta_se <- res_B$theta_se_weight
        results$B_z <- res_B_estimtor$z
        results$B_p_value <- res_B_estimtor$p_value
    }, error = function(e) {
        results$B_error <<- e$message
        warning("Error in MR method B (cml_family): ", e$message, call. = FALSE)
    })
    res_B_time_2 <- Sys.time()
    results$B_duration <- as.numeric(difftime(res_B_time_2, res_B_time_1, units = "secs"))


    # C: cml_family_ver2 方法
    res_C_time_1 <- Sys.time()
    tryCatch({
        res_C <- cml_family_ver_2_cpp( # 假设这是你之前有耗时分析的那个函数
            beta_hat_out = hat_outcome$beta_hat, beta_hat_exp = hat_expose$beta_hat,
            beta_sigma_exp = hat_expose$Sigma_inv, beta_sigma_out = hat_outcome$Sigma_inv)
        # 假设 res_C 返回一个列表，第一个元素是 theta, 第二个是 se
        results$C_theta_point <- res_C[[1]]
        results$C_theta_se <- res_C[[2]]
        res_C_estimtor <- calculate_mr_stats(results$C_theta_point, results$C_theta_se)
        results$C_z <- res_C_estimtor$z
        results$C_p_value <- res_C_estimtor$p_value
    }, error = function(e) {
        results$C_error <<- e$message
        warning("Error in MR method C (cml_family_ver2): ", e$message, call. = FALSE)
    })
    res_C_time_2 <- Sys.time()
    results$C_duration <- as.numeric(difftime(res_C_time_2, res_C_time_1, units = "secs"))


    # D: cml_oracle
    res_D_time_1 <- Sys.time()
    tryCatch({
        res_D <- oracle_function( # 拼写修正：oracle_function
            beta_y_hat = hat_outcome_oracle$beta_hat, beta_x_hat = hat_expose_oracle$beta_hat,
            Sigma_inv_x = hat_expose_oracle$Sigma_inv, Sigma_inv_y = hat_outcome_oracle$Sigma_inv)
        # 假设 res_D 返回一个包含 theta 和 se 的向量
        results$D_theta_point <- res_D[1]
        results$D_theta_se <- res_D[2]
        res_D_estimtor <- calculate_mr_stats(res_D[1], res_D[2])
        results$D_z <- res_D_estimtor$z
        results$D_p_value <- res_D_estimtor$p_value
    }, error = function(e) {
        results$D_error <<- e$message
        warning("Error in MR method D (cml_oracle): ", e$message, call. = FALSE)
    })
    res_D_time_2 <- Sys.time()
    results$D_duration <- as.numeric(difftime(res_D_time_2, res_D_time_1, units = "secs"))

    # 确保输入给 MendelianRandomization 的数据是有效的
    # 检查向量长度是否至少为1，并且不全是NA
    can_run_mr_package_methods <- FALSE
    if (length(hat_expose_trMR$results_of_fgwas_beta) > 0 &&
        length(hat_expose_trMR$beta_se) > 0 &&
        length(hat_outcome_trMR$results_of_fgwas_beta) > 0 &&
        length(hat_outcome_trMR$beta_se) > 0 &&
        !all(is.na(hat_expose_trMR$results_of_fgwas_beta)) &&
        !all(is.na(hat_expose_trMR$beta_se)) &&
        !all(is.na(hat_outcome_trMR$results_of_fgwas_beta)) &&
        !all(is.na(hat_outcome_trMR$beta_se)) &&
        length(hat_expose_trMR$results_of_fgwas_beta) == length(hat_expose_trMR$beta_se) &&
        length(hat_expose_trMR$results_of_fgwas_beta) == length(hat_outcome_trMR$results_of_fgwas_beta) &&
        length(hat_expose_trMR$results_of_fgwas_beta) == length(hat_outcome_trMR$beta_se) ) {
            can_run_mr_package_methods <- TRUE
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = hat_expose_trMR$results_of_fgwas_beta,
                bxse = hat_expose_trMR$beta_se,
                by = hat_outcome_trMR$results_of_fgwas_beta,
                byse = hat_outcome_trMR$beta_se
            )
    } else {
        warning_msg_mr_pkg <- "Skipping MendelianRandomization package methods (E, F, I) due to invalid input data (empty, all NA, or mismatched lengths)."
        warning(warning_msg_mr_pkg, call. = FALSE)
        results$E_error <<- warning_msg_mr_pkg
        results$F_error <<- warning_msg_mr_pkg
        results$I_error <<- warning_msg_mr_pkg
        # G 方法 (MRMix) 也有类似的数据需求
        results$G_error <<- "Skipping MRMix due to invalid input data (same reason as above)."
    }


    # E: IVW
    res_E_time_1 <- Sys.time()
    if(can_run_mr_package_methods){
        tryCatch({
            res_E <- MendelianRandomization::mr_ivw(mr_input_obj)
            results$E_theta_point <- res_E$Estimate
            results$E_theta_se <- res_E$StdError
            res_E_estimtor <- calculate_mr_stats(res_E$Estimate, res_E$StdError)
            results$E_z <- res_E_estimtor$z
            results$E_p_value <- res_E$Pvalue
        }, error = function(e) {
            results$E_error <<- e$message
            warning("Error in MR method E (IVW): ", e$message, call. = FALSE)
        })
    }
    res_E_time_2 <- Sys.time()
    results$E_duration <- as.numeric(difftime(res_E_time_2, res_E_time_1, units = "secs"))


    # F: MR-Egger
    res_F_time_1 <- Sys.time()
    if(can_run_mr_package_methods){
        tryCatch({
            res_F <- MendelianRandomization::mr_egger(mr_input_obj)
            results$F_theta_point <- res_F$Estimate
            results$F_theta_se <- res_F$StdError.Est # 注意这里的列名可能与IVW不同
            res_F_estimtor <- calculate_mr_stats(res_F$Estimate, res_F$StdError.Est)
            results$F_z <- res_F_estimtor$z
            results$F_p_value <- res_F$Pvalue.Est # 注意这里的列名 Pvalue.Est 或 Causal.pval
        }, error = function(e) {
            results$F_error <<- e$message
            warning("Error in MR method F (MR-Egger): ", e$message, call. = FALSE)
        })
    }
    res_F_time_2 <- Sys.time()
    results$F_duration <- as.numeric(difftime(res_F_time_2, res_F_time_1, units = "secs"))


    # G: MRMix
    res_G_time_1 <- Sys.time()
    if(can_run_mr_package_methods){ # MRMix 也需要有效的数据输入
        tryCatch({
            # 确保 MRMix 的输入不包含 NA/Inf，否则它可能会出错
            # 可以在这里添加更严格的检查，或者依赖 MRMix 自身的错误处理
            # 确保向量长度至少为1
             if(length(hat_expose_trMR$results_of_fgwas_beta) < 1) stop("MRMix requires at least one SNP.")

            res_G <- MRMix::MRMix( # 假设 MRMix 包已加载，或者用 MRMix::MRMix
                betahat_x = hat_expose_trMR$results_of_fgwas_beta,
                betahat_y = hat_outcome_trMR$results_of_fgwas_beta,
                sx = hat_expose_trMR$beta_se,
                sy = hat_outcome_trMR$beta_se
            )
            results$G_theta_point <- res_G$theta
            results$G_theta_se <- res_G$SE_theta
            # MRMix 直接提供 z 和 p，不需要 calculate_mr_stats
            results$G_z <- res_G$zstat_theta
            results$G_p_value <- res_G$pvalue_theta
        }, error = function(e) {
            results$G_error <<- e$message
            warning("Error in MR method G (MRMix): ", e$message, call. = FALSE)
        })
    }
    res_G_time_2 <- Sys.time()
    results$G_duration <- as.numeric(difftime(res_G_time_2, res_G_time_1, units = "secs"))


    # H: MRBEE (当前被注释掉了，如果启用，也需要 tryCatch)
    # results$H_error <<- "MRBEE not run (commented out)."


    # I：MR-ConMix
    res_I_time_1 <- Sys.time()
    if(can_run_mr_package_methods){
        tryCatch({
            res_I <- MendelianRandomization::mr_conmix(mr_input_obj)
            results$I_theta_point <- res_I$Estimate
            results$I_p_value <- res_I$Pvalue
            # ConMix 可能不直接提供 SE，如果需要，可能要从其他地方获取或计算
            # results$I_theta_se <- res_I$StdError # 检查 res_I 是否有这个字段
            # res_I_estimtor <- calculate_mr_stats(res_I$Estimate, results$I_theta_se)
            # results$I_z <- res_I_estimtor$z
            # 如果没有 SE 和 Z，可以将它们留为 NA
            results$I_theta_se <- ifelse("StdError" %in% names(res_I), res_I$StdError, NA)
            if (!is.na(results$I_theta_se)) {
                res_I_estimtor <- calculate_mr_stats(results$I_theta_point, results$I_theta_se)
                results$I_z <- res_I_estimtor$z
            } else {
                results$I_z <- NA
            }

        }, error = function(e) {
            results$I_error <<- e$message
            warning("Error in MR method I (MR-ConMix): ", e$message, call. = FALSE)
        })
    }
    res_I_time_2 <- Sys.time()
    results$I_duration <- as.numeric(difftime(res_I_time_2, res_I_time_1, units = "secs"))

    return(results)
}



