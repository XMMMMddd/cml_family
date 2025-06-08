# %% 加载必要的函数

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
# %% 加载必要的包
source("样本生成函数/三联体家庭结构生成函数.R") # 加载数据生成函数 library(dplyr) library(MASS) 里面加载了这些包
source("FGWAS 优化版本/FGWAS 函数.R") # 加载 FGWAS 函数1
source("cml家庭函数/cml家庭函数一代/cml_oracle.r")
source("cml家庭函数/cml家庭函数overlap版本/cml_family_overlap2.r")
# %% 辅助函数

## 计算se和估计值的函数
calculate_mr_stats <- function(estimate, se) {
    # 计算 Z 分数
    z_score <- estimate / se
    # 计算 P 值 (双侧检验)
    p_val <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
    # 返回包含 Z 分数和 P 值的列表
    return(list(z = z_score, p_value = p_val))
}
# %% 辅助代码
n_snps <- 3
n_pleiotropy <- 1
n_independent <- 1000
p_trio <- 0.5
p_exp_out <- 0.5
p_overlap <- 0
p_f <- 0.3
p_m <- 0.3 # p_m 当前未在SNP生成中使用
# 暴露效应
beta_fs_oe_exp <- 0.3
beta_ms_oe_exp <- 0
beta_os_oe_exp <- 0
# 结局效应 (直接多效性 / 遗传叠加效应)
beta_fs_oe_out <- 0
beta_ms_oe_out <- 0
beta_os_oe_out <- 0
p_negative_pleiotropy <- 0
# 因果效应
beta_exp_to_out <- 0.4
# 混杂效应
var_confounding_exp <- 0.2
var_confounding_out <- 0.2
# 其他参数
r_correlation <- 0.2
n_seed <- NULL
# 选型婚配强度(跨性状)
assortative_mating_strength <- 1000
# %% #!写一个数据生成的流程


triplet_family_simulation_once <- function(
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
    beta_exp_to_out = 0.4,
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
        b_p_value = NA, b_duration = NA,
        c_theta_point = NA, c_theta_se = NA, c_z = NA,
        c_p_value = NA, c_duration = NA
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
    phase_two_data_analysis <- perform_fgwas_analysis(phase_two_data_full,
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
    phase_two_data_analysis_out <- list(
        beta_hat =
            phase_two_data_analysis$beta_hat_out,
        sigma_inv =
            phase_two_data_analysis$beta_sigma_out
    )
    phase_two_data_analysis_out <- fgwas_to_mr(
        phase_two_data_analysis_out
    )
    # a: cml_overlap
    res_a_time_1 <- Sys.time()
    res_a <- mr_cML_Overlap(
        phase_two_data_analysis_exp$results_of_fgwas_beta,
        phase_two_data_analysis_out$results_of_fgwas_beta,
        phase_two_data_analysis_exp$beta_se,
        phase_two_data_analysis_out$beta_se,
        n = floor(n_independent * p_trio), # 确认样本量参数
        rho = 0
    )
    res_a_time_2 <- Sys.time()
    res_a_estimtor <- calculate_mr_stats(res_a$MA_BIC_theta, res_a$MA_BIC_se)
    results$a_theta_point <- res_a$MA_BIC_theta
    results$a_theta_se <- res_a$MA_BIC_se
    results$a_z <- res_a_estimtor$z
    results$a_p_value <- res_a$MA_BIC_p
    results$a_duration <- as.numeric(res_a_time_2 - res_a_time_1)

    # b:  cml_family_ver2 方法
    beta_sigma_exp <- as.matrix(phase_two_data_analysis$beta_sigma_exp)
    beta_sigma_out <- as.matrix(phase_two_data_analysis $beta_sigma_out)
    beta_hat_exp <- as.matrix(phase_two_data_analysis $beta_hat_exp)
    beta_hat_out <- as.matrix(phase_two_data_analysis $beta_hat_out)
    beta_sigma_t_exp <- invert_beta_sigma_out_matrices(beta_sigma_exp)
    beta_sigma_t_out <- invert_beta_sigma_out_matrices(beta_sigma_out)
    matrix_big <- create_combined_diagonal_matrices(
        beta_sigma_t_exp,
        beta_sigma_t_out,
        off_diagonal_elements = NULL
    )
    res_b_time_1 <- Sys.time()
    res_b <- cml_family_overlap_2(
        n_snps, beta_hat_exp, beta_hat_out,
        matrix_big
    )
    res_b
    res_b_time_2 <- Sys.time()
    results$b_theta_point <- res_b[[1]]
    results$b_theta_se <- res_b[[2]]
    res_b_estimtor <- calculate_mr_stats(
        results$b_theta_point,
        results$b_theta_se
    )
    results$b_z <- res_b_estimtor$z
    results$b_p_value <- res_b_estimtor$p_value
    results$b_duration <- as.numeric(res_b_time_2 - res_b_time_1)


    # c: IVW
    ## 做成能用MendelianRandomization这个包导入的数据
    mr_input_obj <- MendelianRandomization::mr_input(
        bx = phase_two_data_analysis_exp$results_of_fgwas_beta,
        bxse = phase_two_data_analysis_exp$beta_se,
        by = phase_two_data_analysis_out$results_of_fgwas_beta,
        byse = phase_two_data_analysis_out$beta_se
    )
    res_c_time_1 <- Sys.time()
    res_c <- MendelianRandomization::mr_ivw(mr_input_obj)
    res_c_time_2 <- Sys.time()
    results$c_theta_point <- res_c$Estimate
    results$c_theta_se <- res_c$StdError
    res_c_estimtor <- calculate_mr_stats(res_c$Estimate, res_c$StdError)
    results$c_z <- res_c_estimtor$z
    results$c_p_value <- res_c$Pvalue
    results$c_duration <- as.numeric(res_c_time_2 - res_c_time_1)

    # # f: MR-Egger
    # res_f_time_1 <- Sys.time()
    # res_f <- MendelianRandomization::mr_egger(mr_input_obj)
    # res_f_time_2 <- Sys.time()
    # results$f_theta_point <- res_f$Estimate
    # results$f_theta_se <- res_f$StdError.Est
    # res_f_estimtor <- calculate_mr_stats(res_f$Estimate, res_f$StdError.Est)
    # results$f_z <- res_f_estimtor$z
    # results$f_p_value <- res_f$Causal.pval
    # results$f_duration <- as.numeric(res_f_time_2 - res_f_time_1)


    # H: MRBEE
    # debug(MRBEE.IMRP) # 开启调试模式
    # fit <- MRBEE.IMRP(
    #     by = matrix(hat_outcome_trMR$results_of_fgwas_beta),
    #     bX = matrix(hat_expose_trMR$results_of_fgwas_beta),
    #     byse = matrix(hat_outcome_trMR$beta_se),
    #     bXse = matrix(hat_expose_trMR$beta_se),
    #     Rxy = diag(1, 2)
    # )




    return(results)
}
triplet_family_simulation_once(n_independent = 3000)
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
        a_theta_point = NA, a_theta_se = NA, a_z = NA, a_p_value = NA, a_duration = NA,
        b_theta_point = NA, b_theta_se = NA, b_z = NA, b_p_value = NA, b_duration = NA,
        c_theta_point = NA, c_theta_se = NA, c_z = NA, c_p_value = NA, c_duration = NA,
        d_theta_point = NA, d_theta_se = NA, d_z = NA, d_p_value = NA, d_duration = NA,
        e_theta_point = NA, e_theta_se = NA, e_z = NA, e_p_value = NA, e_duration = NA,
        f_theta_point = NA, f_theta_se = NA, f_z = NA, f_p_value = NA, f_duration = NA,
        g_theta_point = NA, g_theta_se = NA, g_z = NA, g_p_value = NA, g_duration = NA,
        h_theta_point = NA, h_theta_se = NA, h_z = NA, h_p_value = NA, h_duration = NA, # 即使 H 未实现，也保留位置
        i_theta_point = NA, i_theta_se = NA, i_z = NA, i_p_value = NA, i_duration = NA,
        # 增加记录错误的字段
        a_error = NA, b_error = NA, c_error = NA, d_error = NA,
        e_error = NA, f_error = NA, g_error = NA, h_error = NA, i_error = NA
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

    # a: cml_overlap
    res_a_time_1 <- Sys.time()
    tryCatch(
        {
            res_a <- mr_cML_Overlap(hat_expose_trMR$results_of_fgwas_beta,
                hat_outcome_trMR$results_of_fgwas_beta,
                hat_expose_trMR$beta_se,
                hat_outcome_trMR$beta_se,
                n = N_out, rho = 0
            )
            res_a_estimtor <- calculate_mr_stats(res_a$MA_BIC_theta, res_a$MA_BIC_se)
            results$a_theta_point <- res_a$MA_BIC_theta
            results$a_theta_se <- res_a$MA_BIC_se
            results$a_z <- res_a_estimtor$z
            results$a_p_value <- res_a$MA_BIC_p
        },
        error = function(e) {
            results$a_error <<- e$message # 使用 <<- 赋值给父环境的 results
            warning("Error in MR method A (cml_overlap): ", e$message, call. = FALSE)
        }
    )
    res_a_time_2 <- Sys.time()
    results$a_duration <- as.numeric(difftime(res_a_time_2, res_a_time_1, units = "secs"))


    # b: cml_family 方法
    res_b_time_1 <- Sys.time()
    tryCatch(
        {
            res_b <- cml_family_2_cn_cpp(
                beta_y_hat = hat_outcome$beta_hat, beta_x_hat = hat_expose$beta_hat,
                Sigma_inv_x = hat_expose$Sigma_inv, Sigma_inv_y = hat_outcome$Sigma_inv
            )
            res_b_estimtor <- calculate_mr_stats(res_b$theta_weight, res_b$theta_se_weight)
            results$b_theta_point <- res_b$theta_weight
            results$b_theta_se <- res_b$theta_se_weight
            results$b_z <- res_b_estimtor$z
            results$b_p_value <- res_b_estimtor$p_value
        },
        error = function(e) {
            results$b_error <<- e$message
            warning("Error in MR method B (cml_family): ", e$message, call. = FALSE)
        }
    )
    res_b_time_2 <- Sys.time()
    results$b_duration <- as.numeric(difftime(res_b_time_2, res_b_time_1, units = "secs"))


    # c: cml_family_ver2 方法
    res_c_time_1 <- Sys.time()
    tryCatch(
        {
            res_c <- cml_family_ver_2_cpp( # 假设这是你之前有耗时分析的那个函数
                beta_hat_out = hat_outcome$beta_hat, beta_hat_exp = hat_expose$beta_hat,
                beta_sigma_exp = hat_expose$Sigma_inv, beta_sigma_out = hat_outcome$Sigma_inv
            )
            # 假设 res_c 返回一个列表，第一个元素是 theta, 第二个是 se
            results$c_theta_point <- res_c[[1]]
            results$c_theta_se <- res_c[[2]]
            res_c_estimtor <- calculate_mr_stats(results$c_theta_point, results$c_theta_se)
            results$c_z <- res_c_estimtor$z
            results$c_p_value <- res_c_estimtor$p_value
        },
        error = function(e) {
            results$c_error <<- e$message
            warning("Error in MR method C (cml_family_ver2): ", e$message, call. = FALSE)
        }
    )
    res_c_time_2 <- Sys.time()
    results$c_duration <- as.numeric(difftime(res_c_time_2, res_c_time_1, units = "secs"))


    # d: cml_oracle
    res_d_time_1 <- Sys.time()
    tryCatch(
        {
            res_d <- oracle_function( # 拼写修正：oracle_function
                beta_y_hat = hat_outcome_oracle$beta_hat, beta_x_hat = hat_expose_oracle$beta_hat,
                Sigma_inv_x = hat_expose_oracle$Sigma_inv, Sigma_inv_y = hat_outcome_oracle$Sigma_inv
            )
            # 假设 res_d 返回一个包含 theta 和 se 的向量
            results$d_theta_point <- res_d[1]
            results$d_theta_se <- res_d[2]
            res_d_estimtor <- calculate_mr_stats(res_d[1], res_d[2])
            results$d_z <- res_d_estimtor$z
            results$d_p_value <- res_d_estimtor$p_value
        },
        error = function(e) {
            results$d_error <<- e$message
            warning("Error in MR method D (cml_oracle): ", e$message, call. = FALSE)
        }
    )
    res_d_time_2 <- Sys.time()
    results$d_duration <- as.numeric(difftime(res_d_time_2, res_d_time_1, units = "secs"))

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
        length(hat_expose_trMR$results_of_fgwas_beta) == length(hat_outcome_trMR$beta_se)) {
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
        results$e_error <<- warning_msg_mr_pkg
        results$f_error <<- warning_msg_mr_pkg
        results$i_error <<- warning_msg_mr_pkg
        # G 方法 (MRMix) 也有类似的数据需求
        results$g_error <<- "Skipping MRMix due to invalid input data (same reason as above)."
    }


    # e: IVW
    res_e_time_1 <- Sys.time()
    if (can_run_mr_package_methods) {
        tryCatch(
            {
                res_e <- MendelianRandomization::mr_ivw(mr_input_obj)
                results$e_theta_point <- res_e$Estimate
                results$e_theta_se <- res_e$StdError
                res_e_estimtor <- calculate_mr_stats(res_e$Estimate, res_e$StdError)
                results$e_z <- res_e_estimtor$z
                results$e_p_value <- res_e$Pvalue
            },
            error = function(e) {
                results$e_error <<- e$message
                warning("Error in MR method E (IVW): ", e$message, call. = FALSE)
            }
        )
    }
    res_e_time_2 <- Sys.time()
    results$e_duration <- as.numeric(difftime(res_e_time_2, res_e_time_1, units = "secs"))


    # f: MR-Egger
    res_f_time_1 <- Sys.time()
    if (can_run_mr_package_methods) {
        tryCatch(
            {
                res_f <- MendelianRandomization::mr_egger(mr_input_obj)
                results$f_theta_point <- res_f$Estimate
                results$f_theta_se <- res_f$StdError.Est # 注意这里的列名可能与IVW不同
                res_f_estimtor <- calculate_mr_stats(res_f$Estimate, res_f$StdError.Est)
                results$f_z <- res_f_estimtor$z
                results$f_p_value <- res_f$Pvalue.Est # 注意这里的列名 Pvalue.Est 或 Causal.pval
            },
            error = function(e) {
                results$f_error <<- e$message
                warning("Error in MR method F (MR-Egger): ", e$message, call. = FALSE)
            }
        )
    }
    res_f_time_2 <- Sys.time()
    results$f_duration <- as.numeric(difftime(res_f_time_2, res_f_time_1, units = "secs"))


    # g: MRMix
    res_g_time_1 <- Sys.time()
    if (can_run_mr_package_methods) { # MRMix 也需要有效的数据输入
        tryCatch(
            {
                # 确保 MRMix 的输入不包含 NA/Inf，否则它可能会出错
                # 可以在这里添加更严格的检查，或者依赖 MRMix 自身的错误处理
                # 确保向量长度至少为1
                if (length(hat_expose_trMR$results_of_fgwas_beta) < 1) stop("MRMix requires at least one SNP.")

                res_g <- MRMix::MRMix( # 假设 MRMix 包已加载，或者用 MRMix::MRMix
                    betahat_x = hat_expose_trMR$results_of_fgwas_beta,
                    betahat_y = hat_outcome_trMR$results_of_fgwas_beta,
                    sx = hat_expose_trMR$beta_se,
                    sy = hat_outcome_trMR$beta_se
                )
                results$g_theta_point <- res_g$theta
                results$g_theta_se <- res_g$SE_theta
                # MRMix 直接提供 z 和 p，不需要 calculate_mr_stats
                results$g_z <- res_g$zstat_theta
                results$g_p_value <- res_g$pvalue_theta
            },
            error = function(e) {
                results$g_error <<- e$message
                warning("Error in MR method G (MRMix): ", e$message, call. = FALSE)
            }
        )
    }
    res_g_time_2 <- Sys.time()
    results$g_duration <- as.numeric(difftime(res_g_time_2, res_g_time_1, units = "secs"))


    # H: MRBEE (当前被注释掉了，如果启用，也需要 tryCatch)
    # results$H_error <<- "MRBEE not run (commented out)."


    # i：MR-ConMix
    res_i_time_1 <- Sys.time()
    if (can_run_mr_package_methods) {
        tryCatch(
            {
                res_i <- MendelianRandomization::mr_conmix(mr_input_obj)
                results$i_theta_point <- res_i$Estimate
                results$i_p_value <- res_i$Pvalue
                # ConMix 可能不直接提供 SE，如果需要，可能要从其他地方获取或计算
                # results$i_theta_se <- res_i$StdError # 检查 res_i 是否有这个字段
                # res_i_estimtor <- calculate_mr_stats(res_i$Estimate, results$i_theta_se)
                # results$i_z <- res_i_estimtor$z
                # 如果没有 SE 和 Z，可以将它们留为 NA
                results$i_theta_se <- ifelse("StdError" %in% names(res_i), res_i$StdError, NA)
                if (!is.na(results$i_theta_se)) {
                    res_i_estimtor <- calculate_mr_stats(results$i_theta_point, results$i_theta_se)
                    results$i_z <- res_i_estimtor$z
                } else {
                    results$i_z <- NA
                }
            },
            error = function(e) {
                results$i_error <<- e$message
                warning("Error in MR method I (MR-ConMix): ", e$message, call. = FALSE)
            }
        )
    }
    res_i_time_2 <- Sys.time()
    results$i_duration <- as.numeric(difftime(res_i_time_2, res_i_time_1, units = "secs"))

    return(results)
}
