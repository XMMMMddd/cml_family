# %% 加载必要的函数
library(dplyr)
library(MASS)
library(nlme) 
library(parallel)
library(pbapply)
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
source("cml家庭函数/cml家庭函数二代/cml_family_version2.r")
source("cml家庭函数/cml家庭函数overlap版本/cml_family_overlap2.r")
source("cml家庭函数/cml家庭函数一代/cml_friamly.r")
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
    test_1 <- cml_family_2_cn_cpp(
        phase_two_data_analysis$beta_hat_exp, phase_two_data_analysis$beta_hat_out,
        phase_two_data_analysis$beta_sigma_exp, phase_two_data_analysis$beta_sigma_out
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
    beta_sigma_out <- as.matrix(phase_two_data_analysis$beta_sigma_out)
    beta_hat_exp <- as.matrix(phase_two_data_analysis$beta_hat_exp)
    beta_hat_out <- as.matrix(phase_two_data_analysis$beta_hat_out)
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
triplet_family_simulation_once(n_independent = 3000, beta_exp_to_out = 0)

# %%
triplet_family_simulation_once_1 <- function(
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

    # a: cml_1

    res_a_time_1 <- Sys.time()
    res_a <- cml_family_2_cn_cpp(
        phase_two_data_analysis$beta_hat_exp,
        phase_two_data_analysis$beta_hat_out,
        phase_two_data_analysis$beta_sigma_exp,
        phase_two_data_analysis$beta_sigma_out
    )
    res_a_time_2 <- Sys.time()
    results$a_theta_point <- res_a[[1]]
    results$a_theta_se <- res_a[[2]]
    res_b_estimtor <- calculate_mr_stats(
        results$a_theta_point,
        results$a_theta_se
    )
    results$a_z <- res_a_estimtor$z
    results$a_p_value <- res_a_estimtor$p_value
    results$a_duration <- as.numeric(res_a_time_2 - res_a_time_1)


    # b:  cml_family_ver2 方法
    beta_sigma_exp <- as.matrix(phase_two_data_analysis$beta_sigma_exp)
    beta_sigma_out <- as.matrix(phase_two_data_analysis$beta_sigma_out)
    beta_hat_exp <- as.matrix(phase_two_data_analysis$beta_hat_exp)
    beta_hat_out <- as.matrix(phase_two_data_analysis$beta_hat_out)
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
    # c ivw
    test_2 <- summary(lm(
        offspring_expose ~
            offspring_snps.2,
        data =
            phase_two_data_full$data_independent_exp
    ))
    print(test_2)
    test_1 <- lm(
        offspring_outcome ~
            offspring_snps.2,
        data =
            phase_two_data_full$data_independent_out
    )

    alpha <- test_1$coefficient[2] / test_2$coefficient[2]
    res_c <- cml_family_ver_2_cpp(
        phase_two_data_analysis$beta_hat_exp,
        phase_two_data_analysis$beta_hat_out,
        phase_two_data_analysis$beta_sigma_exp,
        phase_two_data_analysis$beta_sigma_out
    )


    return(list(res_a, res_b, alpha, res_c))
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
        beta_fs_oe_out = beta_fs_oe_out,
        beta_ms_oe_out = beta_ms_oe_out,
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
    phase_two_data_analysis <- perform_fgwas_analysis_lmm(
        phase_two_data_full,
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

    # a: cml_1

    res_a_time_1 <- Sys.time()
    res_a <- cml_family_2_cn_cpp(
        phase_two_data_analysis$beta_hat_exp,
        phase_two_data_analysis$beta_hat_out,
        phase_two_data_analysis$beta_sigma_exp,
        phase_two_data_analysis$beta_sigma_out
    )
    res_a_time_2 <- Sys.time()
    results$a_theta_point <- res_a[[1]]
    results$a_theta_se <- res_a[[2]]
    res_a_estimtor <- calculate_mr_stats(
        results$a_theta_point,
        results$a_theta_se
    )
    results$a_z <- res_a_estimtor$z
    results$a_p_value <- res_a_estimtor$p_value
    results$a_duration <- as.numeric(res_a_time_2 - res_a_time_1)


    # b:  cml_family_ver2 方法
    beta_sigma_exp <- as.matrix(phase_two_data_analysis$beta_sigma_exp)
    beta_sigma_out <- as.matrix(phase_two_data_analysis$beta_sigma_out)
    beta_hat_exp <- as.matrix(phase_two_data_analysis$beta_hat_exp)
    beta_hat_out <- as.matrix(phase_two_data_analysis$beta_hat_out)
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


    res_c <- cml_family_ver_2_cpp(
        phase_two_data_analysis$beta_hat_exp,
        phase_two_data_analysis$beta_hat_out,
        phase_two_data_analysis$beta_sigma_exp,
        phase_two_data_analysis$beta_sigma_out
    )
    # d: ivw
    mr_input_obj <- MendelianRandomization::mr_input(
        bx = phase_two_data_analysis_exp$results_of_fgwas_beta,
        bxse = phase_two_data_analysis_exp$beta_se,
        by = phase_two_data_analysis_out$results_of_fgwas_beta,
        byse = phase_two_data_analysis_out$beta_se
    )
    res_d_time_1 <- Sys.time()
    res_d <- MendelianRandomization::mr_ivw(mr_input_obj)


    return(list(res_a, res_b, res_c, res_d))
}
triplet_family_simulation_once_3 <- function(
    n = 10, num_pleiotropic = 1, N_exp = 1000,
    N_out = 1000, overlap_prop = 0, p_f = 0.3,
    p_m = 0.3, beta_FStoOE_exp = 0.1,
    beta_MStoOE_exp = 0.1, beta_OStoOE_exp = 0.3,
    mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0.05,
    mean_beta_OStoOE_out = 0.1, sd_beta_OStoOE_out = 0.05,
    prop_negative_pleiotropy = 0.5, assortative_mating_prob = 0,
    assortative_mating_strength = 0,
    crowd_stratification_differences = 0,
    beta_exp_to_out = 0.4, beta_confounding_exp = 0.2,
    beta_confounding_out = 0.2,
    correlation = 0.2, seed = NULL) {
    # 首先声明结果
    results <- list(
        a_theta_point = NA, a_theta_se = NA, a_z = NA,
        a_p_value = NA, a_duration = NA,
        b_theta_point = NA, b_theta_se = NA, b_z = NA,
        b_p_value = NA, b_duration = NA
    ) # 初始化结果列表

    # 生成数据
    phase_two_data_full <- generate_multiple_datasets_v3(
        n = n, num_pleiotropic = num_pleiotropic, N_exp = N_exp,
        N_out = N_out, overlap_prop = overlap_prop, p_f = p_f,
        p_m = p_m, beta_FStoOE_exp = beta_FStoOE_exp,
        beta_MStoOE_exp = beta_MStoOE_exp, beta_OStoOE_exp = beta_OStoOE_exp,
        mean_beta_FStoOE_out = mean_beta_FStoOE_out,
        sd_beta_FStoOE_out = sd_beta_FStoOE_out,
        mean_beta_MStoOE_out = mean_beta_MStoOE_out,
        sd_beta_MStoOE_out = sd_beta_MStoOE_out,
        mean_beta_OStoOE_out = mean_beta_OStoOE_out,
        sd_beta_OStoOE_out = sd_beta_OStoOE_out,
        prop_negative_pleiotropy = prop_negative_pleiotropy,
        assortative_mating_prob = assortative_mating_prob,
        assortative_mating_strength = assortative_mating_strength,
        crowd_stratification_differences = crowd_stratification_differences,
        beta_exp_to_out = beta_exp_to_out,
        beta_confounding_exp = beta_confounding_exp,
        beta_confounding_out = beta_confounding_out,
        correlation = correlation, seed = seed
    )

    phase_two_data_analysis_exp_trio <- fgwas_for_data_optimized(
        phase_two_data_full$exposure_data
    )

    phase_two_data_analysis_exp <- list(
        beta_hat =
            phase_two_data_analysis_exp_trio$beta_hat,
        sigma_inv =
            phase_two_data_analysis_exp_trio$Sigma_inv
    )
    phase_two_data_analysis_exp <- fgwas_to_mr(
        phase_two_data_analysis_exp
    )
    phase_two_data_analysis_out_trio <- fgwas_for_data_optimized(
        phase_two_data_full$outcome_data
    )

    phase_two_data_analysis_out <- list(
        beta_hat =
            phase_two_data_analysis_out_trio$beta_hat,
        sigma_inv =
            phase_two_data_analysis_out_trio$Sigma_inv
    )
    phase_two_data_analysis_out <- fgwas_to_mr(
        phase_two_data_analysis_out
    )

    # a: cml_1

    res_a_time_1 <- Sys.time()
    res_a <- cml_family_2_cn_cpp(
        phase_two_data_analysis_exp_trio$beta_hat,
        phase_two_data_analysis_out_trio$beta_hat,
        phase_two_data_analysis_exp_trio$Sigma_inv,
        phase_two_data_analysis_out_trio$Sigma_inv
    )
    res_a_time_2 <- Sys.time()
    results$a_theta_point <- res_a[[1]]
    results$a_theta_se <- res_a[[2]]
    res_a_estimtor <- calculate_mr_stats(
        results$a_theta_point,
        results$a_theta_se
    )
    results$a_z <- res_a_estimtor$z
    results$a_p_value <- res_a_estimtor$p_value
    results$a_duration <- as.numeric(res_a_time_2 - res_a_time_1)


    # b:  cml_family_ver2 方法
    beta_sigma_exp <- as.matrix(phase_two_data_analysis_exp_trio$Sigma_inv)
    beta_sigma_out <- as.matrix(phase_two_data_analysis_out_trio$Sigma_inv)
    beta_hat_exp <- as.matrix(phase_two_data_analysis_exp_trio$beta_hat)
    beta_hat_out <- as.matrix(phase_two_data_analysis_out_trio$beta_hat)
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


    res_c <- cml_family_ver_2_cpp(
        phase_two_data_analysis_exp_trio$beta_hat,
        phase_two_data_analysis_out_trio$beta_hat,
        phase_two_data_analysis_exp_trio$Sigma_inv,
        phase_two_data_analysis_out_trio$Sigma_inv
    )
    # d: ivw
    mr_input_obj <- MendelianRandomization::mr_input(
        bx = phase_two_data_analysis_exp$results_of_fgwas_beta,
        bxse = phase_two_data_analysis_exp$beta_se,
        by = phase_two_data_analysis_out$results_of_fgwas_beta,
        byse = phase_two_data_analysis_out$beta_se
    )
    res_d_time_1 <- Sys.time()
    res_d <- MendelianRandomization::mr_ivw(mr_input_obj)


    return(list(res_a, res_b, res_c, res_d))
}

# %% 情况2

# 1. 运行您的模拟函数 (这部分保持不变)
a <- triplet_family_simulation_once_2(
    n_snps = 3, n_pleiotropy = 0,
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
    assortative_mating_strength = 1000
)

# 2. 分别计算五种方法的统计量

# 方法一: Weighted Theta
est_1 <- a[[1]]$theta_weight
se_1 <- a[[1]]$theta_se_weight
var_1 <- se_1^2
pval_1 <- 2 * pnorm(abs(est_1 / se_1), lower.tail = FALSE)

# 方法二: Weighted Alpha
est_2 <- a[[2]]$weighted_alpha
se_2 <- a[[2]]$weighted_se
var_2 <- se_2^2
pval_2 <- 2 * pnorm(abs(est_2 / se_2), lower.tail = FALSE)

# 方法三: "Truth" Theta (来自 BIC 结果的特定行)
theta_truth <- a[[1]]$bic_results %>% filter(k_x == 0, k_par == 0)
if (nrow(theta_truth) > 0) {
    est_3 <- theta_truth$theta
    se_3 <- theta_truth$theta_se
    var_3 <- se_3^2
    pval_3 <- 2 * pnorm(abs(est_3 / se_3), lower.tail = FALSE)
} else {
    est_3 <- NA
    se_3 <- NA
    var_3 <- NA
    pval_3 <- NA
}

# 方法四: New Method (来自 a[[3]])
# (修正了注释，使其与代码 a[[3]] 保持一致)
if (!is.null(a[[3]]) && length(a[[3]]) >= 2) {
    est_4 <- a[[3]][[1]]
    se_4 <- a[[3]][[2]]
    var_4 <- se_4^2
    pval_4 <- 2 * pnorm(abs(est_4 / se_4), lower.tail = FALSE)
} else {
    est_4 <- NA
    se_4 <- NA
    var_4 <- NA
    pval_4 <- NA
}

# *** 方法五: IVW (来自 MendelianRandomization 包的结果) ***
# (增加了稳健性检查)
if (isS4(a[[4]]) && all(c("Estimate", "StdError") %in% slotNames(a[[4]]))) {
    est_5 <- a[[4]]@Estimate
    # *** 修正了笔误：原代码为 se_4，现已更正为 se_5 ***
    se_5 <- a[[4]]@StdError
    var_5 <- se_5^2
    pval_5 <- 2 * pnorm(abs(est_5 / se_5), lower.tail = FALSE)
} else {
    est_5 <- NA
    se_5 <- NA
    var_5 <- NA
    pval_5 <- NA
}


# 3. 将所有结果汇总到一个比较表格中 (***已更新为5种方法***)
comparison_results <- tibble(
    `分析方法 (Method)` = c(
        "1. Weighted Theta", "2. Weighted Alpha",
        "3. Truth Theta (k_x=0, k_par=0)",
        "4. New Method",
        "5. IVW (MendelianRandomization)" # <-- 新增方法五名称
    ),
    `点估计 (Estimate)` = c(est_1, est_2, est_3, est_4, est_5),
    `方差 (Variance)` = c(var_1, var_2, var_3, var_4, var_5),
    `标准误 (SE)` = c(se_1, se_2, se_3, se_4, se_5),
    `P值 (P-value)` = c(pval_1, pval_2, pval_3, pval_4, pval_5)
)

# 4. 增加显著性标记列 (保持不变)
comparison_results_with_flag <- comparison_results %>%
    mutate(
        `是否显著 (p<0.05)` = if_else(`P值 (P-value)` < 0.05, "是 (*)", "否", missing = "否 (NA P值)")
    )

# 5. 打印最终的完整表格 (保持不变)
print(comparison_results_with_flag, width = Inf)
# %% 情况3


# 1. 运行您的模拟函数 (这部分保持不变)
a <- triplet_family_simulation_once_3(
    n = 10, num_pleiotropic = 0, N_exp = 1000,
    N_out = 1000, overlap_prop = 0, p_f = 0.3,
    p_m = 0.3, beta_FStoOE_exp = 0.1,
    beta_MStoOE_exp = 0.1, beta_OStoOE_exp = 0.3,
    mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0.05,
    mean_beta_OStoOE_out = 0.1, sd_beta_OStoOE_out = 0.05,
    prop_negative_pleiotropy = 0.5, assortative_mating_prob = 0,
    assortative_mating_strength = 0,
    crowd_stratification_differences = 0,
    beta_exp_to_out = 0, beta_confounding_exp = 0.2,
    beta_confounding_out = 0.2,
    correlation = 0.2, seed = NULL
)

# 2. 分别计算五种方法的统计量

# 方法一: Weighted Theta
est_1 <- a[[1]]$theta_weight
se_1 <- a[[1]]$theta_se_weight
var_1 <- se_1^2
pval_1 <- 2 * pnorm(abs(est_1 / se_1), lower.tail = FALSE)

# 方法二: Weighted Alpha
est_2 <- a[[2]]$weighted_alpha
se_2 <- a[[2]]$weighted_se
var_2 <- se_2^2
pval_2 <- 2 * pnorm(abs(est_2 / se_2), lower.tail = FALSE)

# 方法三: "Truth" Theta (来自 BIC 结果的特定行)
theta_truth <- a[[1]]$bic_results %>% filter(k_x == 0, k_par == 0)
if (nrow(theta_truth) > 0) {
    est_3 <- theta_truth$theta
    se_3 <- theta_truth$theta_se
    var_3 <- se_3^2
    pval_3 <- 2 * pnorm(abs(est_3 / se_3), lower.tail = FALSE)
} else {
    est_3 <- NA
    se_3 <- NA
    var_3 <- NA
    pval_3 <- NA
}

# 方法四: New Method (来自 a[[3]])
# (修正了注释，使其与代码 a[[3]] 保持一致)
if (!is.null(a[[3]]) && length(a[[3]]) >= 2) {
    est_4 <- a[[3]][[1]]
    se_4 <- a[[3]][[2]]
    var_4 <- se_4^2
    pval_4 <- 2 * pnorm(abs(est_4 / se_4), lower.tail = FALSE)
} else {
    est_4 <- NA
    se_4 <- NA
    var_4 <- NA
    pval_4 <- NA
}

# *** 方法五: IVW (来自 MendelianRandomization 包的结果) ***
# (增加了稳健性检查)
if (isS4(a[[4]]) && all(c("Estimate", "StdError") %in% slotNames(a[[4]]))) {
    est_5 <- a[[4]]@Estimate
    # *** 修正了笔误：原代码为 se_4，现已更正为 se_5 ***
    se_5 <- a[[4]]@StdError
    var_5 <- se_5^2
    pval_5 <- 2 * pnorm(abs(est_5 / se_5), lower.tail = FALSE)
} else {
    est_5 <- NA
    se_5 <- NA
    var_5 <- NA
    pval_5 <- NA
}


# 3. 将所有结果汇总到一个比较表格中 (***已更新为5种方法***)
comparison_results <- tibble(
    `分析方法 (Method)` = c(
        "1. Weighted Theta", "2. Weighted Alpha",
        "3. Truth Theta (k_x=0, k_par=0)",
        "4. New Method",
        "5. IVW (MendelianRandomization)" # <-- 新增方法五名称
    ),
    `点估计 (Estimate)` = c(est_1, est_2, est_3, est_4, est_5),
    `方差 (Variance)` = c(var_1, var_2, var_3, var_4, var_5),
    `标准误 (SE)` = c(se_1, se_2, se_3, se_4, se_5),
    `P值 (P-value)` = c(pval_1, pval_2, pval_3, pval_4, pval_5)
)

# 4. 增加显著性标记列 (保持不变)
comparison_results_with_flag <- comparison_results %>%
    mutate(
        `是否显著 (p<0.05)` = if_else(`P值 (P-value)` < 0.05, "是 (*)", "否", missing = "否 (NA P值)")
    )

# 5. 打印最终的完整表格 (保持不变)
print(comparison_results_with_flag, width = Inf)

# %% 循环20次模拟情形2

# --- 准备工作: 加载所有需要的包 ---
# 如果没有安装，请先运行 install.packages(c("dplyr", "nlme", "tibble", "progress", "purrr"))
suppressPackageStartupMessages({
    library(dplyr)
    library(nlme)
    library(tibble)
    library(progress)
    library(purrr)
    library(lme4)
})


perform_one_simulation_run <- function() {
    # 1. 运行模拟函数，并用 tryCatch 捕获可能发生的错误
    a <- tryCatch(
        {
            triplet_family_simulation_once_2(
                n_snps = 3, n_pleiotropy = 0, n_independent = 4000,
                p_trio = 0.9, p_exp_out = 0.5,
                beta_fs_oe_exp = 0, beta_ms_oe_exp = 0,
                beta_os_oe_exp = 0.6,
                beta_exp_to_out = 0, beta_fs_oe_out = 0,
                beta_ms_oe_out = 0, beta_os_oe_out = 0
            )
        },
        error = function(e) {
            # 如果模拟出错，打印警告信息并返回NULL
            warning("模拟函数内部出错: ", e$message)
            return(NULL)
        }
    )

    # 如果模拟失败 (返回了NULL)，则直接退出本次函数运行
    if (is.null(a)) {
        return(NULL)
    }

    # 2. 分别计算四种方法的统计量
    # 方法一: Weighted Theta
    est_1 <- a[[1]]$theta_weight
    se_1 <- a[[1]]$theta_se_weight
    var_1 <- se_1^2
    pval_1 <- 2 * pnorm(abs(est_1 / se_1), lower.tail = FALSE)

    # 方法二: Weighted Alpha
    est_2 <- a[[2]]$weighted_alpha
    se_2 <- a[[2]]$weighted_se
    var_2 <- se_2^2
    pval_2 <- 2 * pnorm(abs(est_2 / se_2), lower.tail = FALSE)

    # 方法三: "Truth" Theta
    theta_truth <- a[[1]]$bic_results %>% filter(k_x == 0, k_par == 0)
    if (nrow(theta_truth) > 0) {
        est_3 <- theta_truth$theta
        se_3 <- theta_truth$theta_se
        var_3 <- se_3^2
        pval_3 <- 2 * pnorm(abs(est_3 / se_3), lower.tail = FALSE)
    } else {
        est_3 <- NA
        se_3 <- NA
        var_3 <- NA
        pval_3 <- NA
    }

    # 方法四: New Method
    if (!is.null(a[[3]]) && length(a[[3]]) >= 2) {
        est_4 <- a[[3]][[1]]
        se_4 <- a[[3]][[2]]
        var_4 <- se_4^2
        pval_4 <- 2 * pnorm(abs(est_4 / se_4), lower.tail = FALSE)
    } else {
        est_4 <- NA
        se_4 <- NA
        var_4 <- NA
        pval_4 <- NA
    }

    # 3. 将所有结果汇总到一个比较表格中
    comparison_results <- tibble(
        `分析方法 (Method)` = c("1. Weighted Theta", "2. Weighted Alpha", "3. Truth Theta (k_x=0, k_par=0)", "4. New Method"),
        `点估计 (Estimate)` = c(est_1, est_2, est_3, est_4),
        `方差 (Variance)` = c(var_1, var_2, var_3, var_4),
        `标准误 (SE)` = c(se_1, se_2, se_3, se_4),
        `P值 (P-value)` = c(pval_1, pval_2, pval_3, pval_4)
    )

    # 4. 增加显著性标记列
    comparison_results_with_flag <- comparison_results %>%
        mutate(
            `是否显著 (p<0.05)` = if_else(`P值 (P-value)` < 0.05, "是 (*)", "否", missing = "否 (NA P值)")
        )

    return(comparison_results_with_flag)
}


# --- 1. 设置模拟参数 ---
n_simulations <- 1000 # 您可以在这里设定想要运行的总次数

# --- 2. 初始化进度条 ---
pb <- progress_bar$new(
    format = "模拟进度 [:bar] :percent | 耗时: :elapsed | 预计剩余: :eta",
    total = n_simulations,
    clear = FALSE,
    width = 80
)

# --- 3. 运行多次模拟 ---
# lapply 会重复调用函数，并将每次的结果存入列表 all_results
all_results <- lapply(1:n_simulations, function(i) {
    pb$tick() # 进度条前进一格
    result <- perform_one_simulation_run()
    # 如果当次运行成功，则为其添加运行次序号
    if (!is.null(result)) {
        result <- result %>% mutate(`运行次序 (Run)` = i, .before = 1)
    }
    return(result)
})

# --- 4. 筛选并报告结果 ---
# 移除那些因出错而返回NULL的无效结果
valid_results <- all_results[!sapply(all_results, is.null)]

# 检查是否所有运行都失败了
if (length(valid_results) == 0) {
    cat("\n❌ 警告：所有模拟运行均未能成功返回有效结果，请检查 `perform_one_simulation_run` 函数内部。\n")
} else {
    # 从有效结果中，筛选出P值小于0.05的运行
    significant_runs <- Filter(function(tbl) {
        any(tbl$`P值 (P-value)` < 0.05, na.rm = TRUE)
    }, valid_results)

    # 打印最终的总结报告
    if (length(significant_runs) > 0) {
        cat(sprintf("\n✅ 在 %d 次有效模拟中，共发现 %d 次运行出现了显著P值 (p < 0.05)。详情如下：\n\n", length(valid_results), length(significant_runs)))
        # 逐个打印筛选出来的表格
        walk(significant_runs, ~ {
            cat(sprintf("--- 结果来自第 %d 次运行 ---\n", .x$`运行次序 (Run)`[1]))
            print(.x, width = Inf)
            cat("\n")
        })
    } else {
        cat(sprintf("\n✅ 在 %d 次有效模拟中，没有发现任何P值小于0.05的情况。\n", length(valid_results)))
    }
}


52/976
