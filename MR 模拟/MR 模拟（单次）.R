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
source("cml家庭函数/cml家庭函数一代/cml_friamly.r")
source("cml家庭函数/cml家庭函数二代/cml_family_version2.r")
# %% 辅助函数

cor2cov_big_matrix <- function(cor_matrix, param_se) {


}
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
n <- 10
num_pleiotropic <- 0
n_expose_heterogeneity <- 0
N_exp <- 1000
N_out <- 1000
overlap_prop <- 0
p_f <- 0.3
p_m <- 0.3
# --- 暴露效应 (在暴露组中) ---
beta_FStoOE_exp <- 0.1
beta_MStoOE_exp <- 0.1
beta_OStoOE_exp <- 0.3
# --- (新) 异质性暴露效应 (在结局组中，当被指定时) ---
beta_FStoOE_exp_out <- 0.05
beta_MStoOE_exp_out <- 0.05
beta_OStoOE_exp_out <- 0.15
# --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
mean_beta_FStoOE_out <- 0.1
sd_beta_FStoOE_out <- 0.05
mean_beta_MStoOE_out <- 0.1
sd_beta_MStoOE_out <- 0.05
mean_beta_OStoOE_out <- 0.1
sd_beta_OStoOE_out <- 0.05
prop_negative_pleiotropy <- 0.5 # 多效性效应为负的比例
# --- 其他模拟参数 ---
assortative_mating_prob <- 0
assortative_mating_strength <- 1000
crowd_stratification_differences <- 0
beta_exp_to_out <- 0.4
confounding_exp <- 0.2 # 参数名应与底层函数匹配
confounding_out <- 0.2 # 参数名应与底层函数匹配
correlation <- 0.2
seed <- NULL
beta_confounding_exp <- 0.2
beta_confounding_out <- 0.2
# %% #!写一个数据生成的流程



# %% 模拟一次
# triplet_family_simulation_once(
#     n_snps = 10, n_pleiotropy = 0,
#     n_independent = 10000, p_trio = 0.9,
#     p_exp_out = 0.5, p_overlap = 0,
#     p_f = 0.3, p_m = 0.3, # p_m 当前未在SNP生成中使用
#     # 暴露效应
#     beta_fs_oe_exp = 3, beta_ms_oe_exp = 3,
#     beta_os_oe_exp = 3,
#     # 结局效应 (直接多效性 / 遗传叠加效应)
#     beta_fs_oe_out = 0, beta_ms_oe_out = 0,
#     beta_os_oe_out = 0, p_negative_pleiotropy = 0,
#     # 因果效应
#     beta_exp_to_out = 0,
#     # 混杂效应
#     var_confounding_exp = 0.2, var_confounding_out = 0.2,
#     # 其他参数
#     r_correlation = 0.2, n_seed = NULL,
#     # 选型婚配强度(跨性状)
#     assortative_mating_strength = 1000
# )

# %% 不模拟样本重叠了用原来的函数
triplet_family_simulation_once <- function(
    n = 10, # 要生成的总数据集数量
    num_pleiotropic = 1,
    n_expose_heterogeneity = 1,
    N_exp = 1000, N_out = 1000,
    p_f = 0.3, p_m = 0.3,
    # --- 暴露效应 (当非零时的大小) ---
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3,
    # --- 暴露异质性 ---
    h_beta_FStoOE_exp = 0.1, h_beta_MStoOE_exp = 0.1,
    h_beta_OStoOE_exp = 0.3,
    # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
    # 子代基因型 -> 结局 (遗传叠加效应或基因多效性)
    mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0.05,
    # 父母基因型 -> 结局 (应为子代自身SNP对结局的效应，即基因多效性; 参数名可能需对应调整)
    mean_beta_OStoOE_out = 0.1, sd_beta_OStoOE_out = 0.05,
    prop_negative_pleiotropy = 0.5, # 在指定为多效性的SNP中，其效应为负的比例 (0 到 1)
    # 选型婚配
    assortative_mating_prob = 0,
    assortative_mating_strength = 0, # 选型婚配对结局的影响因子
    # 人群分层
    ## 定义人群分层的差异(次等位基因频率差异)
    crowd_stratification_differences = 0, # 用于模拟两个具有不同等位基因频率的亚群
    # --- 其他效应 ---
    beta_exp_to_out = 0.4, # 暴露对结局的真实因果效应
    beta_confounding_exp = 0.2, # 影响暴露的混杂因素的方差 (效应大小为1)
    beta_confounding_out = 0.2, # 影响结局的混杂因素的方差 (效应大小为1)
    correlation = 0.2, # 共享环境因素的方差
    seed = NULL) {
    # 首先声明结果
    results <- list(
        a_theta_point = NA, a_theta_se = NA, a_z = NA,
        a_p_value = NA, a_duration = NA,
        b_theta_point = NA, b_theta_se = NA, b_z = NA,
        b_p_value = NA, b_duration = NA,
        c_theta_point = NA, c_theta_se = NA, c_z = NA,
        c_p_value = NA, c_duration = NA,
        d_theta_point = NA, d_theta_se = NA, d_z = NA,
        d_p_value = NA, d_duration = NA,
        e_theta_point = NA, e_theta_se = NA, e_z = NA,
        e_p_value = NA, e_duration = NA,
        f_theta_point = NA, f_theta_se = NA, f_z = NA,
        f_p_value = NA, f_duration = NA,
        g_theta_point = NA, g_theta_se = NA, g_z = NA,
        g_p_value = NA, g_duration = NA,
        h_theta_point = NA, h_theta_se = NA, h_z = NA,
        h_p_value = NA, h_duration = NA,
        i_theta_point = NA, i_theta_se = NA, i_z = NA,
        i_p_value = NA, i_duration = NA
    ) # 初始化结果列表

    # 生成数据
    ## 生成用于计算MR的样本
    phase_two_data_full <- generate_multiple_datasets_v4(
        n = n, # 要生成的总数据集数量
        num_pleiotropic = num_pleiotropic,
        n_expose_heterogeneity = n_expose_heterogeneity,
        N_exp = N_exp, N_out = N_out,
        p_f = p_f, p_m = p_m,
        # --- 暴露效应 (当非零时的大小) ---
        beta_FStoOE_exp = beta_FStoOE_exp, beta_MStoOE_exp = beta_MStoOE_exp,
        beta_OStoOE_exp = beta_OStoOE_exp,
        # --- 暴露异质性 ---
        h_beta_FStoOE_exp = h_beta_FStoOE_exp,
        h_beta_MStoOE_exp = h_beta_MStoOE_exp,
        h_beta_OStoOE_exp = h_beta_OStoOE_exp,
        # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
        # 子代基因型 -> 结局 (遗传叠加效应或基因多效性)
        mean_beta_FStoOE_out = mean_beta_FStoOE_out,
        sd_beta_FStoOE_out = sd_beta_FStoOE_out,
        mean_beta_MStoOE_out = mean_beta_MStoOE_out,
        sd_beta_MStoOE_out = sd_beta_MStoOE_out,
        # 父母基因型 -> 结局 (应为子代自身SNP对结局的效应，即基因多效性; 参数名可能需对应调整)
        mean_beta_OStoOE_out = mean_beta_OStoOE_out,
        sd_beta_OStoOE_out = sd_beta_OStoOE_out,
        prop_negative_pleiotropy = prop_negative_pleiotropy,
        # 选型婚配
        assortative_mating_prob = assortative_mating_prob,
        assortative_mating_strength = assortative_mating_strength,
        # 人群分层
        crowd_stratification_differences = crowd_stratification_differences,
        # --- 其他效应 ---
        beta_exp_to_out = beta_exp_to_out, # 暴露对结局的真实因果效应
        beta_confounding_exp = beta_confounding_exp,
        beta_confounding_out = beta_confounding_out,
        correlation = correlation,
        seed = NULL
    )



    # 常规纳入分析的数据（FGWAS）

    data_exp <- phase_two_data_full$exposure_data
    data_out <- phase_two_data_full$outcome_data
    data_out <- data_out %>% dplyr::select(
        "Father_SNPs",
        "Mother_SNPs",
        "Offspring_SNPs",
        "Father_outcome",
        "Mother_outcome",
        "Offspring_outcome",
        "dataset_id"
    )
    names(data_out) <- gsub("_outcome$", "_expose", names(data_out))
    phase_two_data_analysis_exp <- fgwas_for_data_optimized(
        data_exp
    )
    phase_two_data_analysis_out <- fgwas_for_data_optimized(
        data_out
    )
    # 基于普通的lmm的数据
    phase_two_data_analysis_exp_gls <- fgwas_for_data_optimized_gls(
        data_exp
    )
    phase_two_data_analysis_out_gls <- fgwas_for_data_optimized_gls(
        data_out
    )
    phase_two_data_analysis_exp_mr <- fgwas_to_mr(phase_two_data_analysis_exp)
    phase_two_data_analysis_out_mr <- fgwas_to_mr(phase_two_data_analysis_out)


    # -------- FGWAS为基础 --------


    # a: cml_overlap
    res_a_time_1 <- Sys.time()
    res_a <- mr_cML_Overlap(
        phase_two_data_analysis_exp_mr$results_of_fgwas_beta,
        phase_two_data_analysis_out_mr$results_of_fgwas_beta,
        phase_two_data_analysis_exp_mr$beta_se,
        phase_two_data_analysis_out_mr$beta_se,
        n = N_exp, # 确认样本量参数
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
    beta_sigma_exp <- as.matrix(phase_two_data_analysis_exp$Sigma_inv)
    beta_sigma_out <- as.matrix(phase_two_data_analysis_out$Sigma_inv)
    beta_hat_exp <- as.matrix(phase_two_data_analysis_exp$beta_hat)
    beta_hat_out <- as.matrix(phase_two_data_analysis_out$beta_hat)

    res_b_time_1 <- Sys.time()
    res_b <- cml_family_2_cn_cpp(
        beta_hat_exp, beta_hat_out, beta_sigma_exp, beta_sigma_out
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

    # c: IVW
    ## 做成能用MendelianRandomization这个包导入的数据
    mr_input_obj <- MendelianRandomization::mr_input(
        bx = as.vector(phase_two_data_analysis_exp_mr$results_of_fgwas_beta),
        bxse = as.vector(phase_two_data_analysis_exp_mr$beta_se),
        by = as.vector(phase_two_data_analysis_out_mr$results_of_fgwas_beta),
        byse = as.vector(phase_two_data_analysis_out_mr$beta_se)
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



    # d: MR-Egger
    res_d_time_1 <- Sys.time()
    res_d <- MendelianRandomization::mr_egger(mr_input_obj)
    res_d_time_2 <- Sys.time()
    results$d_theta_point <- res_d$Estimate
    results$d_theta_se <- res_d$StdError.Est
    res_d_estimtor <- calculate_mr_stats(res_d$Estimate, res_d$StdError.Est)
    results$d_z <- res_d_estimtor$z
    results$d_p_value <- res_d$Causal.pval
    results$d_duration <- as.numeric(res_d_time_2 - res_d_time_1)

    # e: MR-Egger
    res_e_time_1 <- Sys.time()
    res_e <- MendelianRandomization::mr_lasso(mr_input_obj)
    res_e_time_2 <- Sys.time()
    results$e_theta_point <- res_e$Estimate
    results$e_theta_se <- res_e$StdError
    res_e_estimtor <- calculate_mr_stats(res_e$Estimate, res_e$StdError)
    results$e_z <- res_e_estimtor$z
    results$e_p_value <- res_e$Pvalue
    results$e_duration <- as.numeric(res_e_time_2 - res_e_time_1)

    # H: MRBEE
    # debug(MRBEE.IMRP) # 开启调试模式
    # fit <- MRBEE.IMRP(
    #     by = matrix(hat_outcome_trMR$results_of_fgwas_beta),
    #     bX = matrix(hat_expose_trMR$results_of_fgwas_beta),
    #     byse = matrix(hat_outcome_trMR$beta_se),
    #     bXse = matrix(hat_expose_trMR$beta_se),
    #     Rxy = diag(1, 2)
    # )

    # -------- 基础LMM为基础 --------
    # f: cml_overlap

    res_f_time_1 <- Sys.time()
    res_f <- mr_cML_Overlap(
        phase_two_data_analysis_exp_gls$beta_hat,
        phase_two_data_analysis_out_gls$beta_hat,
        phase_two_data_analysis_exp_gls$beta_se,
        phase_two_data_analysis_out_gls$beta_se,
        n = N_exp, # 确认样本量参数
        rho = 0
    )
    res_f_time_2 <- Sys.time()
    res_f_estimtor <- calculate_mr_stats(res_f$MA_BIC_theta, res_f$MA_BIC_se)
    results$f_theta_point <- res_f$MA_BIC_theta
    results$f_theta_se <- res_f$MA_BIC_se
    results$f_z <- res_f_estimtor$z
    results$f_p_value <- res_f$MA_BIC_p
    results$f_duration <- as.numeric(res_f_time_2 - res_f_time_1)


    # g: IVW
    ## 做成能用MendelianRandomization这个包导入的数据
    mr_input_obj <- MendelianRandomization::mr_input(
        bx = as.vector(phase_two_data_analysis_exp_gls$beta_hat),
        bxse = as.vector(phase_two_data_analysis_exp_gls$beta_se),
        by = as.vector(phase_two_data_analysis_out_gls$beta_hat),
        byse = as.vector(phase_two_data_analysis_out_gls$beta_se)
    )
    res_g_time_1 <- Sys.time()
    res_g <- MendelianRandomization::mr_ivw(mr_input_obj)
    res_g_time_2 <- Sys.time()
    results$g_theta_point <- res_g$Estimate
    results$g_theta_se <- res_g$StdError
    res_g_estimtor <- calculate_mr_stats(res_g$Estimate, res_g$StdError)
    results$g_z <- res_g_estimtor$z
    results$g_p_value <- res_g$Pvalue
    results$g_duration <- as.numeric(res_g_time_2 - res_g_time_1)



    # d: MR-Egger
    res_h_time_1 <- Sys.time()
    res_h <- MendelianRandomization::mr_egger(mr_input_obj)
    res_h_time_2 <- Sys.time()
    results$h_theta_point <- res_h$Estimate
    results$h_theta_se <- res_h$StdError.Est
    res_h_estimtor <- calculate_mr_stats(res_h$Estimate, res_h$StdError.Est)
    results$h_z <- res_h_estimtor$z
    results$h_p_value <- res_h$Causal.pval
    results$h_duration <- as.numeric(res_h_time_2 - res_h_time_1)

    # e: MR-Egger
    res_i_time_1 <- Sys.time()
    res_i <- MendelianRandomization::mr_lasso(mr_input_obj)
    res_i_time_2 <- Sys.time()
    results$i_theta_point <- res_i$Estimate
    results$i_theta_se <- res_i$StdError
    res_i_estimtor <- calculate_mr_stats(res_i$Estimate, res_i$StdError)
    results$i_z <- res_i_estimtor$z
    results$i_p_value <- res_i$Pvalue
    results$i_duration <- as.numeric(res_i_time_2 - res_i_time_1)




    return(results)
}
triplet_family_simulation_once_robust_old <- function(
    n = 10, # 要生成的总数据集数量
    num_pleiotropic = 1,
    n_expose_heterogeneity = 1,
    N_exp = 1000, N_out = 1000,
    p_f = 0.3, p_m = 0.3,
    # --- 暴露效应 (当非零时的大小) ---
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3,
    # --- 暴露异质性 ---
    h_beta_FStoOE_exp = 0.1, h_beta_MStoOE_exp = 0.1,
    h_beta_OStoOE_exp = 0.3,
    # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
    # 子代基因型 -> 结局 (遗传叠加效应或基因多效性)
    mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0.05,
    # 父母基因型 -> 结局 (应为子代自身SNP对结局的效应，即基因多效性; 参数名可能需对应调整)
    mean_beta_OStoOE_out = 0.1, sd_beta_OStoOE_out = 0.05,
    prop_negative_pleiotropy = 0.5, # 在指定为多效性的SNP中，其效应为负的比例 (0 到 1)
    # 选型婚配
    assortative_mating_prob = 0,
    assortative_mating_strength = 0, # 选型婚配对结局的影响因子
    # 人群分层
    ## 定义人群分层的差异(次等位基因频率差异)
    crowd_stratification_differences = 0, # 用于模拟两个具有不同等位基因频率的亚群
    # --- 其他效应 ---
    beta_exp_to_out = 0.4, # 暴露对结局的真实因果效应
    beta_confounding_exp = 0.2, # 影响暴露的混杂因素的方差 (效应大小为1)
    beta_confounding_out = 0.2, # 影响结局的混杂因素的方差 (效应大小为1)
    correlation = 0.2, # 共享环境因素的方差
    seed = NULL) {
    # 首先声明结果
    results <- list(
        a_theta_point = NA, a_theta_se = NA, a_z = NA,
        a_p_value = NA, a_duration = NA, a_error = NA,
        b_theta_point = NA, b_theta_se = NA, b_z = NA,
        b_p_value = NA, b_duration = NA, b_error = NA,
        c_theta_point = NA, c_theta_se = NA, c_z = NA,
        c_p_value = NA, c_duration = NA, c_error = NA,
        d_theta_point = NA, d_theta_se = NA, d_z = NA,
        d_p_value = NA, d_duration = NA, d_error = NA,
        e_theta_point = NA, e_theta_se = NA, e_z = NA,
        e_p_value = NA, e_duration = NA, e_error = NA,
        f_theta_point = NA, f_theta_se = NA, f_z = NA,
        f_p_value = NA, f_duration = NA, f_error = NA,
        g_theta_point = NA, g_theta_se = NA, g_z = NA,
        g_p_value = NA, g_duration = NA, g_error = NA,
        h_theta_point = NA, h_theta_se = NA, h_z = NA,
        h_p_value = NA, h_duration = NA, h_error = NA,
        i_theta_point = NA, i_theta_se = NA, i_z = NA,
        i_p_value = NA, i_duration = NA, i_error = NA
    ) # 初始化结果列表

    # 辅助函数：安全执行函数并捕获错误
    safe_execute <- function(expr, method_name, results_list, prefix) {
        tryCatch(
            {
                result <- expr
                results_list[[paste0(prefix, "_error")]] <- NA
                return(result)
            },
            error = function(e) {
                warning(paste("Method", method_name, "failed:", e$message))
                results_list[[paste0(prefix, "_error")]] <- e$message
                return(NULL)
            }
        )
    }

    # 生成数据 - 添加错误处理
    tryCatch(
        {
            phase_two_data_full <- generate_multiple_datasets_v4(
                n = n, # 要生成的总数据集数量
                num_pleiotropic = num_pleiotropic,
                n_expose_heterogeneity = n_expose_heterogeneity,
                N_exp = N_exp, N_out = N_out,
                p_f = p_f, p_m = p_m,
                # --- 暴露效应 (当非零时的大小) ---
                beta_FStoOE_exp = beta_FStoOE_exp, beta_MStoOE_exp = beta_MStoOE_exp,
                beta_OStoOE_exp = beta_OStoOE_exp,
                # --- 暴露异质性 ---
                h_beta_FStoOE_exp = h_beta_FStoOE_exp,
                h_beta_MStoOE_exp = h_beta_MStoOE_exp,
                h_beta_OStoOE_exp = h_beta_OStoOE_exp,
                # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
                # 子代基因型 -> 结局 (遗传叠加效应或基因多效性)
                mean_beta_FStoOE_out = mean_beta_FStoOE_out,
                sd_beta_FStoOE_out = sd_beta_FStoOE_out,
                mean_beta_MStoOE_out = mean_beta_MStoOE_out,
                sd_beta_MStoOE_out = sd_beta_MStoOE_out,
                # 父母基因型 -> 结局 (应为子代自身SNP对结局的效应，即基因多效性; 参数名可能需对应调整)
                mean_beta_OStoOE_out = mean_beta_OStoOE_out,
                sd_beta_OStoOE_out = sd_beta_OStoOE_out,
                prop_negative_pleiotropy = prop_negative_pleiotropy,
                # 选型婚配
                assortative_mating_prob = assortative_mating_prob,
                assortative_mating_strength = assortative_mating_strength,
                # 人群分层
                crowd_stratification_differences = crowd_stratification_differences,
                # --- 其他效应 ---
                beta_exp_to_out = beta_exp_to_out, # 暴露对结局的真实因果效应
                beta_confounding_exp = beta_confounding_exp,
                beta_confounding_out = beta_confounding_out,
                correlation = correlation,
                seed = NULL
            )
        },
        error = function(e) {
            stop(paste("Data generation failed:", e$message))
        }
    )

    # 数据预处理 - 添加错误处理
    tryCatch(
        {
            data_exp <- phase_two_data_full$exposure_data
            data_out <- phase_two_data_full$outcome_data
            data_out <- data_out %>% dplyr::select(
                "Father_SNPs",
                "Mother_SNPs",
                "Offspring_SNPs",
                "Father_outcome",
                "Mother_outcome",
                "Offspring_outcome",
                "dataset_id"
            )
            names(data_out) <- gsub("_outcome$", "_expose", names(data_out))

            phase_two_data_analysis_exp <- fgwas_for_data_optimized(data_exp)
            phase_two_data_analysis_out <- fgwas_for_data_optimized(data_out)

            # 基于普通的lmm的数据
            phase_two_data_analysis_exp_gls <- fgwas_for_data_optimized_gls(data_exp)
            phase_two_data_analysis_out_gls <- fgwas_for_data_optimized_gls(data_out)

            phase_two_data_analysis_exp_mr <- fgwas_to_mr(phase_two_data_analysis_exp)
            phase_two_data_analysis_out_mr <- fgwas_to_mr(phase_two_data_analysis_out)
        },
        error = function(e) {
            stop(paste("Data preprocessing failed:", e$message))
        }
    )

    # -------- FGWAS为基础 --------

    # a: cml_overlap
    res_a_time_1 <- Sys.time()
    res_a <- safe_execute(
        {
            mr_cML_Overlap(
                phase_two_data_analysis_exp_mr$results_of_fgwas_beta,
                phase_two_data_analysis_out_mr$results_of_fgwas_beta,
                phase_two_data_analysis_exp_mr$beta_se,
                phase_two_data_analysis_out_mr$beta_se,
                n = N_exp, # 确认样本量参数
                rho = 0
            )
        },
        "CML_Overlap",
        results,
        "a"
    )
    res_a_time_2 <- Sys.time()

    if (!is.null(res_a)) {
        tryCatch(
            {
                res_a_estimtor <- calculate_mr_stats(res_a$MA_BIC_theta, res_a$MA_BIC_se)
                results$a_theta_point <- res_a$MA_BIC_theta
                results$a_theta_se <- res_a$MA_BIC_se
                results$a_z <- res_a_estimtor$z
                results$a_p_value <- res_a$MA_BIC_p
            },
            error = function(e) {
                results$a_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$a_duration <- as.numeric(res_a_time_2 - res_a_time_1)

    # b: cml_family_ver2 方法
    res_b_time_1 <- Sys.time()
    res_b <- safe_execute(
        {
            beta_sigma_exp <- as.matrix(phase_two_data_analysis_exp$Sigma_inv)
            beta_sigma_out <- as.matrix(phase_two_data_analysis_out$Sigma_inv)
            beta_hat_exp <- as.matrix(phase_two_data_analysis_exp$beta_hat)
            beta_hat_out <- as.matrix(phase_two_data_analysis_out$beta_hat)

            cml_family_2_cn_cpp(
                beta_hat_exp, beta_hat_out, beta_sigma_exp, beta_sigma_out
            )
        },
        "CML_Family_v2",
        results,
        "b"
    )
    res_b_time_2 <- Sys.time()

    if (!is.null(res_b)) {
        tryCatch(
            {
                results$b_theta_point <- res_b[[1]]
                results$b_theta_se <- res_b[[2]]
                res_b_estimtor <- calculate_mr_stats(
                    results$b_theta_point,
                    results$b_theta_se
                )
                results$b_z <- res_b_estimtor$z
                results$b_p_value <- res_b_estimtor$p_value
            },
            error = function(e) {
                results$b_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$b_duration <- as.numeric(res_b_time_2 - res_b_time_1)

    # c: IVW
    res_c_time_1 <- Sys.time()
    res_c <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(phase_two_data_analysis_exp_mr$results_of_fgwas_beta),
                bxse = as.vector(phase_two_data_analysis_exp_mr$beta_se),
                by = as.vector(phase_two_data_analysis_out_mr$results_of_fgwas_beta),
                byse = as.vector(phase_two_data_analysis_out_mr$beta_se)
            )
            MendelianRandomization::mr_ivw(mr_input_obj)
        },
        "IVW",
        results,
        "c"
    )
    res_c_time_2 <- Sys.time()

    if (!is.null(res_c)) {
        tryCatch(
            {
                results$c_theta_point <- res_c$Estimate
                results$c_theta_se <- res_c$StdError
                res_c_estimtor <- calculate_mr_stats(res_c$Estimate, res_c$StdError)
                results$c_z <- res_c_estimtor$z
                results$c_p_value <- res_c$Pvalue
            },
            error = function(e) {
                results$c_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$c_duration <- as.numeric(res_c_time_2 - res_c_time_1)

    # d: MR-Egger
    res_d_time_1 <- Sys.time()
    res_d <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(phase_two_data_analysis_exp_mr$results_of_fgwas_beta),
                bxse = as.vector(phase_two_data_analysis_exp_mr$beta_se),
                by = as.vector(phase_two_data_analysis_out_mr$results_of_fgwas_beta),
                byse = as.vector(phase_two_data_analysis_out_mr$beta_se)
            )
            MendelianRandomization::mr_egger(mr_input_obj)
        },
        "MR-Egger",
        results,
        "d"
    )
    res_d_time_2 <- Sys.time()

    if (!is.null(res_d)) {
        tryCatch(
            {
                results$d_theta_point <- res_d$Estimate
                results$d_theta_se <- res_d$StdError.Est
                res_d_estimtor <- calculate_mr_stats(res_d$Estimate, res_d$StdError.Est)
                results$d_z <- res_d_estimtor$z
                results$d_p_value <- res_d$Causal.pval
            },
            error = function(e) {
                results$d_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$d_duration <- as.numeric(res_d_time_2 - res_d_time_1)

    # e: MR-Lasso
    res_e_time_1 <- Sys.time()
    res_e <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(phase_two_data_analysis_exp_mr$results_of_fgwas_beta),
                bxse = as.vector(phase_two_data_analysis_exp_mr$beta_se),
                by = as.vector(phase_two_data_analysis_out_mr$results_of_fgwas_beta),
                byse = as.vector(phase_two_data_analysis_out_mr$beta_se)
            )
            MendelianRandomization::mr_lasso(mr_input_obj)
        },
        "MR-Lasso",
        results,
        "e"
    )
    res_e_time_2 <- Sys.time()

    if (!is.null(res_e)) {
        tryCatch(
            {
                results$e_theta_point <- res_e$Estimate
                results$e_theta_se <- res_e$StdError
                res_e_estimtor <- calculate_mr_stats(res_e$Estimate, res_e$StdError)
                results$e_z <- res_e_estimtor$z
                results$e_p_value <- res_e$Pvalue
            },
            error = function(e) {
                results$e_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$e_duration <- as.numeric(res_e_time_2 - res_e_time_1)

    # -------- 基础LMM为基础 --------

    # f: cml_overlap (GLS版本)
    res_f_time_1 <- Sys.time()
    res_f <- safe_execute(
        {
            mr_cML_Overlap(
                phase_two_data_analysis_exp_gls$beta_hat,
                phase_two_data_analysis_out_gls$beta_hat,
                phase_two_data_analysis_exp_gls$beta_se,
                phase_two_data_analysis_out_gls$beta_se,
                n = N_exp, # 确认样本量参数
                rho = 0
            )
        },
        "CML_Overlap_GLS",
        results,
        "f"
    )
    res_f_time_2 <- Sys.time()

    if (!is.null(res_f)) {
        tryCatch(
            {
                res_f_estimtor <- calculate_mr_stats(res_f$MA_BIC_theta, res_f$MA_BIC_se)
                results$f_theta_point <- res_f$MA_BIC_theta
                results$f_theta_se <- res_f$MA_BIC_se
                results$f_z <- res_f_estimtor$z
                results$f_p_value <- res_f$MA_BIC_p
            },
            error = function(e) {
                results$f_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$f_duration <- as.numeric(res_f_time_2 - res_f_time_1)

    # g: IVW (GLS版本)
    res_g_time_1 <- Sys.time()
    res_g <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(phase_two_data_analysis_exp_gls$beta_hat),
                bxse = as.vector(phase_two_data_analysis_exp_gls$beta_se),
                by = as.vector(phase_two_data_analysis_out_gls$beta_hat),
                byse = as.vector(phase_two_data_analysis_out_gls$beta_se)
            )
            MendelianRandomization::mr_ivw(mr_input_obj)
        },
        "IVW_GLS",
        results,
        "g"
    )
    res_g_time_2 <- Sys.time()

    if (!is.null(res_g)) {
        tryCatch(
            {
                results$g_theta_point <- res_g$Estimate
                results$g_theta_se <- res_g$StdError
                res_g_estimtor <- calculate_mr_stats(res_g$Estimate, res_g$StdError)
                results$g_z <- res_g_estimtor$z
                results$g_p_value <- res_g$Pvalue
            },
            error = function(e) {
                results$g_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$g_duration <- as.numeric(res_g_time_2 - res_g_time_1)

    # h: MR-Egger (GLS版本)
    res_h_time_1 <- Sys.time()
    res_h <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(phase_two_data_analysis_exp_gls$beta_hat),
                bxse = as.vector(phase_two_data_analysis_exp_gls$beta_se),
                by = as.vector(phase_two_data_analysis_out_gls$beta_hat),
                byse = as.vector(phase_two_data_analysis_out_gls$beta_se)
            )
            MendelianRandomization::mr_egger(mr_input_obj)
        },
        "MR-Egger_GLS",
        results,
        "h"
    )
    res_h_time_2 <- Sys.time()

    if (!is.null(res_h)) {
        tryCatch(
            {
                results$h_theta_point <- res_h$Estimate
                results$h_theta_se <- res_h$StdError.Est
                res_h_estimtor <- calculate_mr_stats(res_h$Estimate, res_h$StdError.Est)
                results$h_z <- res_h_estimtor$z
                results$h_p_value <- res_h$Causal.pval
            },
            error = function(e) {
                results$h_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$h_duration <- as.numeric(res_h_time_2 - res_h_time_1)

    # i: MR-Lasso (GLS版本)
    res_i_time_1 <- Sys.time()
    res_i <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(phase_two_data_analysis_exp_gls$beta_hat),
                bxse = as.vector(phase_two_data_analysis_exp_gls$beta_se),
                by = as.vector(phase_two_data_analysis_out_gls$beta_hat),
                byse = as.vector(phase_two_data_analysis_out_gls$beta_se)
            )
            MendelianRandomization::mr_lasso(mr_input_obj)
        },
        "MR-Lasso_GLS",
        results,
        "i"
    )
    res_i_time_2 <- Sys.time()

    if (!is.null(res_i)) {
        tryCatch(
            {
                results$i_theta_point <- res_i$Estimate
                results$i_theta_se <- res_i$StdError
                res_i_estimtor <- calculate_mr_stats(res_i$Estimate, res_i$StdError)
                results$i_z <- res_i_estimtor$z
                results$i_p_value <- res_i$Pvalue
            },
            error = function(e) {
                results$i_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$i_duration <- as.numeric(res_i_time_2 - res_i_time_1)

    return(results)
}

# %% 新的函数
triplet_family_simulation_once_robust <- function(
    n = 10, # 要生成的总数据集数量
    num_pleiotropic = 1,
    n_expose_heterogeneity = 1,
    N_exp = 1000, N_out = 1000,
    p_f = 0.3, p_m = 0.3,
    # --- 暴露效应 (当非零时的大小) ---
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3,
    # --- 暴露异质性 ---
    h_beta_FStoOE_exp = 0.1, h_beta_MStoOE_exp = 0.1,
    h_beta_OStoOE_exp = 0.3,
    # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
    # 子代基因型 -> 结局 (遗传叠加效应或基因多效性)
    mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0.05,
    # 父母基因型 -> 结局 (应为子代自身SNP对结局的效应，即基因多效性; 参数名可能需对应调整)
    mean_beta_OStoOE_out = 0.1, sd_beta_OStoOE_out = 0.05,
    prop_negative_pleiotropy = 0.5, # 在指定为多效性的SNP中，其效应为负的比例 (0 到 1)
    # 选型婚配
    assortative_mating_prob = 0,
    assortative_mating_strength = 0, # 选型婚配对结局的影响因子
    # 人群分层
    ## 定义人群分层的差异(次等位基因频率差异)
    crowd_stratification_differences = 0, # 用于模拟两个具有不同等位基因频率的亚群
    # --- 其他效应 ---
    beta_exp_to_out = 0.4, # 暴露对结局的真实因果效应
    beta_confounding_exp = 0.2, # 影响暴露的混杂因素的方差 (效应大小为1)
    beta_confounding_out = 0.2, # 影响结局的混杂因素的方差 (效应大小为1)
    correlation = 0.2, # 共享环境因素的方差
    seed = NULL) {
    # 首先声明结果
    results <- list(
        a_theta_point = NA, a_theta_se = NA, a_z = NA,
        a_p_value = NA, a_duration = NA, a_error = NA,
        b_theta_point = NA, b_theta_se = NA, b_z = NA,
        b_p_value = NA, b_duration = NA, b_error = NA,
        c_theta_point = NA, c_theta_se = NA, c_z = NA,
        c_p_value = NA, c_duration = NA, c_error = NA,
        d_theta_point = NA, d_theta_se = NA, d_z = NA,
        d_p_value = NA, d_duration = NA, d_error = NA,
        e_theta_point = NA, e_theta_se = NA, e_z = NA,
        e_p_value = NA, e_duration = NA, e_error = NA,
        f_theta_point = NA, f_theta_se = NA, f_z = NA,
        f_p_value = NA, f_duration = NA, f_error = NA,
        g_theta_point = NA, g_theta_se = NA, g_z = NA,
        g_p_value = NA, g_duration = NA, g_error = NA,
        h_theta_point = NA, h_theta_se = NA, h_z = NA,
        h_p_value = NA, h_duration = NA, h_error = NA,
        i_theta_point = NA, i_theta_se = NA, i_z = NA,
        i_p_value = NA, i_duration = NA, i_error = NA
    ) # 初始化结果列表

    # 辅助函数：安全执行函数并捕获错误
    safe_execute <- function(expr, method_name, results_list, prefix) {
        tryCatch(
            {
                result <- expr
                results_list[[paste0(prefix, "_error")]] <- NA
                return(result)
            },
            error = function(e) {
                warning(paste("Method", method_name, "failed:", e$message))
                results_list[[paste0(prefix, "_error")]] <- e$message
                return(NULL)
            }
        )
    }

    # 生成数据 - 添加错误处理
    tryCatch(
        {
            phase_two_data_full <- generate_multiple_datasets_v4(
                n = n, # 要生成的总数据集数量
                num_pleiotropic = num_pleiotropic,
                n_expose_heterogeneity = n_expose_heterogeneity,
                N_exp = N_exp, N_out = N_out,
                p_f = p_f, p_m = p_m,
                # --- 暴露效应 (当非零时的大小) ---
                beta_FStoOE_exp = beta_FStoOE_exp, beta_MStoOE_exp = beta_MStoOE_exp,
                beta_OStoOE_exp = beta_OStoOE_exp,
                # --- 暴露异质性 ---
                h_beta_FStoOE_exp = h_beta_FStoOE_exp,
                h_beta_MStoOE_exp = h_beta_MStoOE_exp,
                h_beta_OStoOE_exp = h_beta_OStoOE_exp,
                # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
                # 子代基因型 -> 结局 (遗传叠加效应或基因多效性)
                mean_beta_FStoOE_out = mean_beta_FStoOE_out,
                sd_beta_FStoOE_out = sd_beta_FStoOE_out,
                mean_beta_MStoOE_out = mean_beta_MStoOE_out,
                sd_beta_MStoOE_out = sd_beta_MStoOE_out,
                # 父母基因型 -> 结局 (应为子代自身SNP对结局的效应，即基因多效性; 参数名可能需对应调整)
                mean_beta_OStoOE_out = mean_beta_OStoOE_out,
                sd_beta_OStoOE_out = sd_beta_OStoOE_out,
                prop_negative_pleiotropy = prop_negative_pleiotropy,
                # 选型婚配
                assortative_mating_prob = assortative_mating_prob,
                assortative_mating_strength = assortative_mating_strength,
                # 人群分层
                crowd_stratification_differences = crowd_stratification_differences,
                # --- 其他效应 ---
                beta_exp_to_out = beta_exp_to_out, # 暴露对结局的真实因果效应
                beta_confounding_exp = beta_confounding_exp,
                beta_confounding_out = beta_confounding_out,
                correlation = correlation,
                seed = NULL
            )
        },
        error = function(e) {
            stop(paste("Data generation failed:", e$message))
        }
    )

    # 数据预处理 - 添加错误处理
    tryCatch(
        {
            data_exp <- phase_two_data_full$exposure_data %>%
                dplyr::select(-ends_with("outcome"))
            data_out <- phase_two_data_full$outcome_data
            data_out <- data_out %>% dplyr::select(
                "Father_SNPs",
                "Mother_SNPs",
                "Offspring_SNPs",
                "Father_outcome",
                "Mother_outcome",
                "Offspring_outcome",
                "dataset_id"
            )
            names(data_out) <- gsub("_outcome$", "_expose", names(data_out))

# FGWAS
            phase_two_data_analysis_exp <-
                fgwas_for_data_optimized(data_exp,
                    processing_func = FMR_trio_optimized
                )
            phase_two_data_analysis_out <-
                fgwas_for_data_optimized(data_out,
                    processing_func = FMR_trio_optimized
                )

            # 基于普通的lmm的数据（父母没有当协变量）
            phase_two_data_analysis_exp_gls <-
                fgwas_for_data_optimized_single(data_exp,
                    processing_func = FMR_trio_optimized_onlyc3
                )
            phase_two_data_analysis_out_gls <-
                fgwas_for_data_optimized_single(data_out,
                    processing_func = FMR_trio_optimized_onlyc3
                )
            # 处理成能够输入的
            beta_hat_exp_lmm <- phase_two_data_analysis_exp_gls$beta_hat
            se_hat_exp_lmm <- phase_two_data_analysis_exp_gls$beta_se

            beta_hat_out_lmm <- phase_two_data_analysis_out_gls$beta_hat
            se_hat_out_lmm <- phase_two_data_analysis_out_gls$beta_se
            # 基于只分析孩子父母当协变量，只分析孩子
            phase_two_data_analysis_exp_mr <- fgwas_to_mr(phase_two_data_analysis_exp)
            phase_two_data_analysis_out_mr <- fgwas_to_mr(phase_two_data_analysis_out)
        },
        error = function(e) {
            stop(paste("Data preprocessing failed:", e$message))
        }
    )

    # -------- FGWAS为基础 --------

    # a: cml_overlap
    res_a_time_1 <- Sys.time()
    res_a <- safe_execute(
        {
            mr_cML_Overlap(
                phase_two_data_analysis_exp_mr$results_of_fgwas_beta,
                phase_two_data_analysis_out_mr$results_of_fgwas_beta,
                phase_two_data_analysis_exp_mr$beta_se,
                phase_two_data_analysis_out_mr$beta_se,
                n = N_exp, # 确认样本量参数
                rho = 0
            )
        },
        "CML_Overlap",
        results,
        "a"
    )
    res_a_time_2 <- Sys.time()

    if (!is.null(res_a)) {
        tryCatch(
            {
                res_a_estimtor <- calculate_mr_stats(
                    res_a$MA_BIC_theta,
                    res_a$MA_BIC_se
                )
                results$a_theta_point <- res_a$MA_BIC_theta
                results$a_theta_se <- res_a$MA_BIC_se
                results$a_z <- res_a_estimtor$z
                results$a_p_value <- res_a$MA_BIC_p
            },
            error = function(e) {
                results$a_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$a_duration <- as.numeric(res_a_time_2 - res_a_time_1)

    # b: cml_family_ver2 方法
    res_b_time_1 <- Sys.time()
    res_b <- safe_execute(
        {
            beta_sigma_exp <- as.matrix(phase_two_data_analysis_exp$Sigma_inv)
            beta_sigma_out <- as.matrix(phase_two_data_analysis_out$Sigma_inv)
            beta_hat_exp <- as.matrix(phase_two_data_analysis_exp$beta_hat)
            beta_hat_out <- as.matrix(phase_two_data_analysis_out$beta_hat)

            cml_family_multistart(
                beta_hat_exp, beta_hat_out, beta_sigma_exp, beta_sigma_out,
                n = N_exp, cov_index = "normal"
            )
        },
        "CML_Family_v2",
        results,
        "b"
    )
    res_b_time_2 <- Sys.time()

    if (!is.null(res_b)) {
        tryCatch(
            {
                results$b_theta_point <- res_b[[2]]$weighted_alpha
                results$b_theta_se <- res_b[[2]]$weighted_se
                res_b_estimtor <- calculate_mr_stats(
                    results$b_theta_point,
                    results$b_theta_se
                )
                results$b_z <- res_b_estimtor$z
                results$b_p_value <- res_b_estimtor$p_value
            },
            error = function(e) {
                results$b_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$b_duration <- as.numeric(res_b_time_2 - res_b_time_1)

    # c: IVW
    res_c_time_1 <- Sys.time()
    res_c <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(phase_two_data_analysis_exp_mr$results_of_fgwas_beta),
                bxse = as.vector(phase_two_data_analysis_exp_mr$beta_se),
                by = as.vector(phase_two_data_analysis_out_mr$results_of_fgwas_beta),
                byse = as.vector(phase_two_data_analysis_out_mr$beta_se)
            )
            MendelianRandomization::mr_ivw(mr_input_obj)
        },
        "IVW",
        results,
        "c"
    )
    res_c_time_2 <- Sys.time()

    if (!is.null(res_c)) {
        tryCatch(
            {
                results$c_theta_point <- res_c$Estimate
                results$c_theta_se <- res_c$StdError
                res_c_estimtor <- calculate_mr_stats(res_c$Estimate, res_c$StdError)
                results$c_z <- res_c_estimtor$z
                results$c_p_value <- res_c$Pvalue
            },
            error = function(e) {
                results$c_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$c_duration <- as.numeric(res_c_time_2 - res_c_time_1)

    # d: MR-Egger
    res_d_time_1 <- Sys.time()
    res_d <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(phase_two_data_analysis_exp_mr$results_of_fgwas_beta),
                bxse = as.vector(phase_two_data_analysis_exp_mr$beta_se),
                by = as.vector(phase_two_data_analysis_out_mr$results_of_fgwas_beta),
                byse = as.vector(phase_two_data_analysis_out_mr$beta_se)
            )
            MendelianRandomization::mr_egger(mr_input_obj)
        },
        "MR-Egger",
        results,
        "d"
    )
    res_d_time_2 <- Sys.time()

    if (!is.null(res_d)) {
        tryCatch(
            {
                results$d_theta_point <- res_d$Estimate
                results$d_theta_se <- res_d$StdError.Est
                res_d_estimtor <- calculate_mr_stats(res_d$Estimate, res_d$StdError.Est)
                results$d_z <- res_d_estimtor$z
                results$d_p_value <- res_d$Causal.pval
            },
            error = function(e) {
                results$d_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$d_duration <- as.numeric(res_d_time_2 - res_d_time_1)

    # e: MR-Lasso
    res_e_time_1 <- Sys.time()
    res_e <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(phase_two_data_analysis_exp_mr$results_of_fgwas_beta),
                bxse = as.vector(phase_two_data_analysis_exp_mr$beta_se),
                by = as.vector(phase_two_data_analysis_out_mr$results_of_fgwas_beta),
                byse = as.vector(phase_two_data_analysis_out_mr$beta_se)
            )
            MendelianRandomization::mr_lasso(mr_input_obj)
        },
        "MR-Lasso",
        results,
        "e"
    )
    res_e_time_2 <- Sys.time()

    if (!is.null(res_e)) {
        tryCatch(
            {
                results$e_theta_point <- res_e$Estimate
                results$e_theta_se <- res_e$StdError
                res_e_estimtor <- calculate_mr_stats(res_e$Estimate, res_e$StdError)
                results$e_z <- res_e_estimtor$z
                results$e_p_value <- res_e$Pvalue
            },
            error = function(e) {
                results$e_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$e_duration <- as.numeric(res_e_time_2 - res_e_time_1)

    # -------- 基础LMM为基础 --------




    # f: cml_overlap (GLS版本)
    res_f_time_1 <- Sys.time()
    res_f <- safe_execute(
        {
            mr_cML_Overlap(
                beta_hat_exp_lmm,
                beta_hat_out_lmm,
                se_hat_exp_lmm,
                se_hat_out_lmm,
                n = N_exp, # 确认样本量参数
                rho = 0
            )
        },
        "CML_Overlap_GLS",
        results,
        "f"
    )
    res_f_time_2 <- Sys.time()

    if (!is.null(res_f)) {
        tryCatch(
            {
                res_f_estimtor <- calculate_mr_stats(
                    res_f$MA_BIC_theta,
                    res_f$MA_BIC_se
                )
                results$f_theta_point <- res_f$MA_BIC_theta
                results$f_theta_se <- res_f$MA_BIC_se
                results$f_z <- res_f_estimtor$z
                results$f_p_value <- res_f$MA_BIC_p
            },
            error = function(e) {
                results$f_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$f_duration <- as.numeric(res_f_time_2 - res_f_time_1)



    # g: IVW (GLS版本)
    res_g_time_1 <- Sys.time()
    res_g <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(beta_hat_exp_lmm),
                bxse = as.vector(se_hat_exp_lmm),
                by = as.vector(beta_hat_out_lmm),
                byse = as.vector(se_hat_out_lmm)
            )
            MendelianRandomization::mr_ivw(mr_input_obj)
        },
        "IVW_GLS",
        results,
        "g"
    )
    res_g_time_2 <- Sys.time()

    if (!is.null(res_g)) {
        tryCatch(
            {
                results$g_theta_point <- res_g$Estimate
                results$g_theta_se <- res_g$StdError
                res_g_estimtor <- calculate_mr_stats(res_g$Estimate, res_g$StdError)
                results$g_z <- res_g_estimtor$z
                results$g_p_value <- res_g$Pvalue
            },
            error = function(e) {
                results$g_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$g_duration <- as.numeric(res_g_time_2 - res_g_time_1)



    # h: MR-Egger (GLS版本)
    res_h_time_1 <- Sys.time()
    res_h <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(beta_hat_exp_lmm),
                bxse = as.vector(se_hat_exp_lmm),
                by = as.vector(beta_hat_out_lmm),
                byse = as.vector(se_hat_out_lmm)
            )
            MendelianRandomization::mr_egger(mr_input_obj)
        },
        "MR-Egger_GLS",
        results,
        "h"
    )
    res_h_time_2 <- Sys.time()

    if (!is.null(res_h)) {
        tryCatch(
            {
                results$h_theta_point <- res_h$Estimate
                results$h_theta_se <- res_h$StdError.Est
                res_h_estimtor <- calculate_mr_stats(res_h$Estimate, res_h$StdError.Est)
                results$h_z <- res_h_estimtor$z
                results$h_p_value <- res_h$Causal.pval
            },
            error = function(e) {
                results$h_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$h_duration <- as.numeric(res_h_time_2 - res_h_time_1)




    # i: MR-Lasso (GLS版本)
    res_i_time_1 <- Sys.time()
    res_i <- safe_execute(
        {
            mr_input_obj <- MendelianRandomization::mr_input(
                bx = as.vector(beta_hat_exp_lmm),
                bxse = as.vector(se_hat_exp_lmm),
                by = as.vector(beta_hat_out_lmm),
                byse = as.vector(se_hat_out_lmm)
            )
            MendelianRandomization::mr_lasso(mr_input_obj)
        },
        "MR-Lasso_GLS",
        results,
        "i"
    )
    res_i_time_2 <- Sys.time()

    if (!is.null(res_i)) {
        tryCatch(
            {
                results$i_theta_point <- res_i$Estimate
                results$i_theta_se <- res_i$StdError
                res_i_estimtor <- calculate_mr_stats(res_i$Estimate, res_i$StdError)
                results$i_z <- res_i_estimtor$z
                results$i_p_value <- res_i$Pvalue
            },
            error = function(e) {
                results$i_error <<- paste("Post-processing error:", e$message)
            }
        )
    }
    results$i_duration <- as.numeric(res_i_time_2 - res_i_time_1)

    return(results)
}



# %% 模拟
if (FALSE) {
    phase_two_data_full <- generate_multiple_datasets_v4(
        n = 10, num_pleiotropic = 0,
        n_expose_heterogeneity = 0,
        beta_FStoOE_exp = 0,
        beta_MStoOE_exp = 0,
        beta_OStoOE_exp = 0.3,
        h_beta_FStoOE_exp = 0,
        h_beta_MStoOE_exp = 0,
        h_beta_OStoOE_exp = 0.3,
        assortative_mating_strength = 1000,
        crowd_stratification_differences = 0,
        beta_exp_to_out = 1.2,
        beta_confounding_exp = 0.2,
        beta_confounding_out = 0.2,
        correlation = 0.2, seed = NULL
    )

    # 常规纳入分析的数据（FGWAS）

    data_exp <- phase_two_data_full$exposure_data
    data_out <- phase_two_data_full$outcome_data
    data_out <- data_out %>% dplyr::select(
        "Father_SNPs",
        "Mother_SNPs",
        "Offspring_SNPs",
        "Father_outcome",
        "Mother_outcome",
        "Offspring_outcome",
        "dataset_id"
    )
    names(data_out) <- gsub("_outcome$", "_expose", names(data_out))
    phase_two_data_analysis_exp <- fgwas_for_data_optimized(
        data_exp
    )
    phase_two_data_analysis_out <- fgwas_for_data_optimized(
        data_out
    )
    # 基于普通的lmm的数据
    phase_two_data_analysis_exp_gls <- fgwas_for_data_optimized_gls(
        data_exp
    )
    phase_two_data_analysis_out_gls <- fgwas_for_data_optimized_gls(
        data_out
    )
    phase_two_data_analysis_exp_mr <- fgwas_to_mr(phase_two_data_analysis_exp)
    phase_two_data_analysis_out_mr <- fgwas_to_mr(phase_two_data_analysis_out)


    # -------- FGWAS为基础 --------

    # b:  cml_family_ver2 方法
    beta_sigma_exp <- as.matrix(phase_two_data_analysis_exp$Sigma_inv)
    beta_sigma_out <- as.matrix(phase_two_data_analysis_out$Sigma_inv)
    beta_hat_exp <- as.matrix(phase_two_data_analysis_exp$beta_hat)
    beta_hat_out <- as.matrix(phase_two_data_analysis_out$beta_hat)

    res_b <- cml_family_2_cn_cpp(
        beta_hat_exp, beta_hat_out, beta_sigma_exp, beta_sigma_out
    )
    res_b$best_result$theta
    res_b$best_result$theta_se
    round(arrange(res_b$bic_results, bic), 3)
    2 * pnorm(abs(res_b$theta_weight / res_b$theta_se_weight),
        lower.tail = FALSE
    )
}

if (FALSE) {
    a <- triplet_family_simulation_once_robust(
        n = 10, # 要生成的总数据集数量
        num_pleiotropic = 0,
        n_expose_heterogeneity = 0,
        N_exp = 1000, N_out = 1000,
        p_f = 0.3, p_m = 0.3,
        # --- 暴露效应 (当非零时的大小) ---
        beta_FStoOE_exp = 0, beta_MStoOE_exp = 0,
        beta_OStoOE_exp = 0.3,
        # --- 暴露异质性 ---
        h_beta_FStoOE_exp = 0.5, h_beta_MStoOE_exp = 0.5,
        h_beta_OStoOE_exp = 0.3,
        # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
        # 子代基因型 -> 结局 (遗传叠加效应或基因多效性)
        mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0.05,
        mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0.05,
        # 父母基因型 -> 结局 (应为子代自身SNP对结局的效应，即基因多效性; 参数名可能需对应调整)
        mean_beta_OStoOE_out = 0, sd_beta_OStoOE_out = 0,
        prop_negative_pleiotropy = 0, # 在指定为多效性的SNP中，其效应为负的比例 (0 到 1)
        # 选型婚配
        assortative_mating_prob = 0,
        assortative_mating_strength = 1000, # 选型婚配对结局的影响因子
        # 人群分层
        ## 定义人群分层的差异(次等位基因频率差异)
        crowd_stratification_differences = 0, # 用于模拟两个具有不同等位基因频率的亚群
        # --- 其他效应 ---
        beta_exp_to_out = 0, # 暴露对结局的真实因果效应
        beta_confounding_exp = 0.2, # 影响暴露的混杂因素的方差 (效应大小为1)
        beta_confounding_out = 0.2, # 影响结局的混杂因素的方差 (效应大小为1)
        correlation = 0.2, # 共享环境因素的方差
        seed = NULL
    )

    a$b_p_value
    a
}
