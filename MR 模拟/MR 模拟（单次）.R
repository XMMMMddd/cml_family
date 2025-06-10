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
# n_snps <- 3
# n_null_snps <- 10
# n_pleiotropy <- 1
# n_independent <- 1000
# p_trio <- 0.5
# p_exp_out <- 0.5
# p_overlap <- 0
# p_f <- 0.3
# p_m <- 0.3 # p_m 当前未在SNP生成中使用
# # 暴露效应
# beta_fs_oe_exp <- 0.9
# beta_ms_oe_exp <- 0.9
# beta_os_oe_exp <- 0.9
# # 结局效应 (直接多效性 / 遗传叠加效应)
# beta_fs_oe_out <- 0
# beta_ms_oe_out <- 0
# beta_os_oe_out <- 0
# p_negative_pleiotropy <- 0
# # 因果效应
# beta_exp_to_out <- 0.4
# # 混杂效应
# var_confounding_exp <- 0.2
# var_confounding_out <- 0.2
# # 其他参数
# r_correlation <- 0.2
# n_seed <- NULL
# # 选型婚配强度(跨性状)
# assortative_mating_strength <- 1000
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
    phase_two_data_analysis <- perform_fgwas_analysis_lmm(phase_two_data_full,
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
    res_b <- cml_family_overlap_2_robust(
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



# %% 模拟样本重叠
triplet_family_simulation_once_overlap <- function(
    n_snps = 3, n_pleiotropy = 1, n_null_snps = 50,
    n_independent = 3000, p_trio = 0.5,
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
        c_p_value = NA, c_duration = NA,
        d_theta_point = NA, d_theta_se = NA, d_z = NA,
        d_p_value = NA, d_duration = NA
    ) # 初始化结果列表


    # 生成数据
    ## 生成用于计算MR的样本
    phase_two_data_full <- generate_mr_trio_data_ultra_overlap(
        n_snps = n_snps, n_pleiotropy = n_pleiotropy,
        n_null_snps = n_null_snps,
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




    #' @title 从数据框列表中筛选、排序并重命名SNP列（优化版）
    #' @description 根据SNP索引筛选列，并对SNP列重新编号，使其连续。
    #'              此优化版本旨在提高效率和代码清晰度。
    #' @param data_list 一个列表，其中每个元素都是一个数据框。
    #' @param snp_indices 一个数值型向量，包含需要被提取的SNP的编号。
    #' @return 一个新的列表，其中每个数据框都只包含被筛选和重命名后的列。
    select_snps_from_list <- function(data_list, snp_indices) {
        # --- 输入验证 ---
        if (!is.list(data_list)) stop("错误: 'data_list' 必须是一个列表。")
        if (!is.numeric(snp_indices)) stop("错误: 'snp_indices' 必须是一个数值型向量。")
        if (length(data_list) == 0 || length(snp_indices) == 0) {
            warning("输入列表或SNP索引为空，返回原始列表。")
            return(data_list)
        }

        # 构建要匹配的 SNP 列名的正则表达式
        snp_pattern <- paste0("_(snps)\\.", snp_indices, "$", collapse = "|")

        lapply(data_list, function(df) {
            if (!is.data.frame(df)) {
                warning("列表中的一个元素不是数据框，已跳过。")
                return(df)
            }

            all_cols <- colnames(df)

            # 步骤 1: 识别非SNP列和被选中的SNP列
            is_snp_col <- grepl("_snps\\.", all_cols)
            non_snp_cols <- all_cols[!is_snp_col]

            is_selected_snp <- grepl(snp_pattern, all_cols)
            selected_snp_cols <- all_cols[is_selected_snp]

            # 如果没有选中的SNP，则只保留非SNP列
            if (length(selected_snp_cols) == 0) {
                return(df[, non_snp_cols, drop = FALSE])
            }

            # 步骤 2: 对选中的SNP列进行排序以确保重命名的一致性
            # 按前缀（如 "offspring_snps"）和数字后缀排序
            prefixes <- sub("\\..*$", "", selected_snp_cols)
            suffixes <- as.numeric(sub(".*\\.", "", selected_snp_cols))
            sorted_indices <- order(prefixes, suffixes)
            sorted_selected_snp_cols <- selected_snp_cols[sorted_indices]

            # 步骤 3: 筛选数据框，并保持列的原始分组（非SNP列在前）
            df_selected <- df[, c(non_snp_cols, sorted_selected_snp_cols), drop = FALSE]

            # 步骤 4: 高效地重命名排好序的SNP列
            # 使用 rle (run-length encoding) 来处理每个前缀的连续编号
            sorted_prefixes <- sub("\\..*$", "", sorted_selected_snp_cols)
            rle_prefixes <- rle(sorted_prefixes)

            # 为每个前缀生成 1:N 的序列
            new_suffixes <- unlist(lapply(rle_prefixes$lengths, seq_len))

            # 组合成新的列名
            new_snp_names <- paste(sorted_prefixes, new_suffixes, sep = ".")

            # 更新数据框的列名
            colnames(df_selected)[(length(non_snp_cols) + 1):ncol(df_selected)] <- new_snp_names

            return(df_selected)
        })
    }



    # 常规纳入分析的数据
    input_mr <- which(phase_two_data_full$snp_info$type != "Unrelated")
    phase_two_data_input_mr <- select_snps_from_list(
        phase_two_data_full, input_mr
    )
    phase_two_data_analysis <- perform_fgwas_analysis_lmm(
        phase_two_data_input_mr,
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
    # 用于计算相关性的数据
    input_overlap <- which(
        phase_two_data_full$snp_info$type == "Unrelated"
    )
    phase_two_data_overlap <- select_snps_from_list(
        phase_two_data_full, input_overlap
    )

    phase_two_data_overlap <- perform_fgwas_analysis_lmm(
        phase_two_data_overlap,
        n_snps = n_null_snps
    )

    phase_two_data_analysis_exp_overlap <- list(
        beta_hat =
            phase_two_data_overlap$beta_hat_exp,
        sigma_inv =
            phase_two_data_overlap$beta_sigma_exp
    )

    phase_two_data_analysis_out_overlap <- list(
        beta_hat =
            phase_two_data_overlap$beta_hat_out,
        sigma_inv =
            phase_two_data_overlap$beta_sigma_out
    )

    cor_overlap_data <- bind_cols(
        phase_two_data_analysis_exp_overlap$beta_hat,
        phase_two_data_analysis_out_overlap$beta_hat
    )
    # 得到了相关系数矩阵
    overlap_cor <- cor(cor_overlap_data)

    # 优化相关性向量计算：向量化操作替代for循环
    # 提取所有相关的SE向量
    se_exp <- phase_two_data_analysis_exp$beta_se
    se_p_exp <- phase_two_data_analysis_exp$beta_se_p
    se_out <- phase_two_data_analysis_out$beta_se
    se_p_out <- phase_two_data_analysis_out$beta_se_p

    # 向量化计算协方差: cov(X,Y) = cor(X,Y) * se(X) * se(Y)
    # overlap_cor[1,3] 是 cor(beta_exp, beta_out)
    # overlap_cor[2,3] 是 cor(beta_p_exp, beta_out)
    # overlap_cor[1,4] 是 cor(beta_exp, beta_p_out)
    # overlap_cor[2,4] 是 cor(beta_p_exp, beta_p_out)
    cov_13 <- overlap_cor[1, 3] * se_exp * se_out
    cov_23 <- overlap_cor[2, 3] * se_p_exp * se_out
    cov_14 <- overlap_cor[1, 4] * se_exp * se_p_out
    cov_24 <- overlap_cor[2, 4] * se_p_exp * se_p_out

    # 将四个向量交错合并成一个 cor_vector
    # rbind 按行绑定，然后 as.vector 会按列读取，实现交错效果
    cor_vector <- as.vector(rbind(cov_13, cov_23, cov_14, cov_24))


    # a: cml_overlap
    res_a_time_1 <- Sys.time()
    res_a <- mr_cML_Overlap(
        phase_two_data_analysis_exp$results_of_fgwas_beta,
        phase_two_data_analysis_out$results_of_fgwas_beta,
        phase_two_data_analysis_exp$beta_se,
        phase_two_data_analysis_out$beta_se,
        n = floor(n_independent * p_trio), # 确认样本量参数
        rho = overlap_cor[1, 3]
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
        off_diagonal_elements = cor_vector
    )
    res_b_time_1 <- Sys.time()
    res_b <- cml_family_overlap_2_robust(
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



    # f: MR-Egger
    res_d_time_1 <- Sys.time()
    res_d <- MendelianRandomization::mr_egger(mr_input_obj)
    res_d_time_2 <- Sys.time()
    results$d_theta_point <- res_d$Estimate
    results$d_theta_se <- res_d$StdError.Est
    res_d_estimtor <- calculate_mr_stats(res_d$Estimate, res_d$StdError.Est)
    results$d_z <- res_d_estimtor$z
    results$d_p_value <- res_d$Causal.pval
    results$d_duration <- as.numeric(res_d_time_2 - res_d_time_1)


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

# triplet_family_simulation_once_overlap(
#     n_snps = 5, n_pleiotropy = 0, n_null_snps = 10, n_independent = 1000, p_trio = 0.9,
#     p_exp_out = 0.5, p_overlap = 1, p_f = 0.3, p_m = 0.3,
#     beta_fs_oe_exp = 0, beta_ms_oe_exp = 0, beta_os_oe_exp = 1.5, beta_fs_oe_out = 0, beta_ms_oe_out = 0, beta_os_oe_out = 0, p_negative_pleiotropy = 0, beta_exp_to_out = 0,
#     var_confounding_exp = 0.2, var_confounding_out = 0.2, r_correlation = 0.2, n_seed = NULL, assortative_mating_strength = 1000
# )
