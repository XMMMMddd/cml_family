# %%
# 多核模拟函数

library(parallel)
library(doParallel)
library(foreach)
#' @title 并行运行三重体家族MR模拟 (已修复)
#'
#' @param n_simulations 要运行的总模拟次数。
#' @param n_cores 要使用的CPU核心数。
#' @param script_path 包含 `triplet_family_simulation_once_robust` 函数的 R 脚本的路径。
#' @param seed 用于可重复性研究的随机数种子。
#' @param ... 其他要传递给 `triplet_family_simulation_once_robust` 函数的参数。
#'
#' @return 一个数据框，每一行代表一次模拟的结果。
#'
run_simulation_in_parallel <- function(n_simulations,
                                       n_cores,
                                       script_path = "MR 模拟/MR 模拟（单次）.R",
                                       seed = 123,
                                       ...) {
    message(paste("正在启动", n_cores, "个核心..."))
    cl <- makeCluster(n_cores)
    on.exit(stopCluster(cl))

    clusterSetRNGStream(cl, seed)

    # -------------------- 新增的修复代码在这里 --------------------
    # 将 'script_path' 变量从当前环境导出到所有核心
    clusterExport(cl, "script_path", envir = environment())
    # -----------------------------------------------------------

    message(paste("正在为每个核心加载脚本:", script_path))
    clusterEvalQ(cl, {
        # 现在每个核心都知道 'script_path' 是什么了
        library(dplyr)
        library(MendelianRandomization)
        source(script_path)
    })

    message(paste("正在", n_cores, "个核心上运行", n_simulations, "次模拟..."))
    start_time <- Sys.time()

    results_list <- parLapply(cl, 1:n_simulations, function(i, ...) {
        triplet_family_simulation_once_robust(...)
    }, ...)

    end_time <- Sys.time()
    message(paste("所有模拟完成。总耗时:", round(difftime(end_time, start_time, units = "mins"), 2), "分钟。"))

    message("正在合并结果...")
    results_df <- dplyr::bind_rows(results_list)
    results_df$simulation_id <- 1:nrow(results_df)
    return(results_df)
}
calculate_simulation_summary <- function(simulation_results, true_theta) {
    # 定义方法名称
    methods <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")

    # 初始化结果列表
    summary_stats <- list()

    # --- 数据处理逻辑已更新为适配数据框 ---
    # 不再需要 list-based 的过滤
    total_runs <- nrow(simulation_results)
    cat("总运行次数:", total_runs, "\n")

    for (method in methods) {
        # 构造当前方法所需的列名
        theta_col <- paste0(method, "_theta_point")
        se_col <- paste0(method, "_theta_se")
        pval_col <- paste0(method, "_p_value")
        error_col <- paste0(method, "_error")

        # 筛选有效数据
        # 1. 排除该方法本身报告了错误的运行
        # 2. 排除效应值、标准误、p值为 NA 或 Inf 的情况
        # 3. 确保标准误和p值在有效范围内
        valid_indices <- (
            is.na(simulation_results[[error_col]]) &
                is.finite(simulation_results[[theta_col]]) &
                is.finite(simulation_results[[se_col]]) &
                simulation_results[[se_col]] > 0 &
                is.finite(simulation_results[[pval_col]]) &
                simulation_results[[pval_col]] >= 0 &
                simulation_results[[pval_col]] <= 1
        )

        # 提取经过清理的有效数据向量
        theta_valid <- simulation_results[[theta_col]][valid_indices]
        se_valid <- simulation_results[[se_col]][valid_indices]
        p_valid <- simulation_results[[pval_col]][valid_indices]

        n_valid <- length(theta_valid)
        cat("方法 '", method, "' 的有效结果数量:", n_valid, "\n")

        # --- 以下的计算逻辑和输出结构与您的原始函数完全相同 ---
        if (n_valid > 0) {
            mean_theta <- mean(theta_valid, na.rm = TRUE)
            sd_theta <- sd(theta_valid, na.rm = TRUE)
            median_theta <- median(theta_valid, na.rm = TRUE)
            mean_se <- mean(se_valid, na.rm = TRUE)
            power_005 <- mean(p_valid < 0.05, na.rm = TRUE)
            power_001 <- mean(p_valid < 0.01, na.rm = TRUE)
            ci_95_lower <- theta_valid - 1.96 * se_valid
            ci_95_upper <- theta_valid + 1.96 * se_valid
            coverage_95 <- mean(ci_95_lower <= true_theta & ci_95_upper >= true_theta, na.rm = TRUE)
            ci_99_lower <- theta_valid - 2.576 * se_valid
            ci_99_upper <- theta_valid + 2.576 * se_valid
            coverage_99 <- mean(ci_99_lower <= true_theta & ci_99_upper >= true_theta, na.rm = TRUE)
            bias <- mean_theta - true_theta
            mse <- mean((theta_valid - true_theta)^2, na.rm = TRUE)
            rmse <- sqrt(mse)
        } else {
            mean_theta <- sd_theta <- median_theta <- mean_se <- NA
            power_005 <- power_001 <- coverage_95 <- coverage_99 <- NA
            bias <- mse <- rmse <- NA
        }

        # 存储结果 (保持您原始的列表结构)
        summary_stats[[method]] <- list(
            n_valid = n_valid,
            mean_theta = mean_theta,
            sd_theta = sd_theta,
            median_theta = median_theta,
            power_005 = power_005,
            power_001 = power_001,
            coverage_95 = coverage_95,
            coverage_99 = coverage_99,
            mean_se = mean_se,
            bias = bias,
            mse = mse,
            rmse = rmse
        )
    }

    return(summary_stats)
}
# %% 运行
a <- run_simulation_in_parallel(
    n_simulations = 1000, n_cores = 11,
    script_path = "MR 模拟/MR 模拟（单次）.R",
    n = 10, # 要生成的总数据集数量
    num_pleiotropic = 0,
    n_expose_heterogeneity = 3,
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
    mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0,
    mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0,
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
    beta_exp_to_out = 0.02, # 暴露对结局的真实因果效应
    beta_confounding_exp = 0.2, # 影响暴露的混杂因素的方差 (效应大小为1)
    beta_confounding_out = 0.2, # 影响结局的混杂因素的方差 (效应大小为1)
    correlation = 0.2, # 共享环境因素的方差
    seed = NULL
)

b <- calculate_simulation_summary(a, true_theta = 0.02)
b$g
saveRDS(a, "MR简单模拟结果/都是好的真值为0.rds")
# a <- readRDS("MR简单模拟结果/都是好的SNPs真值为0.02.rds")
# b <- calculate_simulation_summary(a, true_theta = 0.02)