# %%
# 多核模拟函数
library(parallel)
library(doParallel)
library(foreach)
# 载入必要的R包
# 请确保已经安装了这些包: install.packages(c("parallel", "dplyr"))
library(parallel)
library(dplyr)

#' @title (已更新) 并行运行三重奏家庭孟德尔随机化模拟
#'
#' @description
#' 该函数使用并行处理来多次运行模拟。它通过在每个核心加载一个指定的 R 脚本文件
#' 来获取所需的函数。此版本更易于维护。
#'
#' @param n_simulations 要运行的模拟总次数。
#' @param source_file 一个字符串，表示包含所有必需函数（如
#'   'triplet_family_simulation_once_robust' 及其依赖项）的 R 脚本文件的路径。
#' @param ... 其他所有要传递给 'triplet_family_simulation_once_robust' 函数的参数。
#' @param n_cores 要使用的CPU核心数。默认为 `detectCores() - 1`。
#'
#' @return 一个数据框 (tibble)，其格式与 'calculate_simulation_summary' 兼容。
#'
run_simulation_parallel <- function(n_simulations, source_file, ..., n_cores = NULL) {
    # 1. 检查源文件是否存在
    if (!file.exists(source_file)) {
        stop("指定的源文件不存在: ", source_file)
    }

    # 2. 设置并行环境
    if (is.null(n_cores)) {
        n_cores <- detectCores() - 1
        cat("未指定核心数，将使用 ", n_cores, " 个核心。\n")
    }
    cl <- makeCluster(n_cores)
    cat("并行集群已启动...\n")

    # 3. 将必要的对象和库导出到每个核心
    tryCatch(
        {
            # 为确保路径的稳健性，获取文件的绝对路径
            absolute_path_to_source <- normalizePath(source_file)

            # 将“绝对路径”这个变量本身导出到每个核心
            clusterExport(cl, "absolute_path_to_source", envir = environment())

            # 现在，让每个核心使用这个绝对路径来 source 文件，并加载库
            clusterEvalQ(cl, {
                # source() 命令现在在每个核心上运行
                source(absolute_path_to_source)
                # 同样加载所需的库
                library(MendelianRandomization)
            })
            cat("源文件 '", basename(absolute_path_to_source), "' 和所需库已在每个核心加载。\n", sep = "")
        },
        error = function(e) {
            stopCluster(cl)
            stop("在并行核心上加载源文件时发生错误: ", e$message)
        }
    )

    # 4. 为每次模拟生成独立的随机种子
    set.seed(Sys.time())
    seeds <- sample.int(1e9, n_simulations)

    # 5. 捕获要传递给模拟函数的通用参数
    sim_args <- list(...)

    # 6. 使用 parLapply 在多核上执行模拟
    cat("开始并行运行 ", n_simulations, " 次模拟...\n")
    results_list <- parLapply(cl, 1:n_simulations, function(i) {
        current_args <- c(sim_args, list(seed = seeds[i]))
        tryCatch(
            {
                do.call(triplet_family_simulation_once_robust, current_args)
            },
            error = function(e) {
                list(simulation_error = as.character(e))
            }
        )
    })
    cat("所有模拟运行完成。\n")

    # 7. 关闭并行集群
    stopCluster(cl)
    cat("并行集群已关闭。\n")

    # 8. 处理并整合结果 (此部分逻辑与之前完全相同)
    cat("正在处理和整合结果...\n")
    simulation_results_df <- bind_rows(results_list)

    methods <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
    for (method in methods) {
        theta_col <- paste0(method, "_theta_point")
        error_col <- paste0(method, "_error")
        if (theta_col %in% names(simulation_results_df)) {
            simulation_results_df[[error_col]] <- ifelse(
                is.na(simulation_results_df[[theta_col]]),
                "Error: Result was NA",
                NA_character_
            )
        } else {
            simulation_results_df[[error_col]] <- "Error: Result column not found"
        }
    }

    cat("结果处理完成，返回数据框。\n")
    return(simulation_results_df)
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
a <- run_simulation_parallel(
    n_simulations = 1000,
    source_file = "MR 模拟/MR 模拟（单次）.R",
    n = 1000,
    n_snps = 5,
    n_pleiotropic = 0,
    n_expose_heterogeneity = 0,
    p_f = 0.3,
    p_m = 0.3,
    beta_fs_to_oe_exp = 0.3,
    beta_ms_to_oe_exp = 0.3,
    beta_os_to_oe_exp = 0.3,
    h_beta_fs_to_oe_exp = 0.2,
    h_beta_ms_to_oe_exp = 0.1,
    h_beta_os_to_oe_exp = 0.3,
    mean_beta_fs_to_oe_out = 0.1,
    sd_beta_fs_to_oe_out = 0.05,
    mean_beta_ms_to_oe_out = 0.1,
    sd_beta_ms_to_oe_out = 0.05,
    mean_beta_os_to_oe_out = 0.1,
    sd_beta_os_to_oe_out = 0.05,
    p_neg_pleiotropy = 0.5,
    assortative_mating_strength = 0,
    crowd_differences = 0,
    beta_exp_to_out = 0,
    confounding_exp = 0.2,
    confounding_out = 0.2,
    correlation = 0.2
)
a$e_p_value
b <- calculate_simulation_summary(a, true_theta = 0)
b$b
b$c
saveRDS(a, "MR简单模拟结果/都是好的真值为0.rds")
# a <- readRDS("MR简单模拟结果/都是好的SNPs真值为0.02.rds")
# b <- calculate_simulation_summary(a, true_theta = 0.02)
