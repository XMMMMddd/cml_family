# %%
# install.packages("doParallel")
library(doParallel)
library(purrr)
source("MR 模拟/MR 模拟（单次）.R")

# %% 模拟函数定义
simulation_cmlfamily_multsnps_cpp <- function(
    n_simulations, param) {
    # 检测可用核数
    num_cores <- detectCores() - 1

    # 设置集群
    cluster <- makeCluster(num_cores)
    clusterEvalQ(cluster, {
        source("MR 模拟/MR 模拟（单次）.R")
    })


    # 并行化计算
    results <- parLapply(
        cl = cluster, # 使用命名参数 cl
        X = 1:n_simulations, # 迭代序列
        # --- 匿名函数仍然需要接收 sim_index ---
        fun = function(sim_index, sim_params_list) {
            # 第一个参数 sim_index 接收迭代序号，但后面不用它
            # 调用内部模拟函数
            triplet_family_simulation_once(
                # --- 从列表解包参数 ---
                n = sim_params_list$n,
                num_pleiotropic = sim_params_list$num_pleiotropic,
                N_out = sim_params_list$N_out,
                # 样本重叠
                overlap_prop = sim_params_list$overlap_prop,
                # 直接效应
                beta_FStoOE_exp = sim_params_list$beta_FStoOE_exp,
                beta_MStoOE_exp = sim_params_list$beta_MStoOE_exp,
                beta_OStoOE_exp = sim_params_list$beta_OStoOE_exp,
                # 水平多效性
                prop_negative_pleiotropy =
                    sim_params_list$prop_negative_pleiotropy, # 不满足inside假设
                mean_beta_FStoOE_out = sim_params_list$mean_beta_FStoOE_out,
                sd_beta_FStoOE_out = sim_params_list$sd_beta_FStoOE_out,
                mean_beta_MStoOE_out = sim_params_list$mean_beta_MStoOE_out,
                sd_beta_MStoOE_out = sim_params_list$sd_beta_MStoOE_out,
                mean_beta_OStoOE_out = sim_params_list$mean_beta_OStoOE_out,
                sd_beta_OStoOE_out = sim_params_list$sd_beta_OStoOE_out,
                # 选型婚配
                ## 选型婚配基因型
                compatibility_selection_prop =
                    sim_params_list$compatibility_selection_prop,
                compatibility_selection_geno =
                    sim_params_list$compatibility_selection_geno,
                correlation_param = sim_params_list$correlation_param,
                ## 选型婚配暴露相关
                compatibility_selection_factor_exp =
                    sim_params_list$compatibility_selection_factor_exp,
                compatibility_selection_factor_out =
                    sim_params_list$compatibility_selection_factor_out,
                # 人群分层（双人群差异）
                crowd_stratification_differences =
                    sim_params_list$crowd_stratification_differences,
                # 其他参数设置
                beta_exp_to_out = sim_params_list$beta_exp_to_out,
                beta_confounding_exp = sim_params_list$beta_confounding_exp,
                beta_confounding_out = sim_params_list$beta_confounding_out,
                correlation = sim_params_list$correlation,
                seed = NULL
            )
        },
        # --- 将 'param' 列表传递给匿名函数的 'sim_params_list' 参数 ---
        sim_params_list = param
    )

    # 关闭集群
    stopCluster(cluster)


    return(results)
}
# list 转 tibble
simulation_ListtoTibble <- function(simulation_output_list) {
    # 1. 过滤掉可能的错误结果 (如果 parLapply 返回了包含 error 的列表)
    #    假设错误结果是带有 'error' 元素的列表，或者干脆就是 NULL/非列表
    valid_results <- Filter(function(res) is.list(res) && is.null(res$error), simulation_output_list)

    n_total <- length(simulation_output_list)
    n_valid <- length(valid_results)

    if (n_valid == 0) {
        warning("在 ", n_total, " 次模拟中没有找到有效的结果来创建 tibble。")
        return(list(tibble = NULL, n_valid = 0))
    }

    if (n_valid < n_total) {
        warning("在 ", n_total, " 次模拟中，有 ", n_total - n_valid, " 个结果无效或包含错误，已被忽略。")
    }

    # 2. 使用 dplyr::bind_rows() 合并列表中的所有命名列表
    #    bind_rows 会自动将列表元素的名称作为列名
    result_df <- tryCatch(
        {
            dplyr::bind_rows(valid_results)
        },
        error = function(e) {
            warning(paste("使用 dplyr::bind_rows 合并结果时出错:", e$message))
            # 尝试备用方法（可能较慢，且对复杂类型如矩阵列处理不同）
            warning("正在尝试使用 do.call(rbind, ...) 作为备选方案...")
            tryCatch(
                {
                    # 需要先将每个列表转换为单行数据框，确保类型兼容
                    valid_results_df <- lapply(valid_results, function(x) {
                        # 注意：如果列表包含矩阵 (如 B_best_r)，这可能会变复杂或失败
                        # 我们暂时假设 B_best_r 可以被忽略或后续处理
                        # 移除矩阵列或其他复杂类型以简化 rbind
                        x_simple <- x[!sapply(x, is.matrix)]
                        as.data.frame(x_simple)
                    })
                    do.call(rbind, valid_results_df)
                },
                error = function(e2) {
                    warning(paste("使用备选方案 do.call(rbind, ...) 也失败:", e2$message))
                    return(NULL) # 最终失败
                }
            )
        }
    )

    # 3. 检查结果并返回
    if (is.null(result_df)) {
        warning("无法将模拟结果列表转换为数据框。")
        return(list(tibble = NULL, n_valid = n_valid)) # 虽然有有效结果，但转换失败
    } else {
        # 确保返回的是 tibble 或 data.frame
        return(list(tibble = tibble::as_tibble(result_df), n_valid = n_valid))
    }
}
# %% 小型模拟
sim_params_list <- list( # --- 从列表解包参数 ---
    n = 10, num_pleiotropic = 0, N_out = 1000,
    # 样本重叠
    overlap_prop = 0,
    # 直接效应
    beta_FStoOE_exp = 0.3,
    beta_MStoOE_exp = 0.3,
    beta_OStoOE_exp = 0.3,
    # 水平多效性
    prop_negative_pleiotropy =
        0, # 不满足inside假设
    mean_beta_FStoOE_out = 0,
    sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0,
    sd_beta_MStoOE_out = 0.05,
    mean_beta_OStoOE_out = 0,
    sd_beta_OStoOE_out = 0.05,
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
    beta_exp_to_out = 0,
    beta_confounding_exp = 0.2,
    beta_confounding_out = 0.2,
    correlation = 0.2,
    seed = NULL
)
small_simulation <- simulation_cmlfamily_multsnps_cpp(500, sim_params_list)

small_simulation_tibble <- simulation_ListtoTibble(small_simulation)
write.csv2(as.data.frame(small_simulation_tibble$tibble), file = "MR 模拟/MR 模拟结果/基础测试.csv")

# %% 大型模拟

beta_true_values <- c(-0.15, -0.10, -0.05, 0, 0.15, 0.10, 0.05)
indirectGeneticEffect_values <- c(0, 0.1)
sampleOverlapParameter_values <- c(0, 0.3, 0.7)
iNSIDEAssumption_index_values <- c(0, 0.5) # Represents prop_negative_pleiotropy
Violations_IVA <- c(0, 1, 3)

# 构建一个笛卡尔集合
param_grid <- expand.grid(
    current_Violations_IVA = Violations_IVA,
    current_prop_negative = iNSIDEAssumption_index_values,
    current_overlap = sampleOverlapParameter_values,
    current_indirect = indirectGeneticEffect_values,
    current_beta = beta_true_values,
    stringsAsFactors = FALSE # 通常建议，避免字符转因子
)

run_simulation_wrapper <- function(current_Violations_IVA, current_prop_negative,
                                   current_overlap, current_indirect, current_beta) {
    sim_params_list <- list(
        n = 10, num_pleiotropic = current_Violations_IVA, N_out = 4000,
        overlap_prop = current_overlap,
        beta_FStoOE_exp = 0.3, beta_MStoOE_exp = 0.3, beta_OStoOE_exp = 0.3,
        prop_negative_pleiotropy = current_prop_negative,
        mean_beta_FStoOE_out = current_indirect,
        sd_beta_FStoOE_out = 0.05,
        mean_beta_MStoOE_out = current_indirect,
        sd_beta_MStoOE_out = 0.05,
        mean_beta_OStoOE_out = current_indirect,
        sd_beta_OStoOE_out = 0.05,
        compatibility_selection_prop = 0,
        compatibility_selection_geno = "independent", correlation_param = 0.5,
        compatibility_selection_factor_exp = 0,
        compatibility_selection_factor_out = 0,
        crowd_stratification_differences = 0,
        beta_exp_to_out = current_beta,
        beta_confounding_exp = 0.2,
        beta_confounding_out = 0.2,
        correlation = 0.2,
        seed = NULL
    )
    simulation_cmlfamily_multsnps_cpp(1, sim_params_list)
}

# param_grid 的列名必须与 run_simulation_wrapper 函数的参数名匹配
results_list_pmap <- pmap(param_grid, run_simulation_wrapper)

# results_list_pmap 也会是一个包含所有模拟结果的列表
