# %%
# install.packages("doParallel")
library(doParallel)
library(purrr)
# install.packages("here")
library(here)
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
# sim_params_list <- list( # --- 从列表解包参数 ---
#     n = 10, num_pleiotropic = 0, N_out = 1000,
#     # 样本重叠
#     overlap_prop = 0,
#     # 直接效应
#     beta_FStoOE_exp = 0.3,
#     beta_MStoOE_exp = 0.3,
#     beta_OStoOE_exp = 0.3,
#     # 水平多效性
#     prop_negative_pleiotropy =
#         0, # 不满足inside假设
#     mean_beta_FStoOE_out = 0,
#     sd_beta_FStoOE_out = 0.05,
#     mean_beta_MStoOE_out = 0,
#     sd_beta_MStoOE_out = 0.05,
#     mean_beta_OStoOE_out = 0,
#     sd_beta_OStoOE_out = 0.05,
#     # 选型婚配
#     ## 选型婚配基因型
#     compatibility_selection_prop = 0,
#     compatibility_selection_geno = "independent", correlation_param = 0.5,
#     ## 选型婚配暴露相关
#     compatibility_selection_factor_exp = 0,
#     compatibility_selection_factor_out = 0,
#     # 人群分层（双人群差异）
#     crowd_stratification_differences = 0,
#     # 其他参数设置
#     beta_exp_to_out = 0,
#     beta_confounding_exp = 0.2,
#     beta_confounding_out = 0.2,
#     correlation = 0.2,
#     seed = NULL
# )
# small_simulation <- simulation_cmlfamily_multsnps_cpp(500, sim_params_list)

# small_simulation_tibble <- simulation_ListtoTibble(small_simulation)
# write.csv2(as.data.frame(small_simulation_tibble$tibble), file = "MR 模拟/MR 模拟结果/基础测试.csv")

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
    simulation_cmlfamily_multsnps_cpp(500, sim_params_list)
}
# 2. 定义输出目录的名称
output_directory_name <- "MR 模拟/MR 模拟结果"

# 3. 创建输出目录 (如果尚不存在)
# 使用 here() 来确保路径相对于项目根目录是正确的
if (!dir.exists(here(output_directory_name))) {
    dir.create(here(output_directory_name), recursive = TRUE)
    message(paste("已创建输出目录:", here(output_directory_name)))
} else {
    message(paste("输出目录已存在:", here(output_directory_name)))
}

# 4. 定义一个新的包装函数，用于运行模拟并保存结果
run_and_save_simulation_csv <- function(current_Violations_IVA, current_prop_negative,
                                        current_overlap, current_indirect, current_beta,
                                        output_dir) {
    # a. 根据参数生成唯一且描述性的文件名
    # 为了文件名安全，替换特殊字符：将点(.)替换为'p'，负号(-)替换为'neg'
    val_iva_str <- as.character(current_Violations_IVA)
    val_prop_neg_str <- gsub("\\.", "p", as.character(current_prop_negative))
    val_overlap_str <- gsub("\\.", "p", as.character(current_overlap))
    val_indirect_str <- gsub("\\.", "p", as.character(current_indirect))
    val_beta_str <- gsub("-", "neg", gsub("\\.", "p", as.character(current_beta)))

    filename <- paste0(
        "sim_IVA-", val_iva_str,
        "_propNeg-", val_prop_neg_str,
        "_overlap-", val_overlap_str,
        "_indirect-", val_indirect_str,
        "_beta-", val_beta_str,
        ".csv"
    )

    filepath <- file.path(output_dir, filename)

    # b. 运行原始的模拟包装函数
    message(paste("正在运行参数组合:", filename)) # 打印当前运行的参数组合信息
    simulation_result <- tryCatch(
        {
            run_simulation_wrapper(
                current_Violations_IVA,
                current_prop_negative,
                current_overlap,
                current_indirect,
                current_beta
            )
        },
        error = function(e) {
            warning(paste("模拟运行失败:", filename, "错误信息:", e$message))
            return(NULL) # 如果模拟失败，返回NULL
        }
    )

    # c. 将结果保存到CSV文件
    if (!is.null(simulation_result)) {
        tryCatch(
            {
                write.csv(simulation_result, filepath, row.names = FALSE)
                message(paste("结果已保存到:", filepath))
                return(filepath) # 返回文件路径表示成功
            },
            error = function(e) {
                warning(paste("保存CSV文件失败:", filepath, "错误信息:", e$message))
                return(NULL) # 如果保存失败，返回NULL
            }
        )
    } else {
        return(NULL) # 如果模拟本身失败，也返回NULL
    }
}

# 5. 使用 pmap 应用新的包装函数
# pmap 会将 param_grid 的每一行作为参数传递给 run_and_save_simulation_csv 函数
# output_dir 参数会作为附加的固定参数传递给函数
message("开始批量运行模拟并保存结果...")
list_of_created_files <- pmap(
    .l = param_grid,
    .f = run_and_save_simulation_csv,
    output_dir = here(output_directory_name) # 将输出目录作为固定参数传递
)
message("所有模拟运行和保存操作已完成。")

# 6. （可选）检查生成的文件列表或处理失败的情况
# list_of_created_files 将包含成功保存的文件路径或NULL（如果失败）
# 您可以检查这个列表来确认哪些成功了，哪些失败了
successful_saves <- Filter(Negate(is.null), list_of_created_files)
failed_saves_count <- length(list_of_created_files) - length(successful_saves)

message(paste("成功保存的文件数量:", length(successful_saves)))
if (failed_saves_count > 0) {
    message(paste("失败的保存/模拟数量:", failed_saves_count))
}

# 您可以打印出前几个成功保存的文件名
# print(head(successful_saves))



