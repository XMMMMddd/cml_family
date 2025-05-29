# %%
# renv::deactivate()
# install.packages("doParallel")
library(doParallel)
library(purrr)
# install.packages("progress")
library(progress)
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
            triplet_family_simulation_once_robust_mr(
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
                assortative_mating_prob =
                    sim_params_list$assortative_mating_prob, # 选型婚配比例
                assortative_mating_strength =
                    sim_params_list$assortative_mating_strength, # 选型婚配强度
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
    assortative_mating_prob = 0.999, # 选型婚配比例
    assortative_mating_strength = 0.1, # 选型婚配强度
    # 人群分层（双人群差异）
    crowd_stratification_differences = 0,
    # 其他参数设置
    beta_exp_to_out = -0.15,
    beta_confounding_exp = 0.2,
    beta_confounding_out = 0.2,
    correlation = 0.2,
    seed = NULL
)
small_simulation <- simulation_cmlfamily_multsnps_cpp(10, sim_params_list)

small_simulation_tibble <- simulation_ListtoTibble(small_simulation)
write.csv2(as.data.frame(small_simulation_tibble$tibble),
    file = "MR 模拟/MR 模拟结果/cml_family_ver2的模拟结果/基础测试.csv"
)

# %% 大型模拟

run_simulations_with_fancy_pb <- function(
    n_simulations_per_combination = 10, # pb for progress bar
    csv_file_base_address, # CSV文件存储的基础路径
    change_parameters, # 包含参数变化的列表
    base_parameters # 包含基础参数的列表
    ) {
    # 确保 progress 包已加载
    if (!requireNamespace("progress", quietly = TRUE)) {
        stop("Package 'progress' is needed for this function to work. Please install it.", call. = FALSE)
    }

    # 1. 生成所有参数组合
    parameter_grid <- expand.grid(change_parameters, stringsAsFactors = FALSE)
    total_simulations <- nrow(parameter_grid)

    cat(paste("Total number of parameter combinations to simulate:", total_simulations, "\n"))

    # 初始化 progress 包的进度条
    # format 参数定义了进度条的样式
    # :bar 是进度条本身
    # :current 是当前迭代次数
    # :total 是总迭代次数
    # :percent 是完成百分比
    # :eta 是预计剩余时间
    # :elapsed 是已用时间
    pb <- progress::progress_bar$new(
        format = "  Simulating [:bar] :current/:total (:percent) | ETA: :eta | Elapsed: :elapsed",
        total = total_simulations,
        width = 80, # 进度条的宽度
        clear = FALSE # 完成后不清除进度条
    )

    # 2. 循环遍历每一种参数组合
    for (i in 1:total_simulations) {
        # 在每次迭代开始时，更新进度条的状态 (如果你想在参数名中显示当前参数)
        # 或者，更简单的方式是只在迭代结束前调用 pb$tick()
        # 为了简单起见，我们主要使用 pb$tick()

        current_sim_params <- base_parameters
        current_file_name_parts <- c()

        for (param_name in names(parameter_grid)) {
            current_value <- parameter_grid[i, param_name]
            current_sim_params[[param_name]] <- current_value
            current_file_name_parts <- c(current_file_name_parts, paste0(param_name, "_", current_value))
        }

        specific_file_name <- paste(current_file_name_parts, collapse = "__")
        specific_file_name <- gsub("\\.", "p", specific_file_name)
        specific_file_name <- gsub("-", "neg", specific_file_name)
        specific_file_name <- paste0(specific_file_name, ".csv")
        full_file_path <- file.path(csv_file_base_address, specific_file_name)

        # (可选) 更新进度条的 token，如果你想显示当前正在处理的参数
        # current_param_display <- paste(names(parameter_grid), "=", parameter_grid[i, ], collapse = ", ")
        # pb$tick(0, tokens = list(what = current_param_display)) # tick(0) 不推进进度，只更新token

        # 3. 运行模拟
        n_simulations_per_combination <- n_simulations_per_combination # 默认值
        if (!is.null(current_sim_params$n_simulations_per_combination)) {
            n_simulations_per_combination <- current_sim_params$n_simulations_per_combination
        } else if (!is.null(current_sim_params$n_simulations)) {
            n_simulations_per_combination <- current_sim_params$n_simulations
        }

        small_simulation <- simulation_cmlfamily_multsnps_cpp(
            n_simulations_per_combination,
            current_sim_params
        )
        small_simulation_tibble <- simulation_ListtoTibble(small_simulation)

        # 4. 保存结果到 CSV 文件
        if (!dir.exists(csv_file_base_address)) {
            dir.create(csv_file_base_address, recursive = TRUE)
        }

        if (!is.null(small_simulation_tibble) && !is.null(small_simulation_tibble$tibble)) {
            write.csv2(as.data.frame(small_simulation_tibble$tibble),
                file = full_file_path,
                row.names = FALSE
            )
        } else {
            warning_message <- paste0(
                "Warning: Simulation ", i, "/", total_simulations,
                " did not return valid data to save for parameters: ",
                paste(names(parameter_grid), parameter_grid[i, ], collapse = ", ")
            )
            warning(warning_message)
        }

        # 推进进度条 (在每次迭代的末尾)
        pb$tick()
    }

    # 进度条完成后，可以打印一条完成信息
    # pb$terminate() # 如果 clear = TRUE, 这会清除进度条。如果 clear = FALSE, 它只是确保光标在新行。
    cat("\nAll simulations completed.\n")
}
csv_file_address <- "MR 模拟/MR 模拟结果/cml_family_ver2的模拟结果"
change_parameters <- list(
    assortative_mating_strength = c(0.1, 1, 10, 100),
    beta_exp_to_out = c(-0.15, -0.10, -0.05, 0, 0.05, 0.10, 0.15),
    num_pleiotropic = c(0, 1, 3),
    crowd_stratification_differences = c(0, 0.05, 0.1)
)
base_parameters <- sim_params_list <- list( # --- 从列表解包参数 ---
    n = 10, num_pleiotropic = 0, N_out = 1000,
    # 样本重叠
    overlap_prop = 0,
    # 直接效应
    beta_FStoOE_exp = 0.1,
    beta_MStoOE_exp = 0.1,
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
    assortative_mating_prob = 0.999, # 选型婚配比例
    assortative_mating_strength = 1000, # 选型婚配强度
    # 人群分层（双人群差异）
    crowd_stratification_differences = 0,
    # 其他参数设置
    beta_exp_to_out = 0,
    beta_confounding_exp = 0.2,
    beta_confounding_out = 0.2,
    correlation = 0.2,
    seed = NULL
)

test <- run_simulations_with_fancy_pb(
    n_simulations_per_combination = 500,
    csv_file_address, change_parameters, base_parameters
)
