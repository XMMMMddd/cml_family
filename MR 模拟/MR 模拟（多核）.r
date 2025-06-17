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
            triplet_family_simulation_once_robust(
                # --- 从列表解包参数 ---
                n = sim_params_list$n, num_pleiotropic = sim_params_list$num_pleiotropic,
                N_exp = sim_params_list$N_exp, N_out = sim_params_list$N_out,
                overlap_prop = sim_params_list$overlap_prop, p_f = sim_params_list$p_f, p_m = sim_params_list$p_m,
                beta_FStoOE_exp = sim_params_list$beta_FStoOE_exp,
                beta_MStoOE_exp = sim_params_list$beta_MStoOE_exp,
                beta_OStoOE_exp = sim_params_list$beta_OStoOE_exp,
                mean_beta_FStoOE_out = sim_params_list$mean_beta_FStoOE_out,
                sd_beta_FStoOE_out = sim_params_list$sd_beta_FStoOE_out,
                mean_beta_MStoOE_out = sim_params_list$mean_beta_MStoOE_out,
                sd_beta_MStoOE_out = sim_params_list$sd_beta_MStoOE_out,
                mean_beta_OStoOE_out = sim_params_list$mean_beta_OStoOE_out,
                sd_beta_OStoOE_out = sim_params_list$sd_beta_OStoOE_out,
                prop_negative_pleiotropy = sim_params_list$prop_negative_pleiotropy,
                assortative_mating_prob = sim_params_list$assortative_mating_prob,
                assortative_mating_strength = sim_params_list$assortative_mating_strength,
                crowd_stratification_differences = sim_params_list$crowd_stratification_differences,
                beta_exp_to_out = sim_params_list$beta_exp_to_out,
                beta_confounding_exp = sim_params_list$beta_confounding_exp,
                beta_confounding_out = sim_params_list$beta_confounding_out,
                correlation = sim_params_list$correlation, seed = NULL
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
    n = 10,
    num_pleiotropic = 0,
    N_exp = 1000, N_out = 1000,
    overlap_prop = 0, p_f = 0.3, p_m = 0.3,
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3, mean_beta_FStoOE_out = 0.1,
    sd_beta_FStoOE_out = 0.05, mean_beta_MStoOE_out = 0.1,
    sd_beta_MStoOE_out = 0.05, mean_beta_OStoOE_out = 0.1,
    sd_beta_OStoOE_out = 0.05, prop_negative_pleiotropy = 0.5,
    assortative_mating_prob = 0, assortative_mating_strength = 0,
    crowd_stratification_differences = 0, beta_exp_to_out = 0,
    beta_confounding_exp = 0.2, beta_confounding_out = 0.2,
    correlation = 0.2, seed = NULL
)
small_simulation <- simulation_cmlfamily_multsnps_cpp(50, sim_params_list)

small_simulation_tibble <- simulation_ListtoTibble(small_simulation)

# write.csv(small_simulation_tibble, "important_test.csv")
# sum(small_simulation_tibble$tibble$c_p_value < 0.05) / 500
# sum(small_simulation_tibble$tibble$b_p_value < 0.05, na.rm = TRUE) /
#     (500 - sum(is.na(small_simulation_tibble$tibble$b_p_value)))

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

# %% 多核模拟
# 多核并行版本的 triplet_family_simulation
triplet_family_simulation_parallel <- function(
    # 仿真参数
    n_simulations = 100,
    n_cores = NULL, # 核心数，NULL为自动检测
    progress = TRUE, # 是否显示进度条
    save_results = TRUE, # 是否保存详细结果
    output_file = NULL, # 输出文件路径
    seed_base = 12345, # 基础种子，每个仿真会有不同的种子
    dependency_file = "MR 模拟/MR 模拟（单次）.R", # 依赖文件路径
    # 原函数的所有参数
    n = 10,
    num_pleiotropic = 0,
    N_exp = 1000, N_out = 1000,
    overlap_prop = 0, p_f = 0.3, p_m = 0.3,
    beta_FStoOE_exp = 0.3, beta_MStoOE_exp = 0.3,
    beta_OStoOE_exp = 0.3, mean_beta_FStoOE_out = 0.1,
    sd_beta_FStoOE_out = 0.05, mean_beta_MStoOE_out = 0.1,
    sd_beta_MStoOE_out = 0.05, mean_beta_OStoOE_out = 0.1,
    sd_beta_OStoOE_out = 0.05, prop_negative_pleiotropy = 0.5,
    assortative_mating_prob = 0, assortative_mating_strength = 0,
    crowd_stratification_differences = 0, beta_exp_to_out = 0,
    beta_confounding_exp = 0.2, beta_confounding_out = 0.2,
    correlation = 0.2) {
    # ===== 加载必需的包 =====
    required_packages <- c("parallel", "foreach", "doParallel")

    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("请安装必需的包:", pkg, "\n可以使用: install.packages('", pkg, "')", sep = ""))
        }
    }

    # ===== 参数验证 =====
    if (!is.numeric(n_simulations) || n_simulations <= 0 || n_simulations != floor(n_simulations)) {
        stop("n_simulations 必须是正整数")
    }

    if (!is.null(n_cores)) {
        if (!is.numeric(n_cores) || n_cores <= 0 || n_cores != floor(n_cores)) {
            stop("n_cores 必须是正整数或NULL")
        }
    }

    # ===== 设置并行环境 =====
    # 检测可用核心数
    max_cores <- parallel::detectCores()
    if (is.null(n_cores)) {
        n_cores <- max(1, max_cores - 1) # 保留一个核心给系统
    } else {
        n_cores <- min(n_cores, max_cores)
    }

    cat(paste("检测到", max_cores, "个CPU核心，将使用", n_cores, "个核心进行并行计算\n"))

    # 检查是否有 triplet_family_simulation_once_robust 函数
    if (!exists("triplet_family_simulation_once_robust", mode = "function")) {
        stop("未找到 triplet_family_simulation_once_robust 函数，请先加载该函数")
    }

    # ===== 准备仿真参数 =====
    # 为每个仿真生成不同的种子
    seeds <- if (is.null(seed_base)) {
        NULL
    } else {
        seed_base + 1:n_simulations
    }

    # 创建参数列表
    simulation_params <- list(
        n = n, num_pleiotropic = num_pleiotropic,
        N_exp = N_exp, N_out = N_out,
        overlap_prop = overlap_prop, p_f = p_f, p_m = p_m,
        beta_FStoOE_exp = beta_FStoOE_exp, beta_MStoOE_exp = beta_MStoOE_exp,
        beta_OStoOE_exp = beta_OStoOE_exp, mean_beta_FStoOE_out = mean_beta_FStoOE_out,
        sd_beta_FStoOE_out = sd_beta_FStoOE_out, mean_beta_MStoOE_out = mean_beta_MStoOE_out,
        sd_beta_MStoOE_out = sd_beta_MStoOE_out, mean_beta_OStoOE_out = mean_beta_OStoOE_out,
        sd_beta_OStoOE_out = sd_beta_OStoOE_out, prop_negative_pleiotropy = prop_negative_pleiotropy,
        assortative_mating_prob = assortative_mating_prob, assortative_mating_strength = assortative_mating_strength,
        crowd_stratification_differences = crowd_stratification_differences, beta_exp_to_out = beta_exp_to_out,
        beta_confounding_exp = beta_confounding_exp, beta_confounding_out = beta_confounding_out,
        correlation = correlation
    )

    # ===== 定义单次仿真包装函数 =====
    run_single_simulation <- function(sim_id, params, seed_val = NULL) {
        tryCatch(
            {
                # 设置种子
                if (!is.null(seed_val)) {
                    params$seed <- seed_val
                }

                # 执行仿真 - 直接调用函数而不是do.call
                result <- triplet_family_simulation_once_robust(
                    n = params$n,
                    num_pleiotropic = params$num_pleiotropic,
                    N_exp = params$N_exp,
                    N_out = params$N_out,
                    overlap_prop = params$overlap_prop,
                    p_f = params$p_f,
                    p_m = params$p_m,
                    beta_FStoOE_exp = params$beta_FStoOE_exp,
                    beta_MStoOE_exp = params$beta_MStoOE_exp,
                    beta_OStoOE_exp = params$beta_OStoOE_exp,
                    mean_beta_FStoOE_out = params$mean_beta_FStoOE_out,
                    sd_beta_FStoOE_out = params$sd_beta_FStoOE_out,
                    mean_beta_MStoOE_out = params$mean_beta_MStoOE_out,
                    sd_beta_MStoOE_out = params$sd_beta_MStoOE_out,
                    mean_beta_OStoOE_out = params$mean_beta_OStoOE_out,
                    sd_beta_OStoOE_out = params$sd_beta_OStoOE_out,
                    prop_negative_pleiotropy = params$prop_negative_pleiotropy,
                    assortative_mating_prob = params$assortative_mating_prob,
                    assortative_mating_strength = params$assortative_mating_strength,
                    crowd_stratification_differences = params$crowd_stratification_differences,
                    beta_exp_to_out = params$beta_exp_to_out,
                    beta_confounding_exp = params$beta_confounding_exp,
                    beta_confounding_out = params$beta_confounding_out,
                    correlation = params$correlation,
                    seed = params$seed
                )

                # 添加仿真ID
                result$simulation_id <- sim_id
                result$seed_used <- seed_val
                result$timestamp <- Sys.time()

                return(result)
            },
            error = function(e) {
                warning(paste("仿真", sim_id, "失败:", e$message))

                # 返回失败结果
                failed_result <- list(
                    simulation_id = sim_id,
                    seed_used = seed_val,
                    timestamp = Sys.time(),
                    error = e$message,
                    # 初始化所有预期的字段为NA
                    a_theta_point = NA, a_theta_se = NA, a_z = NA, a_p_value = NA, a_duration = NA,
                    b_theta_point = NA, b_theta_se = NA, b_z = NA, b_p_value = NA, b_duration = NA,
                    c_theta_point = NA, c_theta_se = NA, c_z = NA, c_p_value = NA, c_duration = NA,
                    d_theta_point = NA, d_theta_se = NA, d_z = NA, d_p_value = NA, d_duration = NA,
                    e_theta_point = NA, e_theta_se = NA, e_z = NA, e_p_value = NA, e_duration = NA,
                    f_theta_point = NA, f_theta_se = NA, f_z = NA, f_p_value = NA, f_duration = NA,
                    g_theta_point = NA, g_theta_se = NA, g_z = NA, g_p_value = NA, g_duration = NA,
                    h_theta_point = NA, h_theta_se = NA, h_z = NA, h_p_value = NA, h_duration = NA,
                    i_theta_point = NA, i_theta_se = NA, i_z = NA, i_p_value = NA, i_duration = NA
                )

                return(failed_result)
            }
        )
    }

    # ===== 执行并行仿真 =====
    cat(paste("开始执行", n_simulations, "次仿真...\n"))
    start_time <- Sys.time()

    # 设置并行集群
    if (n_cores > 1) {
        cl <- parallel::makeCluster(n_cores)
        doParallel::registerDoParallel(cl)

        # 导出dependency_file变量到集群
        parallel::clusterExport(cl, "dependency_file", envir = environment())

        # 在每个核心上加载依赖和设置环境
        parallel::clusterEvalQ(cl, {
            # 加载依赖文件
            source(dependency_file)

            # 确保必要的包被加载
            suppressWarnings({
                library(methods) # 确保methods包可用
                if (exists("library")) {
                    try(library(MendelianRandomization), silent = TRUE)
                    try(library(MASS), silent = TRUE)
                }
            })
        })

        # 检查关键全局变量是否存在，如果存在则导出到集群
        global_vars_to_export <- c("n_independent", "p_trio", "n_snps", "overlap_cor")
        existing_vars <- c()

        for (var in global_vars_to_export) {
            if (exists(var, envir = .GlobalEnv)) {
                existing_vars <- c(existing_vars, var)
            }
        }

        if (length(existing_vars) > 0) {
            cat(paste("导出全局变量到集群:", paste(existing_vars, collapse = ", "), "\n"))
            parallel::clusterExport(cl, existing_vars, envir = .GlobalEnv)
        }

        # 验证函数是否加载成功
        cluster_check <- parallel::clusterEvalQ(cl, {
            exists("triplet_family_simulation_once_robust", mode = "function")
        })

        if (!all(unlist(cluster_check))) {
            parallel::stopCluster(cl)
            stop("部分集群节点未能成功加载依赖文件，请检查文件路径是否正确")
        }

        cat("集群设置完成，所有节点已成功加载依赖\n")

        tryCatch({
            if (progress && requireNamespace("pbapply", quietly = TRUE)) {
                # 使用 pbapply 显示进度条
                results_list <- pbapply::pblapply(1:n_simulations, function(i) {
                    run_single_simulation(i, simulation_params, if (is.null(seeds)) NULL else seeds[i])
                }, cl = cl)
            } else {
                # 使用标准 parLapply
                results_list <- parallel::parLapply(cl, 1:n_simulations, function(i) {
                    run_single_simulation(i, simulation_params, if (is.null(seeds)) NULL else seeds[i])
                })
            }
        }, finally = {
            parallel::stopCluster(cl)
        })
    } else {
        # 单核执行
        if (progress && requireNamespace("pbapply", quietly = TRUE)) {
            results_list <- pbapply::pblapply(1:n_simulations, function(i) {
                run_single_simulation(i, simulation_params, if (is.null(seeds)) NULL else seeds[i])
            })
        } else {
            results_list <- lapply(1:n_simulations, function(i) {
                run_single_simulation(i, simulation_params, if (is.null(seeds)) NULL else seeds[i])
                if (progress && i %% max(1, n_simulations %/% 20) == 0) {
                    cat(paste("完成进度:", round(i / n_simulations * 100, 1), "%\n"))
                }
            })
        }
    }

    end_time <- Sys.time()
    total_duration <- as.numeric(end_time - start_time)

    # ===== 处理结果 =====
    cat("处理结果...\n")

    # 统计成功和失败的仿真
    successful_sims <- sum(sapply(results_list, function(x) is.null(x$error)))
    failed_sims <- n_simulations - successful_sims

    cat(paste("仿真完成! 总用时:", round(total_duration, 2), "秒\n"))
    cat(paste("成功:", successful_sims, "次, 失败:", failed_sims, "次\n"))
    cat(paste("成功率:", round(successful_sims / n_simulations * 100, 1), "%\n"))

    # 转换为数据框
    if (save_results) {
        # 提取所有字段名
        all_fields <- unique(unlist(lapply(results_list, names)))

        # 创建数据框
        results_df <- data.frame(
            simulation_id = sapply(results_list, function(x) x$simulation_id %||% NA),
            stringsAsFactors = FALSE
        )

        # 添加所有其他字段
        for (field in setdiff(all_fields, "simulation_id")) {
            results_df[[field]] <- sapply(results_list, function(x) x[[field]] %||% NA)
        }
    } else {
        results_df <- NULL
    }

    # ===== 汇总统计 =====
   # ===== 汇总统计 =====
method_names <- c("a", "b", "c", "d", "e", "f", "g", "h", "i")
summary_stats <- list()

# 使用仿真参数中的真实因果效应值
true_theta <- beta_exp_to_out

for (method in method_names) {
    theta_values <- sapply(results_list, function(x) x[[paste0(method, "_theta_point")]])
    se_values <- sapply(results_list, function(x) x[[paste0(method, "_theta_se")]])
    p_values <- sapply(results_list, function(x) x[[paste0(method, "_p_value")]])
    
    # 移除NA值
    valid_indices <- !is.na(theta_values) & !is.na(se_values)
    valid_theta <- theta_values[valid_indices]
    valid_se <- se_values[valid_indices]
    valid_p <- p_values[!is.na(p_values)]
    
    if (length(valid_theta) > 0) {
        # 计算95%置信区间的覆盖率
        ci_lower_95 <- valid_theta - 1.96 * valid_se
        ci_upper_95 <- valid_theta + 1.96 * valid_se
        coverage_95 <- mean(ci_lower_95 <= true_theta & true_theta <= ci_upper_95)
        
        # 计算99%置信区间的覆盖率
        ci_lower_99 <- valid_theta - 2.576 * valid_se
        ci_upper_99 <- valid_theta + 2.576 * valid_se
        coverage_99 <- mean(ci_lower_99 <= true_theta & true_theta <= ci_upper_99)
        
        summary_stats[[method]] <- list(
            n_valid = length(valid_theta),
            mean_theta = mean(valid_theta),
            sd_theta = sd(valid_theta),
            median_theta = median(valid_theta),
            power_005 = if (length(valid_p) > 0) mean(valid_p < 0.05) else NA,
            power_001 = if (length(valid_p) > 0) mean(valid_p < 0.01) else NA,
            coverage_95 = coverage_95,
            coverage_99 = coverage_99,
            # 额外统计信息
            mean_se = mean(valid_se),
            bias = mean(valid_theta) - true_theta,
            mse = mean((valid_theta - true_theta)^2),  # 均方误差
            rmse = sqrt(mean((valid_theta - true_theta)^2))  # 均方根误差
        )
    } else {
        summary_stats[[method]] <- list(
            n_valid = 0,
            mean_theta = NA,
            sd_theta = NA,
            median_theta = NA,
            power_005 = NA,
            power_001 = NA,
            coverage_95 = NA,
            coverage_99 = NA,
            mean_se = NA,
            bias = NA,
            mse = NA,
            rmse = NA
        )
    }
}

    # ===== 保存结果 =====
    if (!is.null(output_file) && save_results) {
        tryCatch(
            {
                if (grepl("\\.rds$", output_file)) {
                    saveRDS(list(results = results_df, summary = summary_stats), output_file)
                } else if (grepl("\\.csv$", output_file)) {
                    write.csv(results_df, output_file, row.names = FALSE)
                } else {
                    save(results_df, summary_stats, file = output_file)
                }
                cat(paste("结果已保存到:", output_file, "\n"))
            },
            error = function(e) {
                warning(paste("保存文件失败:", e$message))
            }
        )
    }

    # ===== 返回结果 =====
    final_result <- list(
        summary_statistics = summary_stats,
        simulation_info = list(
            n_simulations = n_simulations,
            successful_simulations = successful_sims,
            failed_simulations = failed_sims,
            success_rate = successful_sims / n_simulations,
            total_duration_seconds = total_duration,
            n_cores_used = n_cores,
            parameters_used = simulation_params
        )
    )

    if (save_results) {
        final_result$detailed_results <- results_df
        final_result$raw_results <- results_list
    }

    return(final_result)
}


# 用于处理NULL值的操作符
`%||%` <- function(x, y) if (is.null(x)) y else x


# 提供一个简化的接口
run_mr_simulation <- function(n_simulations = 100, n_cores = NULL, ...) {
    triplet_family_simulation_parallel(
        n_simulations = n_simulations,
        n_cores = n_cores,
        ...
    )
}

# ===== 使用示例 =====
if (FALSE) { # 防止意外执行
    # 基本用法
    results <- triplet_family_simulation_parallel(
        n_simulations = 1000,
        n_cores = 10,
        n = 10,
        output_file = NULL, # 输出文件路径
        num_pleiotropic = 0,
        N_exp = 3000,
        N_out = 3000,
        overlap_prop = 0,
        p_f = 0.3, p_m = 0.3,
        beta_FStoOE_exp = 0.3,
        beta_MStoOE_exp = 0.3,
        beta_OStoOE_exp = 0.,
        mean_beta_FStoOE_out = 0.1,
        sd_beta_FStoOE_out = 0.05,
        mean_beta_MStoOE_out = 0.1,
        sd_beta_MStoOE_out = 0.05,
        mean_beta_OStoOE_out = 0.1,
        sd_beta_OStoOE_out = 0.05,
        prop_negative_pleiotropy = 0.5,
        assortative_mating_prob = 0,
        assortative_mating_strength = 0,
        crowd_stratification_differences = 0,
        beta_exp_to_out = 0.02,
        beta_confounding_exp = 0.2,
        beta_confounding_out = 0.2,
        correlation = 0.2, seed = NULL
    )

    # 查看汇总结果
    print(results$summary_statistics)

    # 快速启动
    results2 <- run_mr_simulation(
        n_simulations = 50,
        n = 15,
        beta_exp_to_out = 0.1
    )
}

results$detailed_results
