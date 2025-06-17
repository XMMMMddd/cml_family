# %%
# install.packages("RcppArmadillo")
# install.packages("Rcpp")
library(RcppArmadillo)
library(Rcpp)

# %%
# cML的参数估计 CPP 版本
sourceCpp("cml家庭函数/cml家庭函数一代/cML_para_rcpp.cpp")
# cMl的参数估计
cMl_para <- function(
    k_par = 0, k_x = 0, beta_y_hat, beta_x_hat, Sigma_inv_x,
    Sigma_inv_y, value_diff = 1e-9) {
    # 给出迭代初始值
    r_x <- rep(0, nrow(beta_y_hat))
    r_p <- rep(0, nrow(beta_y_hat))
    r <- matrix(c(r_x, r_p), ncol = 2)

    theta <- 0
    theta_new <- theta
    theta_old <- 1
    # 给beta_x一个初始值
    beta_x <- beta_x_hat
    beta_x_old <- beta_x_hat
    # 记录迭代次数
    n_times <- 0
    while ((abs(theta_old - theta_new) > value_diff |
        abs(sum(beta_x_old - beta_x)) > value_diff) & n_times < 1000) {
        n_times <- n_times + 1
        beta_x_old <- beta_x
        theta_old <- theta_new
        # cat("迭代r\n")
        # 迭代 r
        # 首先开始调整r_x的这一部分
        r_p <- -(c(theta) * beta_x[, 2]) + beta_y_hat[, 2]
        r_x <- rep(0, nrow(beta_y_hat))

        # 让par这一项不影响我们的更新
        A <- matrix(0, nrow = nrow(beta_y_hat))
        beta_x_x <- theta * beta_x[, 1] + r_x

        beta_x_par <- theta * beta_x[, 2] + r_p
        beta_x_cal <- matrix(c(beta_x_x, beta_x_par), ncol = 2)

        for (i in 1:nrow(beta_y_hat)) {
            A[i] <- (beta_y_hat[i, ] - beta_x_cal[i, ]) %*% Sigma_inv_y[(2 * i - 1):(2 * i), ] %*%
                (beta_y_hat[i, ] - beta_x_cal[i, ])
        }
        A <- data.frame(
            A = A, r_x = r_x, beta_x_x = beta_x[, 1], beta_y_hat_x = beta_y_hat[, 1],
            SNPs_code = 1:nrow(beta_y_hat)
        )
        sorted_A <- A[order(A$A, decreasing = TRUE), ]
        sorted_A <- sorted_A %>% mutate(r_x_new = -(theta * beta_x_x) + beta_y_hat_x)
        # 前k_x项是不好的要把他压缩到0

        if (k_x == 0) {
            sorted_A$r_p[1:k_x] <- 0
        } else {
            sorted_A$r_x[1:k_x] <- sorted_A$r_x_new[1:k_x]
        }
        sorted_A <- sorted_A[order(sorted_A$SNPs_code), ]
        r_x_new <- sorted_A$r_x

        # 首先开始调整r_p的这一部分
        r_p <- rep(0, nrow(beta_y_hat))
        r_x <- -(theta * beta_x[, 1]) + beta_y_hat[, 1]

        A <- matrix(0, nrow = nrow(beta_y_hat))
        beta_x_x <- theta * beta_x[, 1] + r_x
        beta_par_x <- theta * beta_x[, 2] + r_p
        beta_x_cal <- matrix(c(beta_x_x, beta_par_x), ncol = 2)

        for (i in 1:nrow(beta_y_hat)) {
            A[i] <- (beta_y_hat[i, ] - beta_x_cal[i, ]) %*% Sigma_inv_y[(2 * i - 1):(2 * i), ] %*%
                (beta_y_hat[i, ] - beta_x_cal[i, ])
        }
        A <- data.frame(
            A = A, r_p = r_p, beta_x_par = beta_x[, 2], beta_y_hat_par = beta_y_hat[, 2],
            SNPs_code = 1:nrow(beta_y_hat)
        )
        sorted_A <- A[order(A$A, decreasing = TRUE), ]
        sorted_A <- sorted_A %>% mutate(r_par_new = -(theta * beta_x_par) + beta_y_hat_par)
        # 前k_par项是不好的要把他压缩到0

        if (k_par == 0) {
            sorted_A$r_p[1:k_par] <- 0
        } else {
            sorted_A$r_p[1:k_par] <- sorted_A$r_par_new[1:k_par]
        }

        sorted_A <- sorted_A[order(sorted_A$SNPs_code), ]
        r_par_new <- sorted_A$r_p

        # 更新r
        r_new <- matrix(c(r_x_new, r_par_new), ncol = 2)
        r <- r_new
        # print(r)
        # cat("迭代beta\n")
        # 迭代beta
        # 根据公式更新beta
        beta_x_new <- beta_x
        for (i in 1:nrow(beta_y_hat)) {
            # 求分子
            fenzi <- beta_x_hat[i, ] %*% Sigma_inv_x[(2 * i - 1):(2 * i), ] + theta * beta_y_hat[i, ] %*%
                Sigma_inv_y[(2 * i - 1):(2 * i), ] -
                theta * r[1, ] %*% Sigma_inv_y[(2 * i - 1):(2 * i), ]
            # 求分母
            fenmu <- (Sigma_inv_x[(2 * i - 1):(2 * i), ] + theta^2 * Sigma_inv_y[(2 * i - 1):(2 * i), ])
            # 求逆
            fenmu <- solve(fenmu)
            beta_x_new[i, ] <- fenzi %*% fenmu
        }
        beta_x <- beta_x_new

        # 更新theta
        fenzi <- 0
        fenmu <- 0
        for (i in 1:nrow(beta_y_hat)) {
            fenzi <- (beta_y_hat[i, ] - r_new[i, ]) %*% Sigma_inv_y[(2 * i - 1):(2 * i), ] %*%
                beta_x[i, ] + fenzi
            fenmu <- beta_x[i, ] %*% Sigma_inv_y[(2 * i - 1):(2 * i), ] %*% beta_x[i, ] + fenmu
        }
        theta_new <- fenzi / fenmu
        theta <- c(theta_new)
    }
    if (n_times == 1000) {
        print("迭代失败")
    }
    c <- list(theta = theta, beta_x = beta_x, r = r)


    return(c)
}


# %%
# 完整版的函数
cml_family_2_cn_cpp <- function(beta_x_hat, beta_y_hat, Sigma_inv_x, Sigma_inv_y,
                                multi_start_args = list(n_starts = 5, seed = NULL, min_freq_prop = 0.7, digits_compare = 6),
                                value_diff_optim = 1e-6,
                                verbose = FALSE, # 新增 verbose 参数，默认为 TRUE
                                n_eff_override = 1000 # 新增 n_eff 覆盖参数
) {
    n_eff_for_bic <- n_eff_override # 优先使用覆盖值
    n_datasets <- nrow(beta_x_hat) # 数据集的数量

    # --- 3. 设置网格和结果存储 ---
    k_grid <- expand.grid(k_par = 0:n_datasets, k_x = 0:n_datasets) # k 从 0 开始
    n_combinations <- nrow(k_grid)
    results_list <- vector("list", n_combinations) # 用于存储每次迭代的详细结果


    # --- 4. 顺序网格搜索 ---
    for (i in 1:n_combinations) {
        k_x_current <- k_grid$k_x[i]
        k_par_current <- k_grid$k_par[i]


        # 调用多起点优化函数
        optim_result <- cMl_para_multi_start_rcpp(
            n_starts = 5, sd_theta_start = 0.001, sd_beta_x_start = 0.001, digits_compare = 6, min_freq_prop = 0.5, seed = NULL, verbose = FALSE, k_par = k_par_current, k_x = k_x_current, beta_y_hat = beta_y_hat, beta_x_hat = beta_x_hat, Sigma_inv_x = Sigma_inv_x, Sigma_inv_y = Sigma_inv_y
        )

        current_theta <- NA_real_
        current_theta_se <- NA_real_
        optim_converged <- FALSE # 增加一个优化是否收敛的本地变量
        optim_error <- TRUE # 默认优化有错误
        bic_value <- Inf # 默认为 Inf
        stable_convergence <- FALSE # 默认不稳定


        if (optim_result$converged) {
            bic_value <- BIC_function_optimized(
                theta = optim_result$theta,
                k_x = k_x_current,
                k_par = k_par_current,
                beta_x = optim_result$beta_x,
                beta_x_hat = beta_x_hat,
                Sigma_inv_x = Sigma_inv_x,
                beta_y_hat = beta_y_hat,
                Sigma_inv_y = Sigma_inv_y,
                n_eff = n_eff_for_bic,
                r = optim_result$r
            )
        }


        # --- 存储结果 ---
        current_theta <- NA_real_
        current_theta_se <- NA_real_
        if (!is.null(optim_result) && optim_result$converged) {
            if (!is.null(optim_result$theta)) {
                current_theta <- as.numeric(optim_result$theta[1])
            }
            if (!is.null(optim_result$theta_se) && optim_result$converged) {
                current_theta_se <- as.numeric(optim_result$theta_se)
            }
        }
        results_list[[i]] <- list(
            k_x = k_x_current,
            k_par = k_par_current,
            theta = current_theta,
            theta_se = current_theta_se,
            bic = as.numeric(bic_value),
            converged = isTRUE(optim_result$converged),
            stable = stable_convergence,
            error_in_optim = isTRUE(optim_result$error),
            full_result = optim_result
        )
    }
    # 结束 for 循环

    if (verbose) cat("网格搜索完成。\n")

    # --- 5. 处理和选择结果 ---
    bic_results_df <- data.frame(
        k_x = sapply(results_list, function(x) x$k_x),
        k_par = sapply(results_list, function(x) x$k_par),
        theta = sapply(results_list, function(x) x$theta),
        theta_se = sapply(results_list, function(x) x$theta_se),
        bic = sapply(results_list, function(x) x$bic),
        converged = sapply(results_list, function(x) x$converged),
        stable = sapply(results_list, function(x) x$stable),
        error_in_optim = sapply(results_list, function(x) x$error_in_optim)
    )

    # 按 BIC 排序，将 Inf 排在后面
    bic_results_df_sorted <- bic_results_df[order(bic_results_df$bic), ]

    # 找到最小的有限 BIC 值对应的行
    best_row <- bic_results_df_sorted[is.finite(bic_results_df_sorted$bic), ][1, , drop = FALSE] # 使用 drop=FALSE 保持数据框结构

    best_full_result <- NULL # 初始化

    if (nrow(best_row) == 0 || !is.finite(best_row$bic)) {
        # 检查 best_row 是否存在且 BIC 有限
        if (verbose) {
            warning("所有 (k_x, k_par) 组合的 BIC 值均为 Inf 或计算/优化失败。无法找到最优模型。")
        }
        # 即使找不到最优，仍然计算加权平均（可能基于 NA 或 Inf，结果可能也是 NA）
    } else {
        # 获取最优的 k_x 和 k_par
        optimal_k_x <- best_row$k_x
        optimal_k_par <- best_row$k_par
        best_bic_value <- best_row$bic

        if (verbose) cat(paste0("找到最优组合：k_x = ", optimal_k_x, ", k_par = ", optimal_k_par, " (BIC = ", round(best_bic_value, 3), ")\n"))

        # 从 results_list 中查找对应的完整结果
        best_result_index <- which(sapply(results_list, function(x) x$k_x == optimal_k_x && x$k_par == optimal_k_par))[1]

        if (length(best_result_index) > 0 && !is.null(results_list[[best_result_index]]$full_result)) {
            best_full_result <- results_list[[best_result_index]]$full_result
        } else {
            if (verbose) warning("找到了最优的 k_x/k_par，但在原始结果列表中未能找到对应的完整优化结果。")
        }
    }

    # --- 计算 BIC 权重和加权平均 ---
    # 过滤掉 BIC 为 Inf 或 NA 的行，以及 theta 或 theta_se 为 NA 的行，以进行稳健的加权计算
    bic_results_valid <- bic_results_df_sorted[is.finite(bic_results_df_sorted$bic) &
        !is.na(bic_results_df_sorted$theta) &
        !is.na(bic_results_df_sorted$theta_se), ]

    theta_weight_avg <- NA_real_
    theta_se_weight_avg <- NA_real_

    if (nrow(bic_results_valid) > 0) {
        # 找到有效的最小 BIC
        min_valid_bic <- min(bic_results_valid$bic)
        # 计算权重，对 BIC 值进行移位以避免 exp(-Inf/2) 的数值问题
        bic_results_valid <- tryCatch(
            {
                dplyr::mutate(bic_results_valid, weight = exp(-(bic - min_valid_bic) / 2))
            },
            error = function(e) {
                if (verbose) warning("计算 BIC 权重时出错（可能因为极大的 BIC 值）。")
                return(dplyr::mutate(bic_results_valid, weight = NA_real_)) # 出错则权重设为 NA
            }
        )

        # 检查权重是否有效
        valid_weights <- !is.na(bic_results_valid$weight) & bic_results_valid$weight > 0
        if (any(valid_weights)) {
            sum_weights <- sum(bic_results_valid$weight[valid_weights], na.rm = TRUE)
            if (sum_weights > 0) {
                # 计算加权平均值
                theta_weight_avg <- sum(bic_results_valid$theta[valid_weights] * bic_results_valid$weight[valid_weights], na.rm = TRUE) / sum_weights
                # 计算加权平均标准误（这是一个简化，更复杂的方法可能需要考虑模型间的不确定性）
                # 这里简单地对 SE 进行加权平均，可能不是最佳统计实践，但提供一个粗略估计
                theta_se_weight_avg <- sum(bic_results_valid$theta_se[valid_weights] * bic_results_valid$weight[valid_weights], na.rm = TRUE) / sum_weights
            } else {
                if (verbose) warning("所有有效的 BIC 权重之和为零或负数，无法计算加权平均。")
            }
        } else {
            if (verbose) warning("没有计算出有效的 BIC 权重，无法计算加权平均。")
        }
    } else {
        if (verbose) warning("没有具有有效 BIC、theta 和 theta_se 的结果，无法计算加权平均。")
    }


    # --- 6. 返回结果 ---
    final_output <- list(
        theta_weight = theta_weight_avg, # 使用计算出的加权平均值
        theta_se_weight = theta_se_weight_avg, # 使用计算出的加权平均 SE
        best_result = best_full_result, # 最优模型的完整结果
        bic_results = bic_results_df_sorted # 所有组合的 BIC 结果（排序后）
    )

    return(final_output)
}

# 多起点函数
cMl_para_multi_start_rcpp <- function(n_starts = 10,
                                      sd_theta_start = 0.00001,
                                      sd_beta_x_start = 0.00001,
                                      digits_compare = 6,
                                      min_freq_prop = 0.5,
                                      seed = NULL,
                                      verbose = TRUE,
                                      ...) {
    # --- 1. 设置 ---
    if (!is.null(seed)) {
        set.seed(seed) # 设置随机种子以保证可复现性
    }

    # 捕获通过 '...' 传递给 cMl_para_rcpp 的参数
    rcpp_args <- list(...)

    # 确保必要的参数已提供
    required_args <- c("beta_y_hat", "beta_x_hat", "Sigma_inv_x", "Sigma_inv_y")
    missing_args <- setdiff(required_args, names(rcpp_args))
    if (length(missing_args) > 0) {
        stop("缺少 cMl_para_rcpp 所需的参数: ", paste(missing_args, collapse = ", "))
    }

    # 对 beta_x_hat 维度进行基本检查，用于生成噪声
    if (!is.matrix(rcpp_args$beta_x_hat) || length(dim(rcpp_args$beta_x_hat)) != 2) {
        stop("'beta_x_hat' 必须是一个矩阵。")
    }
    beta_x_hat_dim <- dim(rcpp_args$beta_x_hat)
    P <- beta_x_hat_dim[1] # 工具变量数量
    num_cols_beta_x <- beta_x_hat_dim[2] # beta_x 的列数（应为 2）

    all_results <- vector("list", n_starts) # 预分配存储结果的列表

    if (verbose) {
        cat(paste("开始使用 Rcpp 后端进行", n_starts, "次优化运行...\n"))
    }

    # --- 2. 多起点循环 ---
    for (i in 1:n_starts) {
        current_rcpp_args <- rcpp_args # 复制本次运行的基础参数

        run_desc <- "" # 运行描述
        if (i == 1) {
            # 第一次运行：通过传递 NULL 来使用 C++ 的默认初始值
            current_rcpp_args$theta_init_nullable <- NULL
            current_rcpp_args$beta_x_init_nullable <- NULL
            run_desc <- "(默认起点)"
        } else {
            # 后续运行：生成随机起始点
            # theta 起始值：从正态分布 N(0, sd^2) 中抽取
            current_rcpp_args$theta_init_nullable <- rnorm(1, mean = 0, sd = sd_theta_start)
            # beta_x 起始值：beta_x_hat 加上正态分布 N(0, sd^2) 的噪声
            noise <- matrix(rnorm(P * num_cols_beta_x, mean = 0, sd = sd_beta_x_start),
                nrow = P, ncol = num_cols_beta_x
            )
            current_rcpp_args$beta_x_init_nullable <- rcpp_args$beta_x_hat + noise
            run_desc <- "(随机起点)"
        }
        # print(current_rcpp_args)

        if (verbose) {
            cat(paste("  运行", i, "/", n_starts, run_desc, "... "))
        }

        # 使用 tryCatch 运行 Rcpp 优化，捕获可能的错误
        # print(current_rcpp_args)
        result_i <- tryCatch(
            {
                # 确保 C++ 函数已加载且可调用
                if (!exists("cMl_para_yuanshi_cpp", mode = "function")) {
                    stop("函数 'cMl_para_yuanshi_cpp' 未找到。您是否已编译/加载 C++ 代码？")
                }
                # 使用 do.call 将参数列表传递给 C++ 函数
                do.call(cMl_para_yuanshi_cpp, current_rcpp_args)
            },
            error = function(e) {
                # 捕获 C++ 调用期间或 C++ 代码内部（通过 Rcpp::stop）抛出的错误
                warning(paste("\n运行", i, "失败:", e$message), immediate. = TRUE, call. = FALSE)
                # 返回一个表示失败的结构体
                return(list(converged = FALSE, error = TRUE, message = e$message))
            }
        )

        # 存储本次运行的结果
        all_results[[i]] <- result_i

        # 如果需要，打印本次运行的状态
        if (verbose) {
            if (isTRUE(result_i$error)) {
                cat("失败。\n")
            } else if (!isTRUE(result_i$converged)) {
                cat("未收敛。\n")
            } else {
                # 如果成功收敛，打印 theta 值
                cat("收敛 (theta =", round(result_i$theta, digits_compare), ")。\n")
            }
        }
    }
    # 多起点循环结束

    if (verbose) {
        cat("多起点优化运行完成。\n--- 分析结果 ---\n")
    }

    # --- 3. 分析结果 ---
    # 筛选成功收敛的运行
    # 需要同时检查 error 标志和 converged 标志
    converged_results <- Filter(
        function(res) is.list(res) && !is.null(res$converged) && res$converged && !isTRUE(res$error),
        all_results
    )
    n_converged <- length(converged_results) # 成功收敛的次数
    # 计算失败的两种情况次数
    n_failed_error <- sum(sapply(all_results, function(res) isTRUE(res$error))) # 出错次数
    n_failed_nonconverge <- n_starts - n_converged - n_failed_error # 未收敛次数

    # 初始化摘要信息列表
    summary_info <- list(
        n_total_runs = n_starts, # 总运行次数
        n_converged = n_converged, # 成功收敛次数
        n_failed_error = n_failed_error, # 出错次数
        n_failed_nonconverge = n_failed_nonconverge, # 未收敛次数
        theta_freq_table = NULL, # theta 频率表
        selected_theta = NA, # 选定的 theta 值
        selected_theta_frequency = NA, # 选定 theta 的频率
        selected_theta_proportion = NA, # 选定 theta 的比例
        selected_run_index = NA, # 选定结果对应的运行索引
        message = "" # 摘要消息
    )

    # 处理所有运行都未成功收敛的情况
    if (n_converged == 0) {
        warning("所有优化运行均失败或未收敛。无法选择最优结果。", immediate. = TRUE, call. = FALSE)
        summary_info$message <- "所有运行均失败或未收敛。"
        # 返回一个与失败运行结构一致的占位列表
        final_result <- list(
            theta = NA, beta_x = NA, r = NA, theta_se = NA,
            beta_x_sparse = NA, iterations = NA, converged = FALSE,
            multi_start_summary = summary_info
            # , all_results = all_results # 可选：返回所有运行结果
        )
        return(final_result)
    }

    # 从成功收敛的运行中提取 theta 值，进行四舍五入和制表
    converged_thetas <- sapply(converged_results, function(res) res$theta)
    # 处理收敛运行中可能出现的非有限 theta 值 (理论上应该很少见)
    valid_thetas <- converged_thetas[is.finite(converged_thetas)]
    if (length(valid_thetas) < n_converged) {
        warning(
            paste(
                n_converged - length(valid_thetas),
                "次收敛运行的 theta 值为非有限值，已从频率分析中排除。"
            ),
            immediate. = TRUE, call. = FALSE
        )
        if (length(valid_thetas) == 0) {
            # 处理所有收敛的 theta 都是非有限值的情况
            warning("没有收敛运行产生有效的有限 theta 值。无法选择最优结果。", immediate. = TRUE, call. = FALSE)
            summary_info$message <- "没有收敛运行产生有效的有限 theta 值。"
            final_result <- list(
                theta = NA, beta_x = NA, r = NA, theta_se = NA,
                beta_x_sparse = NA, iterations = NA, converged = FALSE,
                multi_start_summary = summary_info
            )
            return(final_result)
        }
        # 更新用于频率计算的收敛次数和结果列表
        n_converged <- length(valid_thetas)
        converged_results <- converged_results[is.finite(converged_thetas)]
    }
    # 对有效的 theta 值进行四舍五入
    rounded_thetas <- round(valid_thetas, digits = digits_compare)


    # 计算频率表
    theta_freq <- table(rounded_thetas)
    theta_freq_sorted <- sort(theta_freq, decreasing = TRUE) # 按频率降序排序
    summary_info$theta_freq_table <- theta_freq_sorted

    # 找到众数 (最频繁出现的 theta 值)
    most_frequent_theta_rounded_str <- names(theta_freq_sorted)[1] # 频率最高的 theta (字符型)
    most_frequent_theta <- as.numeric(most_frequent_theta_rounded_str) # 转换回数值
    max_freq <- theta_freq_sorted[1] # 最高频率

    summary_info$selected_theta <- most_frequent_theta
    summary_info$selected_theta_frequency <- max_freq
    summary_info$selected_theta_proportion <- max_freq / n_converged # 计算比例

    # 在原始 all_results 列表中查找第一个收敛到众数的运行的索引
    selected_index <- NA
    for (k in 1:n_starts) {
        res_k <- all_results[[k]]
        # 检查是否是有效、收敛且 theta 有限的结果
        if (is.list(res_k) && !is.null(res_k$converged) && res_k$converged &&
            !isTRUE(res_k$error) && !is.null(res_k$theta) && is.finite(res_k$theta)) {
            # 比较四舍五入后的 theta 值
            if (round(res_k$theta, digits = digits_compare) == most_frequent_theta) {
                selected_index <- k # 记录索引
                break # 找到第一个就停止
            }
        }
    }

    # 处理罕见的找不到索引的情况 (理论上在 n_converged > 0 时不应发生)
    if (is.na(selected_index)) {
        warning("内部逻辑错误：未能找到选定的运行索引。将返回第一个有效的收敛结果。", immediate. = TRUE, call. = FALSE)
        # 查找第一个有效收敛结果的索引作为备选
        selected_index <- which(sapply(all_results, function(res) is.list(res) && !is.null(res$converged) && res$converged && !isTRUE(res$error) && !is.null(res$theta) && is.finite(res$theta)))[1]
        if (is.na(selected_index)) {
            # 如果仍然找不到，则为致命错误
            summary_info$message <- "内部错误：未能选择最终结果运行。"
            final_result <- list(
                theta = NA, beta_x = NA, r = NA, theta_se = NA,
                beta_x_sparse = NA, iterations = NA, converged = FALSE,
                multi_start_summary = summary_info
            )
            return(final_result)
        }
        summary_info$message <- "警告：由于索引查找问题，使用了备选选择。"
        # 基于备选运行更新选定的 theta
        summary_info$selected_theta <- all_results[[selected_index]]$theta
    }
    summary_info$selected_run_index <- selected_index

    # --- 4. 报告并返回结果 ---
    proportion_frequent <- summary_info$selected_theta_proportion
    # 检查结果稳定性
    if (proportion_frequent < min_freq_prop) {
        warning(sprintf(
            "结果稳定性低：最频繁的 theta (%.*f) 仅出现在 %.1f%% 的收敛运行中 (阈值 %.1f%%)。",
            digits_compare, most_frequent_theta, proportion_frequent * 100, min_freq_prop * 100
        ), immediate. = TRUE, call. = FALSE)
        summary_info$message <- "警告：稳定性低 - 最频繁 theta 低于阈值。"
    } else if (length(theta_freq_sorted) > 1) {
        # 如果频率达标，但存在其他收敛值
        if (verbose) {
            cat(sprintf(
                "最频繁的 theta (%.*f) 出现在 %.1f%% 的收敛运行中。已选择运行 %d 的结果。\n",
                digits_compare, most_frequent_theta, proportion_frequent * 100, selected_index
            ))
        }
        summary_info$message <- "多数运行收敛到相似结果。"
    } else {
        # 所有收敛运行结果一致
        if (verbose) {
            cat(sprintf(
                "所有 %d 次收敛运行均得到 theta ≈ %.*f。已选择运行 %d 的结果。\n",
                n_converged, digits_compare, most_frequent_theta, selected_index
            ))
        }
        summary_info$message <- "所有收敛运行结果一致。"
    }

    # 获取最终选定的结果列表
    final_result <- all_results[[selected_index]]
    # 将多起点运行的摘要信息添加到结果列表中
    final_result$multi_start_summary <- summary_info
    # 可选：添加所有运行结果以供详细检查
    # final_result$all_results <- all_results

    return(final_result) # 返回最终结果
}
# BIC 计算函数
BIC_function_optimized <- function(
    theta, k_x, k_par, beta_x, beta_x_hat, Sigma_inv_x,
    beta_y_hat, Sigma_inv_y, n_eff, r) {
    n_snps <- nrow(beta_x)
    if (n_snps == 0) {
        return(Inf)
    } # 处理空输入

    BIC_term1 <- 0
    BIC_term2 <- 0

    # 计算拟合优度项 (与 -2 * logLikelihood 相关)
    for (i in 1:n_snps) {
        idx_range <- (2 * i - 1):(2 * i) # 当前 SNP 的索引范围

        # X 部分的残差平方和 (加权)
        diff_x <- beta_x_hat[i, ] - beta_x[i, ]
        BIC_term1 <- BIC_term1 + (diff_x %*% Sigma_inv_x[idx_range, ] %*% diff_x)

        # Y 部分的残差平方和 (加权)
        predicted_y <- theta * beta_x[i, ] + r[i, ]
        diff_y <- beta_y_hat[i, ] - predicted_y
        BIC_term2 <- BIC_term2 + (diff_y %*% Sigma_inv_y[idx_range, ] %*% diff_y)
    }

    # 计算复杂度惩罚项
    # 惩罚与多效性相关的自由参数数量 (k_x 和 k_par)
    # 这个惩罚可以调整很多
    penalty <- log(n_eff) * (k_x + k_par) * 1 / 2

    # 计算总 BIC
    BIC <- BIC_term1 + BIC_term2 + penalty

    return(as.numeric(BIC)) # 确保返回标量数值
}


