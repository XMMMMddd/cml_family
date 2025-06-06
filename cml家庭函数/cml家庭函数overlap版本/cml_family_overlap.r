# %% 加载包
library(Matrix)
library(dplyr)
library(Rcpp) # 需要 Rcpp 包来编译和加载 C++ 代码

# %% 编译和加载 Rcpp 函数

# 确保 iterative_algorithm_rcpp.cpp 文件与此 R 脚本在正确的相对路径或提供绝对路径
cpp_file_path <- file.path("cml家庭函数", "cml家庭函数overlap版本", "iterative_algorithm_rcpp.cpp")

Rcpp::sourceCpp(cpp_file_path)
print(paste("成功加载 Rcpp 函数从:", cpp_file_path))



# %% 旧的R函数和其他辅助函数 (来自用户之前的文件内容)
source("FGWAS 优化版本/FGWAS 函数.R") # 假设此文件包含 perform_fgwas_analysis
source("样本生成函数/三联体家庭结构生成函数.R") # 假设此文件包含 generate_mr_trio_data_ultra_updata
# 对角阵处理算法
create_combined_diagonal_matrices <- function(
    beta_sigma_exp,
    beta_sigma_out,
    off_diagonal_elements = NULL) {
    if (nrow(beta_sigma_exp) != nrow(beta_sigma_out)) {
        stop("beta_sigma_exp 和 beta_sigma_out 的行数必须一致。")
    }
    n_snps <- nrow(beta_sigma_exp) / 2
    if (n_snps %% 1 != 0) {
        stop("beta_sigma_exp 的行数必须是 2 的倍数。")
    }
    if (!is.null(off_diagonal_elements) && length(off_diagonal_elements) != n_snps * 4) {
        stop("off_diagonal_elements 的长度必须是 n_snps * 4。")
    }
    all_diag_matrices <- list()
    if (n_snps > 0) {
        for (i in 1:n_snps) {
            exp_matrix <- matrix(beta_sigma_exp[(2 * i - 1):(2 * i), ], nrow = 2, ncol = 2)
            out_matrix <- matrix(beta_sigma_out[(2 * i - 1):(2 * i), ], nrow = 2, ncol = 2)
            current_block_matrix <- matrix(0, nrow = 4, ncol = 4)
            current_block_matrix[1:2, 1:2] <- exp_matrix
            current_block_matrix[3:4, 3:4] <- out_matrix
            if (!is.null(off_diagonal_elements)) {
                current_off_diag_block_elements <- off_diagonal_elements[((i - 1) * 4 + 1):((i - 1) * 4 + 4)]
                current_block_matrix[1, 3] <- current_off_diag_block_elements[1]
                current_block_matrix[1, 4] <- current_off_diag_block_elements[2]
                current_block_matrix[2, 3] <- current_off_diag_block_elements[3]
                current_block_matrix[2, 4] <- current_off_diag_block_elements[4]
                current_block_matrix[3, 1] <- current_block_matrix[1, 3]
                current_block_matrix[4, 1] <- current_block_matrix[1, 4]
                current_block_matrix[3, 2] <- current_block_matrix[2, 3]
                current_block_matrix[4, 2] <- current_block_matrix[2, 4]
            }
            all_diag_matrices[[i]] <- current_block_matrix
        }
    }
    if (length(all_diag_matrices) == 0) {
        return(matrix(0, 0, 4))
    }
    final_combined_matrix <- do.call(rbind, all_diag_matrices)
    return(final_combined_matrix)
}
# 处理算法
invert_beta_sigma_out_matrices <- function(beta_sigma_out) {
    n_snps <- nrow(beta_sigma_out) / 2
    if (n_snps %% 1 != 0) {
        stop("beta_sigma_out 的行数必须是 2 的倍数。")
    }
    inverted_beta_sigma_out <- matrix(NA, nrow = nrow(beta_sigma_out), ncol = ncol(beta_sigma_out))
    if (n_snps > 0) {
        for (i in 1:n_snps) {
            current_matrix <- matrix(beta_sigma_out[(2 * i - 1):(2 * i), ], nrow = 2, ncol = 2)
            if (det(current_matrix) == 0) {
                warning(paste("第", i, "个 2x2 矩阵是奇异的，无法求逆。将跳过此矩阵。"))
                inverted_beta_sigma_out[(2 * i - 1):(2 * i), ] <- NA
            } else {
                inverted_matrix <- solve(current_matrix)
                inverted_beta_sigma_out[(2 * i - 1):(2 * i), ] <- inverted_matrix
            }
        }
    }
    return(inverted_beta_sigma_out)
}
# 迭代算法
iterative_algorithm_r_version <- function(
    matrix_big, beta_sigma_exp, beta_sigma_out, beta_hat_exp, beta_hat_out, k,
    initial_alpha = 0,
    initial_beta_exp = NULL,
    initial_r = NULL,
    max_iter = 100,
    tolerance = 1e-6) {
    if (nrow(matrix_big) %% 4 != 0) stop("参数 'matrix_big' 的行数必须是 4 的倍数。")
    if (nrow(beta_hat_exp) != nrow(beta_hat_out)) stop("参数 'beta_hat_exp' 和 'beta_hat_out' 的行数必须一致。")
    if (k < 0 || k > nrow(beta_hat_exp)) stop("参数 'k' 必须在 0 和 SNP 总数之间 (包含 0 和 SNP 总数)。")
    if (max_iter <= 0) stop("参数 'max_iter' 必须是一个正整数。")
    if (tolerance <= 0) stop("参数 'tolerance' 必须是一个正数。")

    n_snps <- nrow(matrix_big) / 4
    alpha <- initial_alpha
    if (is.null(initial_beta_exp)) {
        beta_exp <- beta_hat_exp
    } else {
        beta_exp <- initial_beta_exp
    }
    if (is.null(initial_r)) {
        r <- matrix(0, nrow = n_snps, ncol = 2)
    } else {
        r <- initial_r
    }

    matrix_big_inverses <- list()
    if (n_snps > 0) {
        for (i in 1:n_snps) {
            matrix_index_all <- (4 * i - 3):(4 * i)
            matrix_big_i <- matrix_big[matrix_index_all, ]
            matrix_big_inverses[[i]] <- tryCatch(solve(matrix_big_i),
                error = function(e) {
                    warning(paste("R版本：位于 SNP", i, "的 4x4 矩阵是奇异的 - 将使用伪逆。"))
                    if (!requireNamespace("MASS", quietly = TRUE)) stop("需要 'MASS' 包。", call. = FALSE)
                    MASS::ginv(matrix_big_i)
                }
            )
        }
    }

    iterations_taken <- 0
    converged <- FALSE
    alpha_old <- alpha

    if (n_snps > 0) {
        for (iter_count in 1:max_iter) {
            iterations_taken <- iter_count
            alpha_old <- alpha
            d <- numeric(n_snps)
            r_new_iter <- matrix(0, nrow = n_snps, ncol = 2)
            for (i in 1:n_snps) {
                matrix_index_2x2 <- (2 * i - 1):(2 * i)
                beta_hat_exp_i <- beta_hat_exp[i, , drop = FALSE]
                beta_hat_out_i <- beta_hat_out[i, , drop = FALSE]
                beta_exp_i <- beta_exp[i, , drop = FALSE]
                beta_sigma_exp_i <- beta_sigma_exp[matrix_index_2x2, , drop = FALSE]
                beta_sigma_out_i <- beta_sigma_out[matrix_index_2x2, , drop = FALSE]
                beta_sigma_rho_i <- matrix_big[(4 * i - 3):(4 * i - 2), 3:4, drop = FALSE]
                residual_d <- beta_hat_out_i - alpha * beta_exp_i
                d[i] <- residual_d %*% beta_sigma_out_i %*% t(residual_d)
                exposure_residual_r <- beta_hat_exp_i - beta_exp_i
                r_i_calc <- residual_d - exposure_residual_r %*% t(beta_sigma_exp_i) %*% beta_sigma_rho_i
                r_new_iter[i, ] <- r_i_calc
            }
            d_order <- order(d, decreasing = TRUE)
            r_new_sorted <- r_new_iter[d_order, ]
            if (k < n_snps) {
                r_new_sorted[(k + 1):n_snps, ] <- 0
            } else if (k == n_snps) {
                NULL
            }
            r <- r_new_sorted[order(d_order), ]
            beta_exp_new_iter <- matrix(0, nrow = n_snps, ncol = 2)
            for (i in 1:n_snps) {
                omega_11 <- matrix_big_inverses[[i]][1:2, 1:2]
                omega_12 <- matrix_big_inverses[[i]][1:2, 3:4]
                omega_21 <- matrix_big_inverses[[i]][3:4, 1:2]
                omega_22 <- matrix_big_inverses[[i]][3:4, 3:4]
                r_i_current <- r[i, , drop = FALSE]
                fenmu_beta <- omega_11 + alpha * (omega_12 + omega_21) + alpha^2 * omega_22
                fenzi_beta <- (omega_11 + alpha * omega_21) %*% t(beta_hat_exp[i, , drop = FALSE]) +
                    (omega_12 + alpha * omega_22) %*% t(beta_hat_out[i, , drop = FALSE] - r_i_current)
                beta_exp_new_iter[i, ] <- t(solve(fenmu_beta, fenzi_beta))
            }
            beta_exp <- beta_exp_new_iter
            fenzi_alpha_total <- 0
            fenmu_alpha_total <- 0
            for (i in 1:n_snps) {
                omega_12 <- matrix_big_inverses[[i]][1:2, 3:4]
                omega_22 <- matrix_big_inverses[[i]][3:4, 3:4]
                beta_hat_exp_i <- beta_hat_exp[i, , drop = FALSE]
                beta_exp_i_current <- beta_exp[i, , drop = FALSE]
                beta_hat_out_i <- beta_hat_out[i, , drop = FALSE]
                r_i_current <- r[i, , drop = FALSE]
                fenzi_term_alpha <- (beta_hat_exp_i - beta_exp_i_current) %*% omega_12 %*% t(beta_exp_i_current) +
                    beta_exp_i_current %*% omega_22 %*% t(beta_hat_out_i - r_i_current)
                fenmu_term_alpha <- beta_exp_i_current %*% omega_22 %*% t(beta_exp_i_current)
                fenzi_alpha_total <- fenzi_alpha_total + fenzi_term_alpha
                fenmu_alpha_total <- fenmu_alpha_total + fenmu_term_alpha
            }
            if (abs(fenmu_alpha_total) < .Machine$double.eps) {
                warning(paste("R版本：在迭代", iter_count, "中，alpha分母接近零。"))
            } else {
                alpha <- as.numeric(fenzi_alpha_total / fenmu_alpha_total)
            }
            if (abs(alpha - alpha_old) < tolerance) {
                converged <- TRUE
                break
            }
        }
    } else {
        converged <- TRUE
        iterations_taken <- 0
    }
    if (iterations_taken == max_iter && !converged) {
        final_alpha_change <- abs(alpha - alpha_old)
        warning_message <- paste(
            "R版本：算法在", max_iter, "次迭代内未收敛。",
            "最后alpha变化:", format(final_alpha_change, digits = 3)
        )
        warning(warning_message)
    }
    return(list(alpha = alpha, beta_exp = beta_exp, r = r, iterations = iterations_taken, converged = converged))
}

# %% 多起点迭代 Rcpp 函数的 R 封装器
multi_start_iterative_algorithm_rcpp <- function(
    matrix_big, beta_sigma_exp, beta_sigma_out, beta_hat_exp, beta_hat_out, k,
    num_starts = 10, # 新增参数：指定随机起点的数量
    max_iter = 100,
    tolerance = 1e-6,
    alpha_group_tolerance = 1e-4,
    min_converged_fraction = 0.5) {
    if (!exists("iterative_algorithm_rcpp") || !is.function(iterative_algorithm_rcpp)) {
        stop("Rcpp 函数 'iterative_algorithm_rcpp' 未加载或不是一个函数。请检查 Rcpp::sourceCpp() 是否成功执行。")
    }
    if (num_starts <= 0) {
        stop("'num_starts' 必须是一个正整数。")
    }

    all_results <- list()
    n_snps <- nrow(beta_hat_exp) # 确保 n_snps 在循环外定义一次

    for (i in 1:num_starts) {
        # 自动生成起点参数
        initial_alpha_val <- runif(1, -0.5, 0.5) # 随机生成 alpha
        initial_beta_exp_val <- NULL # Rcpp 函数将处理默认初始化
        initial_r_val <- NULL # Rcpp 函数将处理默认初始化

        # 如果需要随机初始化 beta_exp 和 r，可以取消以下注释：
        # if (n_snps > 0) {
        #   initial_beta_exp_val <- matrix(rnorm(n_snps * 2, mean = 0, sd = 0.1), nrow = n_snps, ncol = 2)
        #   initial_r_val <- matrix(rnorm(n_snps * 2, mean = 0, sd = 0.01), nrow = n_snps, ncol = 2)
        # } else {
        #   initial_beta_exp_val <- matrix(0, nrow = 0, ncol = 2)
        #   initial_r_val <- matrix(0, nrow = 0, ncol = 2)
        # }

        cat(paste("运行自动生成的起点", i, "，初始 Alpha:", format(initial_alpha_val, digits = 4), "...\n"))
        res <- tryCatch(
            {
                iterative_algorithm_rcpp(
                    matrix_big = matrix_big,
                    beta_sigma_exp = beta_sigma_exp,
                    beta_sigma_out = beta_sigma_out,
                    beta_hat_exp = beta_hat_exp,
                    beta_hat_out = beta_hat_out,
                    k = as.integer(k),
                    initial_alpha = as.double(initial_alpha_val),
                    initial_beta_exp_r = initial_beta_exp_val,
                    initial_r_r = initial_r_val,
                    max_iter = as.integer(max_iter),
                    tolerance = as.double(tolerance)
                )
            },
            error = function(e) {
                warning(paste("起点", i, "运行出错:", e$message))
                list(alpha = NA, beta_exp = NA, r = NA, iterations = 0, converged = FALSE, error = e$message)
            }
        )
        all_results[[i]] <- res
    }

    converged_alphas <- sapply(all_results, function(res) if (res$converged) res$alpha else NA)
    valid_converged_alphas <- na.omit(converged_alphas)

    overall_success <- FALSE
    consensus_alpha <- NA
    consensus_beta_exp <- NA
    consensus_r <- NA
    num_successful_starts <- 0

    if (length(valid_converged_alphas) > 0) {
        sorted_alphas <- sort(valid_converged_alphas)
        alpha_groups <- list()
        if (length(sorted_alphas) > 0) {
            current_group <- c(sorted_alphas[1])
            for (j in seq_along(sorted_alphas)[-1]) {
                if (abs(sorted_alphas[j] - current_group[length(current_group)]) < alpha_group_tolerance) {
                    current_group <- c(current_group, sorted_alphas[j])
                } else {
                    alpha_groups[[length(alpha_groups) + 1]] <- current_group
                    current_group <- c(sorted_alphas[j])
                }
            }
            alpha_groups[[length(alpha_groups) + 1]] <- current_group
        }

        if (length(alpha_groups) > 0) {
            group_lengths <- sapply(alpha_groups, length)
            max_group_idx <- which.max(group_lengths)
            num_successful_starts <- group_lengths[max_group_idx]

            if (num_successful_starts / num_starts >= min_converged_fraction) {
                overall_success <- TRUE
                consensus_alpha_values_in_group <- alpha_groups[[max_group_idx]]
                consensus_alpha <- mean(consensus_alpha_values_in_group)

                for (res_idx in seq_along(all_results)) {
                    res_detail <- all_results[[res_idx]]
                    if (res_detail$converged && any(abs(res_detail$alpha - consensus_alpha_values_in_group) < alpha_group_tolerance)) {
                        consensus_beta_exp <- res_detail$beta_exp
                        consensus_r <- res_detail$r
                        break
                    }
                }
            }
        }
    }

    return(list(
        overall_success = overall_success,
        consensus_alpha = consensus_alpha,
        consensus_beta_exp = consensus_beta_exp,
        consensus_r = consensus_r,
        num_successful_starts = num_successful_starts,
        total_starts = num_starts,
        details = all_results
    ))
}
# %% 新增：计算 Alpha 方差的函数 (基于Rcpp多起点结果)
calculate_alpha_variance_from_rcpp_results <- function(
    consensus_alpha,
    consensus_beta_exp, # 这是在计算中使用的 beta_exp
    consensus_r,
    matrix_big_input,
    beta_sigma_exp_input,
    beta_sigma_out_input,
    beta_hat_exp_input,
    beta_hat_out_input) {
    # 从输入参数重命名/准备变量以匹配原始脚本逻辑
    alpha <- consensus_alpha
    beta_exp <- consensus_beta_exp
    r_input <- consensus_r # 避免与 base::r 冲突
    matrix_big <- matrix_big_input
    beta_sigma_exp <- beta_sigma_exp_input
    beta_sigma_out <- beta_sigma_out_input
    beta_hat_exp <- beta_hat_exp_input
    beta_hat_out <- beta_hat_out_input

    n_snps <- nrow(beta_hat_exp)

    snps_list_input <- which(r_input[, 1] == 0)
    n_snps_va <- length(snps_list_input)
    omega_alpha_alpha <- matrix(0, nrow = 1, ncol = 1) # 初始化为1x1零矩阵

    # 初始化 omega_alpha_gamma 和 omega_gamma_gamma，即使 n_snps_va 为 0
    omega_alpha_gamma <- matrix(0, ncol = 1, nrow = 2 * n_snps_va)
    omega_gamma_gamma <- matrix(0, ncol = 2 * n_snps_va, nrow = 2 * n_snps_va)

    matrix_big_inverses <- list()
    if (n_snps > 0) {
        for (i_loop_var in 1:n_snps) {
            matrix_index_all <- (4 * i_loop_var - 3):(4 * i_loop_var)
            matrix_big_i <- matrix_big[matrix_index_all, , drop = FALSE]
            matrix_big_inverses[[i_loop_var]] <- tryCatch(solve(matrix_big_i),
                error = function(e) {
                    warning(paste("R版本 (方差函数内)：位于 SNP", i_loop_var, "的 4x4 矩阵是奇异的 - 将使用伪逆。"))
                    if (!requireNamespace("MASS", quietly = TRUE)) stop("需要 'MASS' 包。", call. = FALSE)
                    MASS::ginv(matrix_big_i)
                }
            )
        }
    }

    if (n_snps_va > 0) {
        temp_omega_alpha_gamma_list <- vector("list", n_snps_va)
        temp_omega_gamma_gamma_list <- vector("list", n_snps_va)

        for (j_idx in seq_along(snps_list_input)) {
            i_snp_actual_idx <- snps_list_input[j_idx]

            beta_hat_exp_i_val <- t(beta_hat_exp[i_snp_actual_idx, , drop = FALSE])
            beta_hat_out_i_val <- t(beta_hat_out[i_snp_actual_idx, , drop = FALSE])
            beta_exp_i_val <- t(beta_exp[i_snp_actual_idx, , drop = FALSE])
            matrix_index_2x2_actual <- (2 * i_snp_actual_idx - 1):(2 * i_snp_actual_idx)
            beta_sigma_exp_i_val <- beta_sigma_exp[matrix_index_2x2_actual, , drop = FALSE]
            beta_sigma_out_i_val <- beta_sigma_out[matrix_index_2x2_actual, , drop = FALSE]
            beta_sigma_rho_i_val <- matrix_big[(4 * i_snp_actual_idx - 3):(4 * i_snp_actual_idx - 2), 3:4, drop = FALSE]

            current_matrix_big_inv <- matrix_big_inverses[[i_snp_actual_idx]]
            omega_11_val <- current_matrix_big_inv[1:2, 1:2, drop = FALSE]
            omega_12_val <- current_matrix_big_inv[1:2, 3:4, drop = FALSE]
            omega_21_val <- current_matrix_big_inv[3:4, 1:2, drop = FALSE]
            omega_22_val <- current_matrix_big_inv[3:4, 3:4, drop = FALSE]

            omega_alpha_alpha <- t(beta_exp_i_val) %*% omega_22_val %*% beta_exp_i_val + omega_alpha_alpha

            current_omega_alpha_gamma_i <- omega_12_val %*% beta_exp_i_val +
                omega_21_val %*% (beta_exp_i_val - beta_hat_exp_i_val) -
                omega_22_val %*% (beta_hat_out_i_val - 2 * alpha * beta_exp_i_val)
            temp_omega_alpha_gamma_list[[j_idx]] <- current_omega_alpha_gamma_i

            current_omega_gamma_gamma_i <- omega_11_val + alpha * (omega_12_val +
                t(omega_12_val)) + alpha^2 * omega_22_val
            temp_omega_gamma_gamma_list[[j_idx]] <- current_omega_gamma_gamma_i
        }

        if (length(temp_omega_alpha_gamma_list) > 0 && !all(sapply(temp_omega_alpha_gamma_list, is.null))) { # Check if list is not empty or full of NULLs
            omega_alpha_gamma <- do.call(rbind, temp_omega_alpha_gamma_list)
        }
        if (length(temp_omega_gamma_gamma_list) > 0 && !all(sapply(temp_omega_gamma_gamma_list, is.null))) { # Check if list is not empty or full of NULLs
            omega_gamma_gamma <- as.matrix(Matrix::bdiag(temp_omega_gamma_gamma_list))
        }
    }
    # 如果 n_snps_va == 0, omega_alpha_alpha 保持为 matrix(0,1,1),
    # omega_alpha_gamma 和 omega_gamma_gamma 保持为零维度矩阵

    H <- matrix(0, nrow = 2 * n_snps_va + 1, ncol = 2 * n_snps_va + 1)
    H[1, 1] <- omega_alpha_alpha
    if (n_snps_va > 0) {
        # Ensure dimensions match before assignment
        if (nrow(H[2:(2 * n_snps_va + 1), 1, drop = FALSE]) == nrow(omega_alpha_gamma) && ncol(H[2:(2 * n_snps_va + 1), 1, drop = FALSE]) == ncol(omega_alpha_gamma)) {
            H[2:(2 * n_snps_va + 1), 1] <- omega_alpha_gamma
        }
        if (nrow(H[1, 2:(2 * n_snps_va + 1), drop = FALSE]) == nrow(t(omega_alpha_gamma)) && ncol(H[1, 2:(2 * n_snps_va + 1), drop = FALSE]) == ncol(t(omega_alpha_gamma))) {
            H[1, 2:(2 * n_snps_va + 1)] <- t(omega_alpha_gamma)
        }
        if (all(dim(H[2:(2 * n_snps_va + 1), 2:(2 * n_snps_va + 1), drop = FALSE]) == dim(omega_gamma_gamma))) {
            H[2:(2 * n_snps_va + 1), 2:(2 * n_snps_va + 1)] <- omega_gamma_gamma
        }
    }

    alpha_var <- tryCatch(
        {
            solve(H)[1, 1]
        },
        error = function(e) {
            warning(paste("Error in solving H matrix for alpha_var:", e$message, "- H may be singular. Returning NA."))
            NA
        }
    )
    return(alpha_var)
}
# %% 新增：BIC 计算函数

# 函数功能：根据迭代算法的结果和其他输入计算BIC值。
# 参数：
#   test_results: multi_start_iterative_algorithm_rcpp 函数的输出列表，
#                 应包含 consensus_alpha, consensus_beta_exp, consensus_r。
#   n_snps_input: SNP 的数量 (数值型，正整数)。
#   matrix_big_input: 组合的 (n_snps_input * 4) x 4 矩阵，通常是协方差矩阵的逆的块组合。
#   beta_hat_exp_input: n_snps_input x 2 的暴露效应估计矩阵。
#   beta_hat_out_input: n_snps_input x 2 的结局效应估计矩阵。
#   k_input: 迭代算法中使用的 k 值 (数值型，非负整数，通常 <= n_snps_input)。
# 返回：
#   计算得到的 BIC 值 (单一数值型)。
calculate_bic_value <- function(test_results, n_snps_input, matrix_big_input, beta_hat_exp_input, beta_hat_out_input, k_input) {
    # 从 test_results 中提取所需的值
    alpha <- test_results$consensus_alpha
    beta_expose <- test_results$consensus_beta_exp
    r_values <- test_results$consensus_r # 使用 r_values 避免与 base::r 冲突

    # --- 输入验证 ---
    if (is.null(alpha) || !is.numeric(alpha) || length(alpha) != 1) {
        stop("test_results$consensus_alpha 必须是一个单一数值。")
    }
    if (is.null(beta_expose) || !is.matrix(beta_expose)) {
        stop("test_results$consensus_beta_exp 必须是一个矩阵。")
    }
    if (is.null(r_values) || !is.matrix(r_values)) {
        stop("test_results$consensus_r 必须是一个矩阵。")
    }
    if (!is.numeric(n_snps_input) || n_snps_input <= 0 || (n_snps_input %% 1 != 0)) {
        stop("n_snps_input 必须是一个正整数。")
    }
    if (nrow(beta_expose) != n_snps_input || ncol(beta_expose) != 2) {
        stop(paste("beta_expose 的维度应为", n_snps_input, "x 2. 当前维度:", nrow(beta_expose), "x", ncol(beta_expose)))
    }
    if (nrow(r_values) != n_snps_input || ncol(r_values) != 2) {
        stop(paste("r_values 的维度应为", n_snps_input, "x 2. 当前维度:", nrow(r_values), "x", ncol(r_values)))
    }
    if (!is.matrix(beta_hat_exp_input) || nrow(beta_hat_exp_input) != n_snps_input || ncol(beta_hat_exp_input) != 2) {
        stop(paste("beta_hat_exp_input 必须是一个", n_snps_input, "x 2 的矩阵. 当前维度:", nrow(beta_hat_exp_input), "x", ncol(beta_hat_exp_input)))
    }
    if (!is.matrix(beta_hat_out_input) || nrow(beta_hat_out_input) != n_snps_input || ncol(beta_hat_out_input) != 2) {
        stop(paste("beta_hat_out_input 必须是一个", n_snps_input, "x 2 的矩阵. 当前维度:", nrow(beta_hat_out_input), "x", ncol(beta_hat_out_input)))
    }
    if (!is.matrix(matrix_big_input) || nrow(matrix_big_input) != n_snps_input * 4 || ncol(matrix_big_input) != 4) {
        stop(paste("matrix_big_input 必须是一个 (", n_snps_input, "* 4) x 4 的矩阵. 当前维度:", nrow(matrix_big_input), "x", ncol(matrix_big_input)))
    }
    if (!is.numeric(k_input) || k_input < 0 || (k_input %% 1 != 0) || k_input > n_snps_input) {
        stop(paste("k_input (", k_input, ") 必须是一个介于 0 和 n_snps_input (", n_snps_input, ") 之间的整数。"))
    }
    # --- 输入验证结束 ---

    l_k_val <- 0 # 初始化 l_k (似然函数的组成部分)
    for (i in 1:n_snps_input) {
        matrix_index_all <- (4 * i - 3):(4 * i)
        matrix_big_i <- matrix_big_input[matrix_index_all, , drop = FALSE] # 应为 4x4 矩阵

        # 尝试求逆，如果奇异则使用广义逆
        matrix_big_i_inv <- tryCatch(
            {
                solve(matrix_big_i)
            },
            error = function(e) {
                warning(paste("matrix_big_i 在 SNP", i, "处是奇异的 (错误:", e$message, "). 尝试使用 MASS::ginv()."))
                if (!requireNamespace("MASS", quietly = TRUE)) {
                    stop("需要 'MASS' 包来处理奇异矩阵 (ginv)。请先执行 install.packages('MASS')")
                }
                MASS::ginv(matrix_big_i)
            }
        )

        beta_hat_exp_i <- beta_hat_exp_input[i, ]
        beta_hat_out_i <- beta_hat_out_input[i, ]
        beta_exp_i <- beta_expose[i, ]
        r_i <- r_values[i, ]

        d_exp_i <- beta_hat_exp_i - beta_exp_i
        d_out_i <- beta_hat_out_i - alpha * beta_exp_i - r_i
        d_i <- c(d_exp_i, d_out_i) # 这是一个长度为4的向量

        if (!is.numeric(d_i) || length(d_i) != 4) {
            stop(paste("d_i 在 SNP", i, "处计算错误，不是长度为4的数值向量。当前值:", paste(d_i, collapse = ", ")))
        }

        term_value <- tryCatch(
            {
                as.numeric(t(d_i) %*% matrix_big_i_inv %*% d_i) # 确保结果是单个数值
            },
            error = function(e_mult) {
                stop(paste("在 SNP", i, "处计算 t(d_i) %*% matrix_big_i_inv %*% d_i 时出错:", e_mult$message))
            }
        )

        if (!is.numeric(term_value) || length(term_value) != 1) {
            stop(paste("在 SNP", i, "处，项 t(d_i) %*% matrix_big_i_inv %*% d_i 未返回单个数值。当前值:", term_value))
        }
        l_k_val <- l_k_val + term_value
    }

    # BIC 计算遵循用户提供的原始公式: (1 + 2*n_snps + 2*k) + l_k
    # 参数数量的解释:
    # 1: for alpha
    # 2 * n_snps_input: for beta_expose (n_snps_input 个 SNP, 每个 SNP 有2个 beta_exp 参数)
    # 2 * k_input: for non-zero r_values (k_input 个 SNP 的 r_values 是非零的, 每个 r_value 有2个参数)
    # 总参数数量 = 1 + 2 * n_snps_input + 2 * k_input

    bic_value <- (1 + 2 * n_snps_input + 2 * k_input) + l_k_val

    return(as.numeric(bic_value)) # 确保最终返回的是单个数值
}

# %% BIC计算 代码版
test_data <- generate_mr_trio_data_ultra_updata(n_snps = 3)
test_fgwas <- perform_fgwas_analysis(test_data, 3)
beta_sigma_exp <- test_fgwas$beta_sigma_exp
beta_sigma_out <- test_fgwas$beta_sigma_out
beta_hat_exp <- test_fgwas$beta_hat_exp
beta_hat_out <- test_fgwas$beta_hat_out
beta_sigma_exp_inv <- invert_beta_sigma_out_matrices(beta_sigma_exp)
beta_sigma_out_inv <- invert_beta_sigma_out_matrices(beta_sigma_out)
matrix_big <- create_combined_diagonal_matrices(beta_sigma_exp_inv, beta_sigma_out_inv)

test <- multi_start_iterative_algorithm_rcpp(
    matrix_big, beta_sigma_exp, beta_sigma_out, beta_hat_exp, beta_hat_out,
    k = 0,
    num_starts = 3
)
k <- 0
alpha <- test$consensus_alpha
beta_expose <- test$consensus_beta_exp
r <- test$consensus_r
matrix_big
beta_sigma_exp
beta_sigma_out
beta_hat_exp
beta_hat_out
l_k <- 0 # 似然函数
for (i in 1:n_snps) {
    matrix_index_all <- (4 * i - 3):(4 * i)
    matrix_big_i <- matrix_big[matrix_index_all, ]
    matrix_big_i_inv <- solve(matrix_big_i)
    beta_hat_exp_i <- beta_hat_exp[i, ]
    beta_hat_out_i <- beta_hat_out[i, ]
    beta_exp_i <- beta_expose[i, ]
    r_i <- r[i, ]
    d_exp_i <- beta_hat_exp_i - beta_exp_i
    d_out_i <- beta_hat_out_i - alpha * beta_exp_i - r_i
    d_i <- c(d_exp_i, d_out_i)
    l_k <- t(d_i) %*% matrix_big_i_inv %*% d_i + l_k
}
bic <- (1 + 2 * n_snps + 2 * k) + l_k

# %% 方差计算函数

# test <- multi_start_iterative_algorithm_rcpp(
#     matrix_big, beta_sigma_exp, beta_sigma_out, beta_hat_exp, beta_hat_out, k,
#     num_starts = 3
# )
# alpha <- test$consensus_alpha
# beta_expose <- test$consensus_beta_exp
# r <- test$consensus_r
# matrix_big
# beta_sigma_exp
# beta_sigma_out
# beta_hat_exp
# beta_hat_out
# # 声明协方差矩阵
# snps_list_input <- which(r[, 1] == 0) # 拿到计算方差的snps序列
# n_snps_va <- length(snps_list_input)
# omega_alpha_alpha <- matrix(0)
# omega_alpha_gamma <- matrix(0, ncol = 1, nrow = 2 * n_snps_va)
# omega_gamma_gamma <- matrix(0, ncol = 2 * n_snps_va, nrow = 2 * n_snps_va)
# matrix_big_inverses <- list() # 拿到大矩阵的逆
# if (n_snps > 0) {
#     for (i in 1:n_snps) {
#         matrix_index_all <- (4 * i - 3):(4 * i)
#         matrix_big_i <- matrix_big[matrix_index_all, ]
#         matrix_big_inverses[[i]] <- tryCatch(solve(matrix_big_i),
#             error = function(e) {
#                 # %% 调用新封装的方差计算函数示例
#                 # 注意：此示例依赖于在前面的 '# %% 方差计算函数' 块中定义的 'test' 对象
#                 # 以及在全局环境中定义的 matrix_big, beta_sigma_exp 等变量。
#                 # 这些全局变量通常在脚本末尾的 "脚本部分：数据准备和函数调用示例" 中被定义为 *_ex 版本。
#                 # 为确保此示例可独立运行或与该部分协调，请确保变量名和来源正确。

#                 if (exists("test") && is.list(test) && !is.null(test$consensus_alpha) &&
#                     exists("matrix_big") && exists("beta_sigma_exp") && exists("beta_sigma_out") &&
#                     exists("beta_hat_exp") && exists("beta_hat_out")) {
#                     print("--- 调用新封装的 calculate_alpha_variance_from_rcpp_results 函数 ---")
#                     # 确保从 'test' 对象和全局变量中获取正确的参数
#                     alpha_var_from_function <- calculate_alpha_variance_from_rcpp_results(
#                         consensus_alpha = test$consensus_alpha,
#                         consensus_beta_exp = test$consensus_beta_exp, # 在原始脚本中，beta_expose 被赋值为此值
#                         consensus_r = test$consensus_r,
#                         matrix_big_input = matrix_big,
#                         beta_sigma_exp_input = beta_sigma_exp,
#                         beta_sigma_out_input = beta_sigma_out,
#                         beta_hat_exp_input = beta_hat_exp,
#                         beta_hat_out_input = beta_hat_out
#                     )
#                     print(paste("使用新函数计算得到的 Alpha 方差:", alpha_var_from_function))

#                     # 与原始脚本块计算得到的 alpha_var (如果存在) 进行比较
#                     if (exists("alpha_var")) {
#                         print(paste("原始脚本块直接计算得到的 Alpha 方差:", alpha_var))
#                     }
#                 } else {
#                     print("必要的变量 (如 'test' 对象或相关的矩阵数据) 未在当前作用域定义，跳过新方差函数的示例调用。")
#                     print("请确保 '# %% 方差计算函数' 块已成功执行，并且所需的数据矩阵已加载。")
#                 }
#                 warning(paste("R版本：位于 SNP", i, "的 4x4 矩阵是奇异的 - 将使用伪逆。"))
#                 if (!requireNamespace("MASS", quietly = TRUE)) stop("需要 'MASS' 包。", call. = FALSE)
#                 MASS::ginv(matrix_big_i)
#             }
#         )
#     }
# }

# for (i in snps_list_input) {
#     matrix_index_2x2 <- (2 * i - 1):(2 * i)
#     beta_hat_exp_i <- t(beta_hat_exp[i, , drop = FALSE])
#     beta_hat_out_i <- t(beta_hat_out[i, , drop = FALSE])
#     beta_exp_i <- t(beta_exp[i, , drop = FALSE])
#     beta_sigma_exp_i <- beta_sigma_exp[matrix_index_2x2, , drop = FALSE]
#     beta_sigma_out_i <- beta_sigma_out[matrix_index_2x2, , drop = FALSE]
#     beta_sigma_rho_i <- matrix_big[(4 * i - 3):(4 * i - 2), 3:4, drop = FALSE]
#     omega_11 <- matrix_big_inverses[[i]][1:2, 1:2]
#     omega_12 <- matrix_big_inverses[[i]][1:2, 3:4]
#     omega_21 <- matrix_big_inverses[[i]][3:4, 1:2]
#     omega_22 <- matrix_big_inverses[[i]][3:4, 3:4]
#     omega_alpha_alpha <- t(beta_exp_i) %*% omega_22 %*%
#         beta_exp_i + omega_alpha_alpha
#     omega_alpha_gamma_i <- omega_12 %*% beta_exp_i +
#         omega_21 %*% (beta_exp_i - beta_hat_exp_i) -
#         omega_22 %*% (beta_hat_out_i - 2 * alpha * beta_exp_i)
#     omega_alpha_gamma[matrix_index_2x2, ] <- omega_alpha_gamma_i
#     omega_gamma_gamma_i <- omega_11 + alpha * (omega_12 +
#         t(omega_12)) + alpha^2 * omega_22
#     omega_gamma_gamma[matrix_index_2x2, matrix_index_2x2] <- omega_gamma_gamma_i
# }
# H <- matrix(0, nrow = 2 * n_snps_va + 1, ncol = 2 * n_snps_va + 1)
# H[1, 1] <- omega_alpha_alpha
# H[2:(2 * n_snps_va + 1), 1] <- omega_alpha_gamma
# H[1, 2:(2 * n_snps_va + 1)] <- t(omega_alpha_gamma)
# H[2:(2 * n_snps_va + 1), 2:(2 * n_snps_va + 1)] <- omega_gamma_gamma
# alpha_var <- solve(H)[1, 1]

# test_var <- calculate_alpha_variance_from_rcpp_results(
#     alpha,
#     beta_exp, # 这是在计算中使用的 beta_exp
#     r,
#     matrix_big,
#     beta_sigma_exp,
#     beta_sigma_out,
#     beta_hat_exp,
#     beta_hat_out
# )

# %% 测试方差函数print("--- 开始脚本执行 ---")

n_snps_example <- 5
print(paste("示例 SNP 数量 (n_snps_example):", n_snps_example))

if (exists("generate_mr_trio_data_ultra_updata") && is.function(generate_mr_trio_data_ultra_updata)) {
    data_set_example <- generate_mr_trio_data_ultra_updata(beta_os_oe_exp = 0.5, n_snps = n_snps_example)
} else {
    stop("函数 'generate_mr_trio_data_ultra_updata' 未定义。")
}

if (exists("perform_fgwas_analysis") && is.function(perform_fgwas_analysis)) {
    data_set_fgwas_example <- perform_fgwas_analysis(data_set_example, n_snps = n_snps_example)
} else {
    stop("函数 'perform_fgwas_analysis' 未定义。")
}

beta_sigma_exp_ex <- as.matrix(data_set_fgwas_example$beta_sigma_exp)
beta_sigma_out_ex <- as.matrix(data_set_fgwas_example$beta_sigma_out)
beta_hat_exp_ex <- as.matrix(data_set_fgwas_example$beta_hat_exp)
beta_hat_out_ex <- as.matrix(data_set_fgwas_example$beta_hat_out)
k_ex <- 0

beta_sigma_t_exp_ex <- invert_beta_sigma_out_matrices(beta_sigma_exp_ex)
beta_sigma_t_out_ex <- invert_beta_sigma_out_matrices(beta_sigma_out_ex)
matrix_big_ex <- create_combined_diagonal_matrices(
    beta_sigma_t_exp_ex,
    beta_sigma_t_out_ex,
    off_diagonal_elements = NULL
)
print("示例数据准备完毕。")

default_r_ex <- matrix(0, nrow = n_snps_example, ncol = 2)

# # 构造起点列表，分解复杂行 (现在由函数内部自动生成)
# start_point_1 <- list(alpha = 0.0, beta_exp = NULL, r = NULL)
# start_point_2 <- list(alpha = 0.2, beta_exp = NULL, r = NULL)
# start_point_3 <- list(alpha = -0.2, beta_exp = NULL, r = NULL)
#
# perturbed_beta_exp_val <- beta_hat_exp_ex + matrix(rnorm(n_snps_example * 2, 0, 0.01), nrow = n_snps_example, ncol = 2)
# perturbed_r_val <- default_r_ex + matrix(rnorm(n_snps_example * 2, 0, 0.01), nrow = n_snps_example, ncol = 2)
# start_point_4 <- list(alpha = 0.5, beta_exp = perturbed_beta_exp_val, r = perturbed_r_val)
#
# start_point_5 <- list(alpha = -0.5, beta_exp = NULL, r = NULL)
#
# start_points_example <- list(
# start_point_1,
# start_point_2,
# start_point_3,
# start_point_4,
# start_point_5
# )
# print(paste("定义了", length(start_points_example), "个多起点。")) # 不再需要手动定义

print("--- 调用多起点 Rcpp 迭代函数 (自动生成起点) ---")
if (exists("iterative_algorithm_rcpp") && is.function(iterative_algorithm_rcpp) && !identical(body(iterative_algorithm_rcpp), body(function(...) {
    stop("Rcpp 函数 iterative_algorithm_rcpp未能成功加载。请检查路径和编译。")
}))) {
    multi_start_results <- multi_start_iterative_algorithm_rcpp(
        matrix_big = matrix_big_ex,
        beta_sigma_exp = beta_sigma_exp_ex,
        beta_sigma_out = beta_sigma_out_ex,
        beta_hat_exp = beta_hat_exp_ex,
        beta_hat_out = beta_hat_out_ex,
        k = k_ex,
        num_starts = 5, # 指定自动生成5个起点
        max_iter = 200,
        tolerance = 1e-7,
        alpha_group_tolerance = 1e-5,
        min_converged_fraction = 0.6
    )

    print("--- 多起点迭代结果 (自动生成起点) ---")
    print(paste("整体是否成功:", multi_start_results$overall_success))
    if (multi_start_results$overall_success) {
        print(paste("共识 Alpha:", multi_start_results$consensus_alpha))
    }
    print(paste("收敛到共识点的起点数量:", multi_start_results$num_successful_starts, "总起点数:", multi_start_results$total_starts))
} else {
    print("Rcpp 函数 iterative_algorithm_rcpp 未能加载，跳过多起点 Rcpp 测试。")
}

print("--- 脚本执行完毕 ---")
