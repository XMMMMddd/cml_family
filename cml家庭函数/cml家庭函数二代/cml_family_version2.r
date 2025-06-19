# %% 加载必要的包
library(data.table)
# %% 函数结果

run_iterative_mle_v2 <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, a, b,
    max_iter = 100, tol = 1e-6, alpha_init = 0) {
    num_snps <- nrow(beta_hat_exp)

    # --- 1. 数据预处理 ---
    gamma_hat <- t(beta_hat_exp)
    beta_hat <- t(beta_hat_out)

    Sigma_gamma_inv_list <- vector("list", num_snps)
    Sigma_beta_inv_list <- vector("list", num_snps)
    Sigma_gamma_list <- vector("list", num_snps)
    Sigma_beta_list <- vector("list", num_snps)

    for (m in 1:num_snps) {
        Sigma_gamma_inv_list[[m]] <- beta_sigma_exp[(2 * m - 1):(2 * m), ]
        Sigma_beta_inv_list[[m]] <- beta_sigma_out[(2 * m - 1):(2 * m), ]
        Sigma_gamma_list[[m]] <- solve(Sigma_gamma_inv_list[[m]])
        Sigma_beta_list[[m]] <- solve(Sigma_beta_inv_list[[m]])
    }

    # --- 2. 初始化参数 ---
    alpha_current <- alpha_init
    gamma_current <- gamma_hat

    converged <- FALSE
    final_sets <- list()

    # ✨ 新增：初始化一个向量来追踪 alpha 的历史值
    alpha_history <- rep(NA, max_iter)


    # --- 3. 迭代主循环 ---
    for (iter in 1:max_iter) {
        alpha_prev <- alpha_current

        # ==========================================================
        # 步骤 1: 选择 SNPs 集合 A 与 B
        # ==========================================================
        var_beta_o <- sapply(Sigma_beta_list, function(mat) mat[1, 1])
        var_beta_p <- sapply(Sigma_beta_list, function(mat) mat[2, 2])

        t_m <- (beta_hat[1, ] - alpha_current * gamma_current[1, ])^2 / var_beta_o
        f_m <- (beta_hat[2, ] - alpha_current * gamma_current[2, ])^2 / var_beta_p

        a <- min(a, num_snps)
        b <- min(b, num_snps)

        set_A_indices <- order(t_m)[1:a]
        set_B_indices <- order(f_m)[1:b]

        A_intersect_B <- intersect(set_A_indices, set_B_indices)
        A_diff_B <- setdiff(set_A_indices, set_B_indices)
        B_diff_A <- setdiff(set_B_indices, set_A_indices)
        others <- setdiff(1:num_snps, union(set_A_indices, set_B_indices))

        # ==========================================================
        # 步骤 2: 更新 gamma
        # ==========================================================
        gamma_new <- matrix(0, nrow = 2, ncol = num_snps)

        if (length(A_intersect_B) > 0) {
            for (m in A_intersect_B) {
                term1_inv <- Sigma_gamma_inv_list[[m]]
                term2_inv <- alpha_current^2 * Sigma_beta_inv_list[[m]]
                combined_inv <- term1_inv + term2_inv
                term1_val <- term1_inv %*% gamma_hat[, m]
                term2_val <- alpha_current * Sigma_beta_inv_list[[m]] %*% beta_hat[, m]
                gamma_new[, m] <- solve(combined_inv) %*% (term1_val + term2_val)
            }
        }

        if (length(A_diff_B) > 0) {
            for (m in A_diff_B) {
                num_go <- gamma_hat[1, m] * var_beta_o[m] - alpha_current * Sigma_gamma_list[[m]][1, 1] * beta_hat[1, m]
                den_go <- var_beta_o[m] - alpha_current^2 * Sigma_gamma_list[[m]][1, 1]
                gamma_o_new <- num_go / den_go
                term_gp <- Sigma_gamma_list[[m]][1, 2] * (-alpha_current / var_beta_o[m]) * (beta_hat[1, m] - alpha_current * gamma_o_new)
                gamma_p_new <- gamma_hat[2, m] + term_gp
                gamma_new[, m] <- c(gamma_o_new, gamma_p_new)
            }
        }

        if (length(B_diff_A) > 0) {
            for (m in B_diff_A) {
                num_gp <- gamma_hat[2, m] * var_beta_p[m] - alpha_current * Sigma_gamma_list[[m]][2, 2] * beta_hat[2, m]
                den_gp <- var_beta_p[m] - alpha_current^2 * Sigma_gamma_list[[m]][2, 2]
                gamma_p_new <- num_gp / den_gp
                term_go <- Sigma_gamma_list[[m]][1, 2] * (-alpha_current / var_beta_p[m]) * (beta_hat[2, m] - alpha_current * gamma_p_new)
                gamma_o_new <- gamma_hat[1, m] + term_go
                gamma_new[, m] <- c(gamma_o_new, gamma_p_new)
            }
        }

        if (length(others) > 0) {
            gamma_new[, others] <- gamma_hat[, others]
        }

        gamma_current <- gamma_new

        # ==========================================================
        # 步骤 3: 更新 alpha
        # ==========================================================
        num_alpha <- 0
        den_alpha <- 0

        if (length(A_intersect_B) > 0) {
            for (m in A_intersect_B) {
                num_alpha <- num_alpha + t(gamma_current[, m]) %*% Sigma_beta_inv_list[[m]] %*% beta_hat[, m]
                den_alpha <- den_alpha + t(gamma_current[, m]) %*% Sigma_beta_inv_list[[m]] %*% gamma_current[, m]
            }
        }
        if (length(A_diff_B) > 0) {
            for (m in A_diff_B) {
                num_alpha <- num_alpha + (gamma_current[1, m] * beta_hat[1, m] / var_beta_o[m])
                den_alpha <- den_alpha + (gamma_current[1, m]^2 / var_beta_o[m])
            }
        }
        if (length(B_diff_A) > 0) {
            for (m in B_diff_A) {
                num_alpha <- num_alpha + (gamma_current[2, m] * beta_hat[2, m] / var_beta_p[m])
                den_alpha <- den_alpha + (gamma_current[2, m]^2 / var_beta_p[m])
            }
        }

        if (abs(den_alpha) > 1e-9) {
            alpha_current <- as.numeric(num_alpha / den_alpha)
        } else {
            warning(paste("Iteration", iter, ": Denominator for alpha is near zero. Halting."))
            break
        }

        # ✨ 新增：在每次迭代后记录 alpha 的值
        alpha_history[iter] <- alpha_current

        # --- 4. 检查收敛 ---
        if (abs(alpha_current - alpha_prev) < tol) {
            converged <- TRUE
            final_sets <- list(A = set_A_indices, B = set_B_indices)
            message(paste("✅ Algorithm converged after", iter, "iterations."))
            break
        }
    }

    if (!converged) {
        warning(paste("⚠️ Algorithm did not converge after", max_iter, "iterations."))
        final_sets <- list(A = set_A_indices, B = set_B_indices)
    }

    # ✨ 新增：清理历史记录，移除未使用的 NA 值
    alpha_history <- alpha_history[!is.na(alpha_history)]

    return(list(
        alpha_final = as.numeric(alpha_current),
        gamma_final = t(gamma_current),
        iterations = iter,
        converged = converged,
        sets = final_sets,
        alpha_history = alpha_history[!is.na(alpha_history)],
        # ✨ 将所有预处理的列表一并返回
        data = list(
            beta_hat_exp = beta_hat_exp,
            beta_hat_out = beta_hat_out,
            Sigma_gamma_list = Sigma_gamma_list,
            Sigma_beta_list = Sigma_beta_list,
            Sigma_gamma_inv_list = Sigma_gamma_inv_list,
            Sigma_beta_inv_list = Sigma_beta_inv_list
        )
    ))
}

if (!require(Rcpp)) install.packages("Rcpp")
if (!require(RcppArmadillo)) install.packages("RcppArmadillo")

# 编译函数
Rcpp::sourceCpp("cml家庭函数/cml家庭函数二代/cml_family_rcpp.cpp")

Rcpp::sourceCpp("cml家庭函数/cml家庭函数二代/cml_family_rcpp_gamma.cpp")

run_iterative_mle_v2_rcpp <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, a, b,
    max_iter = 100, tol = 1e-6, alpha_init = 0) {
    # 输入验证
    if (!is.matrix(beta_hat_exp) || !is.matrix(beta_hat_out)) {
        stop("beta_hat_exp and beta_hat_out must be matrices")
    }

    if (nrow(beta_hat_exp) != nrow(beta_hat_out)) {
        stop("beta_hat_exp and beta_hat_out must have the same number of rows (SNPs)")
    }

    if (ncol(beta_hat_exp) != 2 || ncol(beta_hat_out) != 2) {
        stop("beta_hat_exp and beta_hat_out must have exactly 2 columns")
    }

    num_snps <- nrow(beta_hat_exp)

    if (a <= 0 || b <= 0) {
        stop("Parameters a and b must be positive integers")
    }

    if (max_iter <= 0) {
        stop("max_iter must be a positive integer")
    }

    if (tol <= 0) {
        stop("tol must be a positive number")
    }



    # --- 1. 数据预处理 (R语言部分) ---
    start_time <- Sys.time()

    gamma_hat <- t(beta_hat_exp)
    beta_hat <- t(beta_hat_out)

    Sigma_gamma_inv_list <- vector("list", num_snps)
    Sigma_beta_inv_list <- vector("list", num_snps)
    Sigma_gamma_list <- vector("list", num_snps)
    Sigma_beta_list <- vector("list", num_snps)



    for (m in 1:num_snps) {
        # 提取对应的2x2协方差矩阵
        Sigma_gamma_inv_list[[m]] <- beta_sigma_exp[(2 * m - 1):(2 * m), ]
        Sigma_beta_inv_list[[m]] <- beta_sigma_out[(2 * m - 1):(2 * m), ]

        # 计算逆矩阵，添加数值稳定性检查
        tryCatch(
            {
                Sigma_gamma_list[[m]] <- solve(Sigma_gamma_inv_list[[m]])
                Sigma_beta_list[[m]] <- solve(Sigma_beta_inv_list[[m]])
            },
            error = function(e) {
                # 如果矩阵奇异，尝试添加小的正则化项
                warning(paste("Singular matrix encountered at SNP", m, ". Adding regularization."))
                reg_term <- 1e-8 * diag(2)
                Sigma_gamma_list[[m]] <- solve(Sigma_gamma_inv_list[[m]] + reg_term)
                Sigma_beta_list[[m]] <- solve(Sigma_beta_inv_list[[m]] + reg_term)
            }
        )
    }

    preprocessing_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))


    # --- 2. 调用C++核心迭代函数 ---

    core_start_time <- Sys.time()

    tryCatch(
        {
            result <- iterative_mle_core(
                gamma_hat = gamma_hat,
                beta_hat = beta_hat,
                Sigma_gamma_list = Sigma_gamma_list,
                Sigma_beta_list = Sigma_beta_list,
                Sigma_gamma_inv_list = Sigma_gamma_inv_list,
                Sigma_beta_inv_list = Sigma_beta_inv_list,
                a = a, b = b,
                max_iter = max_iter,
                tol = tol,
                alpha_init = alpha_init
            )
        },
        error = function(e) {
            stop(paste("Error in C++ core function:", e$message))
        }
    )

    core_time <- as.numeric(difftime(Sys.time(), core_start_time, units = "secs"))



    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))


    final_sets <- list(
        A = as.integer(result$set_A_indices),
        B = as.integer(result$set_B_indices)
    )

    # 输出关键结果


    return(list(
        alpha_final = as.numeric(result$alpha_final),
        gamma_final = result$gamma_final,
        iterations = result$iterations,
        converged = result$converged,
        sets = final_sets,
        alpha_history = as.numeric(result$alpha_history),
        # 返回预处理的数据，便于进一步分析
        data = list(
            beta_hat_exp = beta_hat_exp,
            beta_hat_out = beta_hat_out,
            Sigma_gamma_list = Sigma_gamma_list,
            Sigma_beta_list = Sigma_beta_list,
            Sigma_gamma_inv_list = Sigma_gamma_inv_list,
            Sigma_beta_inv_list = Sigma_beta_inv_list
        ),
        # 添加计算时间信息
        timing = list(
            preprocessing_time = preprocessing_time,
            core_time = core_time,
            total_time = total_time
        )
    ))
}
run_iterative_mle_v2_rcpp_gamma <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, a, b,
    max_iter = 100, tol = 1e-6, alpha_init = 0) {
    # 输入验证
    if (!is.matrix(beta_hat_exp) || !is.matrix(beta_hat_out)) {
        stop("beta_hat_exp and beta_hat_out must be matrices")
    }

    if (nrow(beta_hat_exp) != nrow(beta_hat_out)) {
        stop("beta_hat_exp and beta_hat_out must have the same number of rows (SNPs)")
    }

    if (ncol(beta_hat_exp) != 2 || ncol(beta_hat_out) != 2) {
        stop("beta_hat_exp and beta_hat_out must have exactly 2 columns")
    }

    num_snps <- nrow(beta_hat_exp)

    if (a <= 0 || b <= 0) {
        stop("Parameters a and b must be positive integers")
    }

    if (max_iter <= 0) {
        stop("max_iter must be a positive integer")
    }

    if (tol <= 0) {
        stop("tol must be a positive number")
    }



    # --- 1. 数据预处理 (R语言部分) ---
    start_time <- Sys.time()

    gamma_hat <- t(beta_hat_exp)
    beta_hat <- t(beta_hat_out)

    Sigma_gamma_inv_list <- vector("list", num_snps)
    Sigma_beta_inv_list <- vector("list", num_snps)
    Sigma_gamma_list <- vector("list", num_snps)
    Sigma_beta_list <- vector("list", num_snps)



    for (m in 1:num_snps) {
        # 提取对应的2x2协方差矩阵
        Sigma_gamma_inv_list[[m]] <- beta_sigma_exp[(2 * m - 1):(2 * m), ]
        Sigma_beta_inv_list[[m]] <- beta_sigma_out[(2 * m - 1):(2 * m), ]

        # 计算逆矩阵，添加数值稳定性检查
        tryCatch(
            {
                Sigma_gamma_list[[m]] <- solve(Sigma_gamma_inv_list[[m]])
                Sigma_beta_list[[m]] <- solve(Sigma_beta_inv_list[[m]])
            },
            error = function(e) {
                # 如果矩阵奇异，尝试添加小的正则化项
                warning(paste("Singular matrix encountered at SNP", m, ". Adding regularization."))
                reg_term <- 1e-8 * diag(2)
                Sigma_gamma_list[[m]] <- solve(Sigma_gamma_inv_list[[m]] + reg_term)
                Sigma_beta_list[[m]] <- solve(Sigma_beta_inv_list[[m]] + reg_term)
            }
        )
    }

    preprocessing_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))


    # --- 2. 调用C++核心迭代函数 ---

    core_start_time <- Sys.time()

    tryCatch(
        {
            result <- iterative_mle_core_with_history(
                gamma_hat = gamma_hat,
                beta_hat = beta_hat,
                Sigma_gamma_list = Sigma_gamma_list,
                Sigma_beta_list = Sigma_beta_list,
                Sigma_gamma_inv_list = Sigma_gamma_inv_list,
                Sigma_beta_inv_list = Sigma_beta_inv_list,
                a = a, b = b,
                max_iter = max_iter,
                tol = tol,
                alpha_init = alpha_init
            )
        },
        error = function(e) {
            stop(paste("Error in C++ core function:", e$message))
        }
    )

    core_time <- as.numeric(difftime(Sys.time(), core_start_time, units = "secs"))



    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))


    final_sets <- list(
        A = as.integer(result$set_A_indices),
        B = as.integer(result$set_B_indices)
    )

    # 输出关键结果


    return(list(
        alpha_final = as.numeric(result$alpha_final),
        gamma_final = result$gamma_final,
        iterations = result$iterations,
        converged = result$converged,
        sets = final_sets,
        alpha_history = as.numeric(result$alpha_history),
        set_A_history = result$set_A_history,
        set_B_history = result$set_B_history,
        # 返回预处理的数据，便于进一步分析
        data = list(
            beta_hat_exp = beta_hat_exp,
            beta_hat_out = beta_hat_out,
            Sigma_gamma_list = Sigma_gamma_list,
            Sigma_beta_list = Sigma_beta_list,
            Sigma_gamma_inv_list = Sigma_gamma_inv_list,
            Sigma_beta_inv_list = Sigma_beta_inv_list
        ),
        # 添加计算时间信息
        timing = list(
            preprocessing_time = preprocessing_time,
            core_time = core_time,
            total_time = total_time
        )
    ))
}

calculate_alpha_variance_ivw_2nd_order <- function(estimation_result) {
    # --- 1. 提取所需变量 ---
    alpha <- estimation_result$alpha_final
    gamma_final <- t(estimation_result$gamma_final)
    sets <- estimation_result$sets

    # 注意：二阶近似需要gamma的方差矩阵
    Sigma_beta_list <- estimation_result$data$Sigma_beta_list
    Sigma_gamma_list <- estimation_result$data$Sigma_gamma_list

    all_selected_snps <- sort(union(sets$A, sets$B))
    K <- length(all_selected_snps)

    if (K == 0) {
        warning("没有被选中的SNPs，无法计算方差。")
        return(list(variance = NA, std_error = NA))
    }

    # --- 2. 使用二阶IVW求和策略计算V* ---
    V_star_total <- 0

    for (m in all_selected_snps) {
        gamma_m <- gamma_final[, m]

        if (m %in% intersect(sets$A, sets$B)) {
            # --- 双变量情况 (A ∩ B) ---
            Sigma_beta_m <- Sigma_beta_list[[m]]
            Sigma_gamma_m <- Sigma_gamma_list[[m]]

            # 计算有效方差矩阵
            effective_variance_matrix <- Sigma_beta_m + alpha^2 * Sigma_gamma_m

            # 求逆并计算二次型 (信息量)
            # 使用 tryCatch 避免因矩阵不可逆而中断
            inv_effective_var <- tryCatch(solve(effective_variance_matrix), error = function(e) NULL)

            if (is.null(inv_effective_var)) {
                V_m <- 0 # 如果不可逆，则此SNP不提供信息
                warning(paste("SNP", m, "的有效方差矩阵不可逆。"))
            } else {
                V_m <- t(gamma_m) %*% inv_effective_var %*% gamma_m
            }
        } else if (m %in% setdiff(sets$A, sets$B)) {
            # --- 单变量情况 (A \ B) ---
            var_beta_o_m <- Sigma_beta_list[[m]][1, 1]
            var_gamma_o_m <- Sigma_gamma_list[[m]][1, 1]

            effective_variance <- var_beta_o_m + alpha^2 * var_gamma_o_m

            if (effective_variance > 0) {
                V_m <- gamma_m[1]^2 / effective_variance
            } else {
                V_m <- 0
            }
        } else if (m %in% setdiff(sets$B, sets$A)) {
            # --- 单变量情况 (B \ A) ---
            var_beta_p_m <- Sigma_beta_list[[m]][2, 2]
            var_gamma_p_m <- Sigma_gamma_list[[m]][2, 2]

            effective_variance <- var_beta_p_m + alpha^2 * var_gamma_p_m

            if (effective_variance > 0) {
                V_m <- gamma_m[2]^2 / effective_variance
            } else {
                V_m <- 0
            }
        }

        # 累加每个SNP提供的信息量
        if (!is.na(V_m) && V_m > 0) {
            V_star_total <- V_star_total + V_m
        }
    }

    # --- 3. 计算最终方差和标准误 ---
    if (V_star_total <= 1e-8) { # 使用一个小的阈值以增加数值稳定性
        warning("计算出的总信息量(V*)为非正数或接近零，无法计算方差。")
        return(list(variance = NA, std_error = NA))
    }

    variance_alpha <- 1 / V_star_total
    se_alpha <- sqrt(variance_alpha)

    return(list(variance = variance_alpha, std_error = se_alpha))
}

calculate_alpha_variance_final <- function(estimation_result) {
    # --- 1. 从结果对象中优雅地提取所有所需变量 ---
    alpha <- estimation_result$alpha_final
    gamma <- t(estimation_result$gamma_final)
    sets <- estimation_result$sets

    # 从嵌套的 data 列表中提取
    beta_hat <- t(estimation_result$data$beta_hat_out)
    Sigma_gamma_list <- estimation_result$data$Sigma_gamma_list
    Sigma_beta_list <- estimation_result$data$Sigma_beta_list
    Sigma_gamma_inv_list <- estimation_result$data$Sigma_gamma_inv_list
    Sigma_beta_inv_list <- estimation_result$data$Sigma_beta_inv_list

    all_selected_snps <- sort(union(sets$A, sets$B))
    K <- length(all_selected_snps)

    if (K == 0) {
        warning("没有被选中的SNPs，无法计算方差。")
        return(list(variance = NA, std_error = NA, hessian_eigenvalues = NA))
    }

    total_dim <- 1 + 2 * K
    H <- matrix(0, nrow = total_dim, ncol = total_dim)

    # --- 2. 填充矩阵 H ---
    H_alpha_alpha_total <- 0

    for (i in 1:K) {
        m <- all_selected_snps[i]
        row_start <- 2 + 2 * (i - 1)

        # H_αα 和 H_γm,α 的计算不变
        if (m %in% intersect(sets$A, sets$B)) {
            H_alpha_alpha_total <- H_alpha_alpha_total +
                t(gamma[, m]) %*% Sigma_beta_inv_list[[m]] %*% gamma[, m]
            H_gmgm <- alpha^2 * Sigma_beta_inv_list[[m]] +
                Sigma_gamma_inv_list[[m]]
            H_gma <- -Sigma_beta_inv_list[[m]] %*% (beta_hat[, m] - 2 * alpha * gamma[, m])
        } else if (m %in% setdiff(sets$A, sets$B)) {
            var_beta_o_m <- Sigma_beta_list[[m]][1, 1]
            var_gamma_o_m <- Sigma_gamma_list[[m]][1, 1] # 获取gamma_o的方差

            H_alpha_alpha_total <- H_alpha_alpha_total + gamma[1, m]^2 / var_beta_o_m
            # ✨ 核心修正：根据最终似然函数更新 H_gmgm
            H_gmgm <- matrix(c(alpha^2 / var_beta_o_m + 1 / var_gamma_o_m, 0, 0, 0), nrow = 2)
            H_gma <- c((-beta_hat[1, m] + 2 * alpha * gamma[1, m]) / var_beta_o_m, 0)
        } else if (m %in% setdiff(sets$B, sets$A)) {
            var_beta_p_m <- Sigma_beta_list[[m]][2, 2]
            var_gamma_p_m <- Sigma_gamma_list[[m]][2, 2] # 获取gamma_p的方差

            H_alpha_alpha_total <- H_alpha_alpha_total + gamma[2, m]^2 / var_beta_p_m
            # ✨ 核心修正：根据最终似然函数更新 H_gmgm
            H_gmgm <- matrix(c(0, 0, 0, alpha^2 / var_beta_p_m + 1 / var_gamma_p_m), nrow = 2)
            H_gma <- c(0, (-beta_hat[2, m] + 2 * alpha * gamma[2, m]) / var_beta_p_m)
        }

        # 填充矩阵
        H[row_start:(row_start + 1), row_start:(row_start + 1)] <- H_gmgm
        H[row_start:(row_start + 1), 1] <- H_gma
        H[1, row_start:(row_start + 1)] <- t(H_gma)
    }

    H[1, 1] <- H_alpha_alpha_total
    H <- H[rowSums(H != 0) > 0, colSums(H != 0) > 0]
    # --- 3. 求逆并提取结果（带正则化） ---
    eigenvals <- eigen(H, only.values = TRUE)$values
    min_eigenval <- min(eigenvals)
    regularization_applied <- FALSE

    # 检查特征值并应用正则化
    if (any(eigenvals <= 1e-8)) { # 使用更严格的阈值
        warning("Hessian 矩阵非正定或接近奇异。正在应用正则化。")

        # 计算需要添加的正则化项
        reg_value <- max(1e-6, abs(min_eigenval) + 1e-6)

        # 应用正则化：添加到对角线
        H_regularized <- H + reg_value * diag(nrow(H))

        # 验证正则化后的矩阵
        eigenvals_reg <- eigen(H_regularized, only.values = TRUE)$values

        if (any(eigenvals_reg <= 0)) {
            # 如果还是有问题，使用更强的正则化
            reg_value <- max(1e-4, abs(min(eigenvals_reg)) + 1e-4)
            H_regularized <- H + reg_value * diag(nrow(H))
            warning(paste("应用强正则化，正则化值:", reg_value))
        }

        H <- H_regularized
    }
    H_inv <- tryCatch(solve(H), error = function(e) NULL)

    if (is.null(H_inv)) {
        return(list(variance = NA, std_error = NA, hessian_eigenvalues = eigenvals))
    }

    var_alpha <- H_inv[1, 1]

    if (var_alpha < 0) {
        warning("计算出的 Alpha 方差为负数。")
        return(list(variance = var_alpha, std_error = NA, hessian_eigenvalues = eigenvals))
    }

    se_alpha <- sqrt(var_alpha)

    return(list(variance = var_alpha, std_error = se_alpha, hessian_eigenvalues = eigenvals))
}
weighted_estimation_robust <- function(data_table) {
    # 调试信息

    # 确保数据是data.frame格式
    if (!is.data.frame(data_table)) {
        cat("转换为data.frame...\n")
        data_table <- as.data.frame(data_table)
    }

    # 检查必要的列是否存在
    required_cols <- c("alpha", "alpha_se", "bic")
    missing_cols <- required_cols[!required_cols %in% names(data_table)]
    if (length(missing_cols) > 0) {
        stop("缺少必要的列: ", paste(missing_cols, collapse = ", "))
    }



    # 手动检查NA值
    alpha_na <- is.na(data_table$alpha)
    alpha_se_na <- is.na(data_table$alpha_se)
    bic_na <- is.na(data_table$bic)

    # 找出完整的行
    valid_rows <- !alpha_na & !alpha_se_na & !bic_na



    if (sum(valid_rows) == 0) {
        warning("没有完整的数据行用于计算")
        return(list(
            weighted_alpha = NA,
            weighted_variance = NA,
            weighted_se = NA,
            n_valid = 0
        ))
    }

    # 提取有效数据
    alpha <- data_table$alpha[valid_rows]
    alpha_se <- data_table$alpha_se[valid_rows]
    bic <- data_table$bic[valid_rows]

    # 确保数据是数值型
    alpha <- as.numeric(alpha)
    alpha_se <- as.numeric(alpha_se)
    bic <- as.numeric(bic)

    # 检查是否还有NA值
    if (any(is.na(alpha)) || any(is.na(alpha_se)) || any(is.na(bic))) {
        stop("数据转换后仍有NA值")
    }


    # 计算基于BIC的权重
    min_bic <- min(bic)
    weights <- exp(-0.5 * (bic - min_bic))
    weights <- weights / sum(weights)



    # 计算加权点估计
    weighted_alpha <- sum(weights * alpha)

    # 计算加权标准误
    se_components <- weights * sqrt(alpha_se^2 + (alpha - weighted_alpha)^2)
    weighted_se <- sum(se_components)

    # 加权方差
    weighted_variance <- weighted_se^2

    # 返回结果
    return(list(
        weighted_alpha = weighted_alpha,
        weighted_variance = weighted_variance,
        weighted_se = weighted_se,
        weights = weights,
        n_valid = sum(valid_rows),
        valid_indices = which(valid_rows)
    ))
}
calculate_bic_cml <- function(input_model, a, b, all_snps, n) {
    # 获取 SNP 的总数 m
    m <- length(all_snps)
    A_hat <- input_model$sets$A
    B_hat <- input_model$sets$B
    alpha <- input_model$alpha_final
    gamma_model <- input_model$gamma_final
    beta_exp_sigma_inv <- input_model$data$Sigma_gamma_inv_list
    beta_out_sigma_inv <- input_model$data$Sigma_beta_inv_list
    beta_exp_hat <- input_model$data$beta_hat_exp
    beta_out_hat <- input_model$data$beta_hat_out
    # 1. 确定各个 SNP 集合的索引
    # -------------------------------------
    snps_intersect <- intersect(A_hat, B_hat)
    snps_A_only <- setdiff(A_hat, B_hat)
    snps_B_only <- setdiff(B_hat, A_hat)
    snps_union <- union(A_hat, B_hat)
    snps_other <- setdiff(all_snps, snps_union)

    # 将 SNP ID 转换为在矩阵/列表中的行索引
    idx_intersect <- match(snps_intersect, all_snps)
    idx_A_only <- match(snps_A_only, all_snps)
    idx_B_only <- match(snps_B_only, all_snps)
    idx_other <- match(snps_other, all_snps)

    # 初始化各部分计算结果
    term1 <- 0
    term2 <- 0
    term3 <- 0
    term4 <- 0

    # 2. 计算第一项 (m ∈ Â ∩ B̂)
    # -------------------------------------
    if (length(idx_intersect) > 0) {
        term1_values <- sapply(idx_intersect, function(i) {
            # 提取对应 SNP 的数据
            beta_hat_i <- beta_out_hat[i, ]
            gamma_hat_i <- beta_exp_hat[i, ]
            gamma_mod_i <- gamma_model[i, ]
            Sigma_beta_inv_i <- beta_out_sigma_inv[[i]]
            Sigma_gamma_inv_i <- beta_exp_sigma_inv[[i]]

            # 计算差值向量
            diff_beta <- beta_hat_i - alpha * gamma_mod_i
            diff_gamma <- gamma_hat_i - gamma_mod_i

            # 计算二次型 (t(x) %*% M %*% x)
            val_beta <- t(diff_beta) %*% Sigma_beta_inv_i %*% diff_beta
            val_gamma <- t(diff_gamma) %*% Sigma_gamma_inv_i %*% diff_gamma

            return(val_beta + val_gamma)
        })
        term1 <- sum(term1_values)
    }

    # 3. 计算第二项 (m ∈ Â \ B̂)
    # -------------------------------------
    if (length(idx_A_only) > 0) {
        term2_values <- sapply(idx_A_only, function(i) {
            # 提取 'o' 分量 (第一列)
            beta_o_hat <- beta_out_hat[i, 1]
            gamma_o_hat <- beta_exp_hat[i, 1]
            gamma_o_mod <- gamma_model[i, 1]

            # 从逆协方差矩阵计算方差 σ²
            # σ²_o = (Σ⁻¹)_pp / det(Σ⁻¹)
            Sigma_beta_inv_i <- beta_out_sigma_inv[[i]]
            Sigma_gamma_inv_i <- beta_exp_sigma_inv[[i]]

            sigma2_beta_o <- Sigma_beta_inv_i[2, 2] / det(Sigma_beta_inv_i)
            sigma2_gamma_o <- Sigma_gamma_inv_i[2, 2] / det(Sigma_gamma_inv_i)

            # 计算该项的值
            val_beta <- (beta_o_hat - alpha * gamma_o_mod)^2 / sigma2_beta_o
            val_gamma <- (gamma_o_hat - gamma_o_mod)^2 / sigma2_gamma_o

            return(val_beta + val_gamma)
        })
        term2 <- sum(term2_values)
    }

    # 4. 计算第三项 (m ∈ B̂ \ Â)
    # -------------------------------------
    if (length(idx_B_only) > 0) {
        term3_values <- sapply(idx_B_only, function(i) {
            # 提取 'p' 分量 (第二列)
            beta_p_hat <- beta_out_hat[i, 2]
            gamma_p_hat <- beta_exp_hat[i, 2]
            gamma_p_mod <- gamma_model[i, 2]

            # 从逆协方差矩阵计算方差 σ²
            # σ²_p = (Σ⁻¹)_oo / det(Σ⁻¹)
            Sigma_beta_inv_i <- beta_out_sigma_inv[[i]]
            Sigma_gamma_inv_i <- beta_exp_sigma_inv[[i]]

            sigma2_beta_p <- Sigma_beta_inv_i[1, 1] / det(Sigma_beta_inv_i)
            sigma2_gamma_p <- Sigma_gamma_inv_i[1, 1] / det(Sigma_gamma_inv_i)

            # 计算该项的值
            val_beta <- (beta_p_hat - alpha * gamma_p_mod)^2 / sigma2_beta_p
            val_gamma <- (gamma_p_hat - gamma_p_mod)^2 / sigma2_gamma_p

            return(val_beta + val_gamma)
        })
        term3 <- sum(term3_values)
    }

    # 5. 计算第四项 (其他 m)
    # -------------------------------------
    if (length(idx_other) > 0) {
        term4_values <- sapply(idx_other, function(i) {
            gamma_hat_i <- beta_exp_hat[i, ]
            gamma_mod_i <- gamma_model[i, ]
            Sigma_gamma_inv_i <- beta_exp_sigma_inv[[i]]

            diff_gamma <- gamma_hat_i - gamma_mod_i
            val_gamma <- t(diff_gamma) %*% Sigma_gamma_inv_i %*% diff_gamma

            return(val_gamma)
        })
        term4 <- sum(term4_values)
    }

    # 6. 计算第五项 (惩罚项)
    # -------------------------------------
    penalty <- 1 / 2 * log(n) * (m - a + m - b)

    # 7. 加总所有项得到最终 BIC 值
    # -------------------------------------
    bic_total <- term1 + term2 + term3 + term4 + penalty

    return(bic_total)
}
cml_family_ultral <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, n = 1000, cov_index = "var") {
    all_snps <- 1:nrow(beta_hat_exp)
    data_table <- expand.grid(a = all_snps, b = all_snps)
    data_table$alpha <- NA
    data_table$alpha_se <- NA
    data_table$bic <- NA
    for (i in 1:nrow(data_table)) {
        data_any <- run_iterative_mle_v2_rcpp(
            beta_hat_exp, beta_hat_out, beta_sigma_exp,
            beta_sigma_out,
            a = data_table[i, 1],
            b = data_table[i, 2]
        )

        data_table$alpha[i] <- data_any$alpha_final
        if (cov_index == "normal") {
            alpha_se <- calculate_alpha_variance_final(data_any)
        } else {
            alpha_se <- calculate_alpha_variance_ivw_2nd_order(data_any)
        }

        data_table$alpha_se[i] <- alpha_se$std_error
        data_table$bic[i] <- calculate_bic_cml(data_any,
            a = data_table[i, 1], b = data_table[i, 2],
            all_snps = all_snps, n = n
        )
    }
    result <- weighted_estimation_robust(data_table)

    return(list(data_table, result, data_any))
}
cml_family_ultral_gamma <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, n = 1000) {
    all_snps <- 1:nrow(beta_hat_exp)
    data_table <- expand.grid(a = all_snps, b = all_snps)
    data_table$alpha <- 0
    data_table$alpha_se <- 0
    data_table$bic <- 0
    for (i in 1:nrow(data_table)) {
        data_any <- run_iterative_mle_v2_rcpp_gamma(
            beta_hat_exp, beta_hat_out, beta_sigma_exp,
            beta_sigma_out,
            a = data_table[i, 1],
            b = data_table[i, 2]
        )
        data_table$alpha[i] <- data_any$alpha_final
        alpha_se <- calculate_alpha_variance_final(data_any)
        data_table$alpha_se[i] <- alpha_se$std_error
        data_table$bic[i] <- calculate_bic_cml(data_any,
            a = data_table[i, 1], b = data_table[i, 2],
            all_snps = all_snps, n = n
        )
    }
    result <- weighted_estimation_robust(data_table)

    return(list(data_table, result, data_any))
}

# %% 多起点优化函数

run_iterative_mle_multistart <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, a, b,
    max_iter = 100, tol = 1e-6,
    num_starts = 5, alpha_init_range = c(-0.0001, 0.0001),
    consensus_method = "majority", consensus_threshold = 0.6,
    seed = NULL) {
    # 设置随机种子以保证可重复性
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # 输入验证（保持原有验证逻辑）
    if (!is.matrix(beta_hat_exp) || !is.matrix(beta_hat_out)) {
        stop("beta_hat_exp and beta_hat_out must be matrices")
    }

    if (nrow(beta_hat_exp) != nrow(beta_hat_out)) {
        stop("beta_hat_exp and beta_hat_out must have the same number of rows (SNPs)")
    }

    if (ncol(beta_hat_exp) != 2 || ncol(beta_hat_out) != 2) {
        stop("beta_hat_exp and beta_hat_out must have exactly 2 columns")
    }

    num_snps <- nrow(beta_hat_exp)

    if (a <= 0 || b <= 0) {
        stop("Parameters a and b must be positive integers")
    }

    if (max_iter <= 0) {
        stop("max_iter must be a positive integer")
    }

    if (tol <= 0) {
        stop("tol must be a positive number")
    }

    # 新增参数验证
    if (num_starts <= 0) {
        stop("num_starts must be a positive integer")
    }

    if (length(alpha_init_range) != 2 || alpha_init_range[1] >= alpha_init_range[2]) {
        stop("alpha_init_range must be a vector of length 2 with first element < second element")
    }

    if (!consensus_method %in% c("majority", "best_likelihood")) {
        stop("consensus_method must be one of: 'majority', 'best_likelihood'")
    }

    if (consensus_threshold <= 0 || consensus_threshold > 1) {
        stop("consensus_threshold must be between 0 and 1")
    }

    # --- 1. 数据预处理 (R语言部分，保持不变) ---
    start_time <- Sys.time()

    gamma_hat <- t(beta_hat_exp)
    beta_hat <- t(beta_hat_out)

    Sigma_gamma_inv_list <- vector("list", num_snps)
    Sigma_beta_inv_list <- vector("list", num_snps)
    Sigma_gamma_list <- vector("list", num_snps)
    Sigma_beta_list <- vector("list", num_snps)

    for (m in 1:num_snps) {
        # 提取对应的2x2协方差矩阵
        Sigma_gamma_inv_list[[m]] <- beta_sigma_exp[(2 * m - 1):(2 * m), ]
        Sigma_beta_inv_list[[m]] <- beta_sigma_out[(2 * m - 1):(2 * m), ]

        # 计算逆矩阵，添加数值稳定性检查
        tryCatch(
            {
                Sigma_gamma_list[[m]] <- solve(Sigma_gamma_inv_list[[m]])
                Sigma_beta_list[[m]] <- solve(Sigma_beta_inv_list[[m]])
            },
            error = function(e) {
                # 如果矩阵奇异，尝试添加小的正则化项
                warning(paste("Singular matrix encountered at SNP", m, ". Adding regularization."))
                reg_term <- 1e-8 * diag(2)
                Sigma_gamma_list[[m]] <- solve(Sigma_gamma_inv_list[[m]] + reg_term)
                Sigma_beta_list[[m]] <- solve(Sigma_beta_inv_list[[m]] + reg_term)
            }
        )
    }

    preprocessing_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    # --- 2. 多起点优化 ---
    core_start_time <- Sys.time()

    # 生成多个初始值
    alpha_inits <- runif(num_starts, alpha_init_range[1], alpha_init_range[2])

    # 存储所有起点的结果
    multistart_results <- vector("list", num_starts)
    convergence_status <- logical(num_starts)
    likelihood_values <- numeric(num_starts)

    cat("Running multi-start optimization with", num_starts, "starting points...\n")

    for (i in 1:num_starts) {
        cat("Starting point", i, "/ ", num_starts, "(alpha_init =", round(alpha_inits[i], 4), ")...\n")

        tryCatch(
            {
                result <- iterative_mle_core(
                    gamma_hat = gamma_hat,
                    beta_hat = beta_hat,
                    Sigma_gamma_list = Sigma_gamma_list,
                    Sigma_beta_list = Sigma_beta_list,
                    Sigma_gamma_inv_list = Sigma_gamma_inv_list,
                    Sigma_beta_inv_list = Sigma_beta_inv_list,
                    a = a, b = b,
                    max_iter = max_iter,
                    tol = tol,
                    alpha_init = alpha_inits[i]
                )

                multistart_results[[i]] <- result
                convergence_status[i] <- result$converged

                # 计算似然值（如果C++函数提供的话，否则需要单独计算）
                # 这里假设C++函数返回final_likelihood，如果没有需要另外计算
                if (!is.null(result$final_likelihood)) {
                    likelihood_values[i] <- result$final_likelihood
                } else {
                    likelihood_values[i] <- NA # 需要单独计算似然值
                }
            },
            error = function(e) {
                warning(paste("Starting point", i, "failed:", e$message))
                multistart_results[[i]] <- NULL
                convergence_status[i] <- FALSE
                likelihood_values[i] <- -Inf
            }
        )
    }

    core_time <- as.numeric(difftime(Sys.time(), core_start_time, units = "secs"))

    # --- 3. 共识分析 ---
    consensus_start_time <- Sys.time()

    # 筛选成功收敛的结果
    converged_results <- multistart_results[convergence_status]
    converged_likelihoods <- likelihood_values[convergence_status]

    if (length(converged_results) == 0) {
        return(NA)
    }

    cat("Successfully converged:", length(converged_results), "out of", num_starts, "starting points\n")

    # 获取共识结果
    consensus_result <- get_consensus_result(
        converged_results,
        converged_likelihoods,
        consensus_method,
        consensus_threshold,
        num_snps
    )

    consensus_time <- as.numeric(difftime(Sys.time(), consensus_start_time, units = "secs"))
    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

    # --- 4. 构建最终输出（保持原有格式） ---
    final_sets <- list(
        A = as.integer(consensus_result$set_A_indices),
        B = as.integer(consensus_result$set_B_indices)
    )

    # 返回结果，保持原有接口不变，并添加多起点信息
    return(list(
        # 原有输出格式（保持不变）
        alpha_final = as.numeric(consensus_result$alpha_final),
        gamma_final = consensus_result$gamma_final,
        iterations = consensus_result$iterations,
        converged = consensus_result$converged,
        sets = final_sets,
        alpha_history = as.numeric(consensus_result$alpha_history),

        # 原有数据和时间信息
        data = list(
            beta_hat_exp = beta_hat_exp,
            beta_hat_out = beta_hat_out,
            Sigma_gamma_list = Sigma_gamma_list,
            Sigma_beta_list = Sigma_beta_list,
            Sigma_gamma_inv_list = Sigma_gamma_inv_list,
            Sigma_beta_inv_list = Sigma_beta_inv_list
        ),
        timing = list(
            preprocessing_time = preprocessing_time,
            core_time = core_time,
            consensus_time = consensus_time,
            total_time = total_time
        ),

        # 新增：多起点优化信息
        multistart_info = list(
            num_starts = num_starts,
            num_converged = length(converged_results),
            convergence_rate = length(converged_results) / num_starts,
            alpha_inits = alpha_inits,
            convergence_status = convergence_status,
            likelihood_values = likelihood_values,
            consensus_method = consensus_method,
            consensus_threshold = consensus_threshold,
            consensus_confidence = consensus_result$confidence
        )
    ))
}

# 辅助函数：获取共识结果
get_consensus_result <- function(results, likelihoods, method, threshold, num_snps) {
    if (method == "best_likelihood") {
        # 选择最佳似然值的结果
        if (all(is.na(likelihoods))) {
            # 如果没有似然值，选择第一个结果
            best_idx <- 1
        } else {
            best_idx <- which.max(likelihoods)
        }

        consensus <- results[[best_idx]]
        consensus$confidence <- 1.0 # 单点选择，置信度为1

        return(consensus)
    } else if (method == "majority") {
        # 基于集合成员的多数投票
        return(get_majority_consensus(results, threshold, num_snps))
    }
}

# 多数投票共识
get_majority_consensus <- function(results, threshold, num_snps) {
    n_results <- length(results)

    # 统计每个SNP在集合A和B中出现的频率
    setA_votes <- integer(num_snps)
    setB_votes <- integer(num_snps)

    for (result in results) {
        if (!is.null(result$set_A_indices)) {
            setA_votes[result$set_A_indices] <- setA_votes[result$set_A_indices] + 1
        }
        if (!is.null(result$set_B_indices)) {
            setB_votes[result$set_B_indices] <- setB_votes[result$set_B_indices] + 1
        }
    }

    # 基于阈值确定共识集合
    consensus_setA <- which(setA_votes / n_results >= threshold)
    consensus_setB <- which(setB_votes / n_results >= threshold)

    # 计算共识置信度
    max_vote_rate_A <- if (length(setA_votes) > 0) max(setA_votes) / n_results else 0
    max_vote_rate_B <- if (length(setB_votes) > 0) max(setB_votes) / n_results else 0
    confidence <- max(max_vote_rate_A, max_vote_rate_B)

    # 检查是否满足共识阈值
    if (confidence < threshold) {
        # 不满足阈值，返回NA结果
        return(list(
            alpha_final = NA_real_,
            gamma_final = matrix(NA_real_,
                nrow = nrow(results[[1]]$gamma_final),
                ncol = ncol(results[[1]]$gamma_final)
            ),
            iterations = NA_integer_,
            converged = FALSE,
            set_A_indices = integer(0),
            set_B_indices = integer(0),
            alpha_history = NA_real_,
            confidence = confidence
        ))
    }

    # 满足阈值，从满足共识的结果中随机选择一个
    # 找出包含共识集合的结果
    valid_results <- list()
    for (i in seq_along(results)) {
        result <- results[[i]]
        # 检查该结果是否包含共识集合
        contains_consensus_A <- all(consensus_setA %in% result$set_A_indices)
        contains_consensus_B <- all(consensus_setB %in% result$set_B_indices)

        if (contains_consensus_A && contains_consensus_B) {
            valid_results <- append(valid_results, list(result), length(valid_results))
        }
    }

    # 如果没有完全包含共识集合的结果，则从所有结果中随机选择
    if (length(valid_results) == 0) {
        valid_results <- results
    }

    # 随机选择一个结果
    selected_idx <- sample(length(valid_results), 1)
    selected_result <- valid_results[[selected_idx]]

    # 但是集合要使用共识结果
    selected_result$set_A_indices <- consensus_setA
    selected_result$set_B_indices <- consensus_setB
    selected_result$confidence <- confidence

    return(selected_result)
}



# 为了保持向后兼容，提供原函数的别名

cml_family_multistart <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, n = 1000, cov_index = "var",
    # 新增多起点优化参数
    num_starts = 5, alpha_init_range = c(-0.001, 0.001),
    consensus_method = "majority", consensus_threshold = 0.6,
    seed = NULL,
    # 控制参数
    max_iter = 100, tol = 1e-6,
    verbose = FALSE) {
    # 输入验证
    if (!cov_index %in% c("var", "normal")) {
        stop("cov_index must be either 'var' or 'normal'")
    }

    if (!consensus_method %in% c("majority", "best_likelihood")) {
        stop("consensus_method must be one of: 'majority', 'best_likelihood'")
    }

    all_snps <- 1:nrow(beta_hat_exp)
    data_table <- expand.grid(a = all_snps, b = all_snps)
    data_table$alpha <- NA
    data_table$alpha_se <- NA
    data_table$bic <- NA

    # 新增列记录多起点优化信息
    data_table$converged <- NA
    data_table$consensus_confidence <- NA
    data_table$num_converged_starts <- NA
    data_table$convergence_rate <- NA

    total_combinations <- nrow(data_table)
    failed_combinations <- 0
    consensus_failures <- 0

    if (verbose) {
        cat("Starting CML family analysis with", total_combinations, "combinations...\n")
        cat(
            "Multi-start parameters: num_starts =", num_starts,
            ", consensus_method =", consensus_method,
            ", threshold =", consensus_threshold, "\n\n"
        )
    }

    # 记录最后一个成功的结果用于兼容性
    last_successful_result <- NULL

    for (i in 1:nrow(data_table)) {
        if (verbose && i %% 10 == 0) {
            cat(
                "Processing combination", i, "/", total_combinations,
                "(a =", data_table[i, 1], ", b =", data_table[i, 2], ")\n"
            )
        }

        tryCatch(
            {
                # 使用多起点优化
                data_any <- run_iterative_mle_multistart(
                    beta_hat_exp, beta_hat_out, beta_sigma_exp,
                    beta_sigma_out,
                    a = data_table[i, 1],
                    b = data_table[i, 2],
                    num_starts = num_starts,
                    alpha_init_range = alpha_init_range,
                    consensus_method = consensus_method,
                    consensus_threshold = consensus_threshold,
                    seed = seed,
                    max_iter = max_iter,
                    tol = tol
                )

                # 记录多起点优化信息
                data_table$converged[i] <- data_any$converged
                data_table$consensus_confidence[i] <- data_any$multistart_info$consensus_confidence
                data_table$num_converged_starts[i] <- data_any$multistart_info$num_converged
                data_table$convergence_rate[i] <- data_any$multistart_info$convergence_rate

                # 检查是否达成共识
                if (data_any$converged &&
                    data_any$multistart_info$consensus_confidence >= consensus_threshold) {
                    # 达成共识，记录结果
                    data_table$alpha[i] <- data_any$alpha_final

                    # 计算标准误
                    if (cov_index == "normal") {
                        alpha_se <- calculate_alpha_variance_final(data_any)
                    } else {
                        alpha_se <- calculate_alpha_variance_ivw_2nd_order(data_any)
                    }

                    data_table$alpha_se[i] <- alpha_se$std_error
                    data_table$bic[i] <- calculate_bic_cml(data_any,
                        a = data_table[i, 1], b = data_table[i, 2],
                        all_snps = all_snps, n = n
                    )

                    # 记录成功结果
                    last_successful_result <- data_any
                } else {
                    # 未达成共识，记录NA
                    data_table$alpha[i] <- NA
                    data_table$alpha_se[i] <- NA
                    data_table$bic[i] <- NA
                    consensus_failures <- consensus_failures + 1

                    if (verbose) {
                        cat(
                            "  Warning: No consensus reached for combination (",
                            data_table[i, 1], ",", data_table[i, 2],
                            "), confidence =",
                            round(data_any$multistart_info$consensus_confidence, 3), "\n"
                        )
                    }
                }
            },
            error = function(e) {
                # 处理错误情况
                data_table$alpha[i] <<- NA
                data_table$alpha_se[i] <<- NA
                data_table$bic[i] <<- NA
                data_table$converged[i] <<- FALSE
                data_table$consensus_confidence[i] <<- 0
                data_table$num_converged_starts[i] <<- 0
                data_table$convergence_rate[i] <<- 0

                failed_combinations <<- failed_combinations + 1

                if (verbose) {
                    cat(
                        "  Error in combination (", data_table[i, 1], ",",
                        data_table[i, 2], "):", e$message, "\n"
                    )
                }
            }
        )
    }

    # 输出汇总信息
    successful_combinations <- sum(!is.na(data_table$alpha))

    if (verbose) {
        cat("\n=== Analysis Summary ===\n")
        cat("Total combinations:", total_combinations, "\n")
        cat("Successful combinations:", successful_combinations, "\n")
        cat("Failed combinations:", failed_combinations, "\n")
        cat("Consensus failures:", consensus_failures, "\n")
        cat("Success rate:", round(successful_combinations / total_combinations * 100, 1), "%\n")

        if (successful_combinations > 0) {
            avg_confidence <- mean(data_table$consensus_confidence, na.rm = TRUE)
            avg_convergence_rate <- mean(data_table$convergence_rate, na.rm = TRUE)
            cat("Average consensus confidence:", round(avg_confidence, 3), "\n")
            cat("Average convergence rate:", round(avg_convergence_rate, 3), "\n")
        }
    }

    # 进行加权估计（只使用成功的结果）
    valid_data <- data_table[!is.na(data_table$alpha), ]

    if (nrow(valid_data) == 0) {
        warning("No valid results obtained. Cannot perform weighted estimation.")
        result <- list(
            estimate = NA,
            se = NA,
            ci_lower = NA,
            ci_upper = NA,
            p_value = NA,
            n_valid = 0
        )
    } else {
        result <- weighted_estimation_robust(valid_data)
    }

    # 构建多起点优化的汇总信息
    multistart_summary <- list(
        total_combinations = total_combinations,
        successful_combinations = successful_combinations,
        failed_combinations = failed_combinations,
        consensus_failures = consensus_failures,
        success_rate = successful_combinations / total_combinations,
        average_consensus_confidence = if (successful_combinations > 0) {
            mean(data_table$consensus_confidence, na.rm = TRUE)
        } else {
            NA
        },
        average_convergence_rate = if (successful_combinations > 0) {
            mean(data_table$convergence_rate, na.rm = TRUE)
        } else {
            NA
        },
        parameters = list(
            num_starts = num_starts,
            alpha_init_range = alpha_init_range,
            consensus_method = consensus_method,
            consensus_threshold = consensus_threshold
        )
    )

    return(list(
        data_table = data_table,
        result = result,
        data_any = last_successful_result, # 保持兼容性
        multistart_summary = multistart_summary # 新增汇总信息
    ))
}


# 辅助函数：分析多起点优化结果
analyze_multistart_results <- function(cml_result) {
    if (!"multistart_summary" %in% names(cml_result)) {
        stop("This function requires results from cml_family_multistart")
    }

    data_table <- cml_result$data_table
    summary_info <- cml_result$multistart_summary

    cat("=== Multi-start Optimization Analysis ===\n\n")

    # 基本统计
    cat("Basic Statistics:\n")
    cat("- Total combinations tested:", summary_info$total_combinations, "\n")
    cat("- Successful combinations:", summary_info$successful_combinations, "\n")
    cat("- Success rate:", round(summary_info$success_rate * 100, 1), "%\n")
    cat("- Consensus failures:", summary_info$consensus_failures, "\n\n")

    # 共识质量分析
    if (summary_info$successful_combinations > 0) {
        cat("Consensus Quality:\n")
        cat(
            "- Average consensus confidence:",
            round(summary_info$average_consensus_confidence, 3), "\n"
        )
        cat(
            "- Average convergence rate:",
            round(summary_info$average_convergence_rate, 3), "\n"
        )

        # 置信度分布
        confidence_dist <- table(cut(data_table$consensus_confidence,
            breaks = c(0, 0.5, 0.7, 0.8, 0.9, 1.0),
            include.lowest = TRUE, right = FALSE
        ))
        cat("- Confidence distribution:\n")
        print(confidence_dist)
        cat("\n")
    }

    # 参数设置
    cat("Parameters Used:\n")
    cat("- Number of starting points:", summary_info$parameters$num_starts, "\n")
    cat("- Initial value range:", paste(summary_info$parameters$alpha_init_range, collapse = " to "), "\n")
    cat("- Consensus method:", summary_info$parameters$consensus_method, "\n")
    cat("- Consensus threshold:", summary_info$parameters$consensus_threshold, "\n\n")

    # 推荐
    if (summary_info$success_rate < 0.5) {
        cat("Recommendation: Success rate is low. Consider:\n")
        cat("- Increasing number of starting points\n")
        cat("- Widening initial value range\n")
        cat("- Lowering consensus threshold\n")
    } else if (summary_info$average_consensus_confidence < 0.7) {
        cat("Recommendation: Consensus confidence is low. Consider:\n")
        cat("- Increasing number of starting points\n")
        cat("- Using 'best_likelihood' consensus method\n")
    } else {
        cat("Results look good! High success rate and consensus confidence.\n")
    }
}
# %%

if (FALSE) {
    type_one_1 <- 0
    type_one_2 <- 0
    type_1 <- 0
    for (i in 1:200) {
        phase_two_data_full <- generate_multiple_datasets_v4(
            n = 10, # 要生成的总数据集数量
            num_pleiotropic = 0,
            n_expose_heterogeneity = 0,
            N_exp = 1000, N_out = 1000,
            p_f = 0.3, p_m = 0.3,
            # --- 暴露效应 (当非零时的大小) ---
            beta_FStoOE_exp = 0.3, beta_MStoOE_exp = 0.3,
            beta_OStoOE_exp = 0.3,
            # --- 暴露异质性 ---
            h_beta_FStoOE_exp = 0, h_beta_MStoOE_exp = 0,
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

        # 数据预处理 - 添加错误处理

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

        phase_two_data_analysis_exp <-
            fgwas_for_data_optimized(data_exp,
                processing_func = FMR_trio_optimized_onlyc2
            )
        phase_two_data_analysis_out <-
            fgwas_for_data_optimized(data_out,
                processing_func = FMR_trio_optimized_onlyc2
            )


        beta_sigma_exp <- as.matrix(phase_two_data_analysis_exp$Sigma_inv)
        beta_sigma_out <- as.matrix(phase_two_data_analysis_out$Sigma_inv)
        beta_hat_exp <- as.matrix(phase_two_data_analysis_exp$beta_hat)
        beta_hat_out <- as.matrix(phase_two_data_analysis_out$beta_hat)


        res_b <- cml_family_multistart(
            beta_hat_exp, beta_hat_out, beta_sigma_exp, beta_sigma_out,
            n = 1000, alpha_init_range = c(-0.01, 0.02), cov_index = "normal"
        )
        p_value <- 2 * pnorm(abs(res_b[[2]]$weighted_alpha /
            res_b[[2]]$weighted_se), lower.tail = FALSE)
        p_value

        arrange(res_b[[1]], bic)
        if (p_value < 0.05) {
            type_1 <- type_1 + 1
        }
    }
    type_1 / 200

    arrange(res_b[[1]], bic)
}
