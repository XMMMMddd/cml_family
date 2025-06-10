# %% 加载包
library(Rcpp) # 需要 Rcpp 包来编译和加载 C++ 代码
# %% 加载 C++ 代码
source("cml家庭函数/cml家庭函数overlap版本/source_all_rcpp.R")
# %% 对角阵处理算法
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
                current_block_matrix[1, 4] <- current_off_diag_block_elements[3]
                current_block_matrix[2, 3] <- current_off_diag_block_elements[2]
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
# %% 迭代算法函数封装
perform_iterative_update <- function(
    n_snps, beta_hat_exp, beta_hat_out,
    alpha_input, matrix_big, a_legal, b_legal,
    max_iter = 100, tol = 1e-6) { # 新增 max_iter 和 tol 参数

    # 使用输入参数初始化函数内部变量
    beta_exp <- beta_hat_exp
    alpha <- alpha_input

    iterations <- 0
    converged <- FALSE

    for (iter_count in 1:max_iter) {
        iterations <- iter_count
        alpha_old <- alpha # 保存上一轮的 alpha

        # 更新r
        ## 构造异质性统计量
        t <- matrix(0, ncol = 1, nrow = n_snps)
        f <- matrix(0, ncol = 1, nrow = n_snps)

        for (i in 1:n_snps) {
            matrix_indicator <- 4 * i - 1
            # 计算关于子代的统计量
            beta_hat_exp_i <- beta_hat_exp[i, 1]
            beta_exp_i <- beta_exp[i, 1]
            beta_out_se <- matrix_big[matrix_indicator, 3]
            t[i] <- (beta_hat_exp_i - alpha * beta_exp_i)^2 / beta_out_se

            # 计算关于父母的统计量
            beta_hat_exp_i <- beta_hat_exp[i, 2]
            beta_exp_i <- beta_exp[i, 2]
            beta_out_se <- matrix_big[matrix_indicator + 1, 4]
            f[i] <- (beta_hat_exp_i - alpha * beta_exp_i)^2 / beta_out_se
        }

        a_snps_set <- arrange(data.frame(snps_id = 1:n_snps, statistic = t), t)
        b_snps_set <- arrange(data.frame(snps_id = 1:n_snps, statistic = f), f)

        a_legal_set <- a_snps_set[1:a_legal, 1]
        b_legal_set <- b_snps_set[1:b_legal, 1]
        complete_set <- intersect(a_legal_set, b_legal_set)
        a_reduce_b_set <- setdiff(a_legal_set, b_legal_set)

        # 更新gamma
        beta_exp_new <- beta_exp # 注意：这里 beta_exp_new 应该在循环开始前基于 beta_exp 初始化，或者确保 beta_exp 被正确传递和更新
        ## 更新全数据集中的gamma
        if (length(complete_set) > 0) { # 添加检查以避免 complete_set 为空时出错
            for (i in complete_set) {
                matrix_indicator <- (4 * i - 3):(4 * i)
                matrix_big_i <- matrix_big[matrix_indicator, ]
                matrix_big_i_inv <- solve(matrix_big_i)
                omega_11 <- matrix_big_i_inv[1:2, 1:2]
                omega_21 <- matrix_big_i_inv[3:4, 1:2]
                omega_12 <- matrix_big_i_inv[1:2, 3:4]
                omega_22 <- matrix_big_i_inv[3:4, 3:4]
                beta_hat_exp_i <- beta_hat_exp[i, ]
                beta_hat_out_i <- beta_hat_out[i, ]
                gamma_fenmu <- (omega_11 + alpha * (omega_21 + omega_12) +
                    alpha^2 * omega_22)
                gamma_fenzi <- ((omega_11 + alpha * omega_21) %*% beta_hat_exp_i +
                    (omega_12 + alpha * omega_22) %*% beta_hat_out_i)

                beta_exp_new[i, ] <- t(solve(gamma_fenmu) %*% gamma_fenzi)
            }
        }
        ## 更新辅数据集的gamma
        if (length(a_reduce_b_set) > 0) { # 添加检查以避免 a_reduce_b_set 为空时出错
            for (i in a_reduce_b_set) {
                matrix_indicator <- (4 * i - 3):(4 * i)
                matrix_big_i <- matrix_big[matrix_indicator, ]
                matrix_o_i <- matrix(0, nrow = 2, ncol = 2)
                matrix_o_i[1, 1] <- matrix_big_i[1, 1]
                matrix_o_i[2, 2] <- matrix_big_i[3, 3]
                matrix_o_i[1, 2] <- matrix_big_i[1, 3]
                matrix_o_i[2, 1] <- matrix_big_i[1, 3]
                matrix_o_i_inv <- solve(matrix_o_i)
                beta_hat_exp_i <- beta_hat_exp[i, 1]
                beta_hat_out_i <- beta_hat_out[i, 1]
                beta_exp_new[i, 2] <- beta_hat_exp[i, 2]

                fenzi <- (matrix_o_i_inv[1, 1] + alpha * matrix_o_i_inv[2, 1]) * beta_hat_exp_i +
                    (matrix_o_i_inv[1, 2] + alpha * matrix_o_i_inv[2, 2]) * beta_hat_out_i
                fenmu <- matrix_o_i_inv[1, 1] + alpha * matrix_o_i_inv[2, 1] +
                    alpha * (matrix_o_i_inv[1, 2] + alpha * matrix_o_i_inv[2, 2])
                beta_exp_new[i, 1] <- fenzi / fenmu
            }
        }
        beta_exp <- beta_exp_new # 更新 beta_exp 以供下一轮迭代或最终返回

        # 更新alpha
        fenzi_complete <- 0
        fenmu_complete <- 0
        fenzi_reduce <- 0
        fenmu_reduce <- 0
        if (length(complete_set) > 0) { # 添加检查
            for (i in complete_set) {
                matrix_indicator <- (4 * i - 3):(4 * i)
                matrix_big_i <- matrix_big[matrix_indicator, ]
                matrix_big_i_inv <- solve(matrix_big_i)
                omega_11 <- matrix_big_i_inv[1:2, 1:2]
                omega_21 <- matrix_big_i_inv[3:4, 1:2]
                omega_12 <- matrix_big_i_inv[1:2, 3:4]
                omega_22 <- matrix_big_i_inv[3:4, 3:4]
                beta_hat_exp_i <- beta_hat_exp[i, ]
                beta_hat_out_i <- beta_hat_out[i, ]
                beta_exp_i_current_iter <- beta_exp[i, ] # 使用当前迭代的 beta_exp
                fenzi_complete <- t(beta_hat_exp_i) %*% (omega_21 %*%
                    (beta_hat_exp_i - beta_exp_i_current_iter) + omega_22 %*% beta_hat_out_i) + fenzi_complete
                fenmu_complete <- t(beta_exp_i_current_iter) %*% omega_22 %*% beta_exp_i_current_iter + fenmu_complete
            }
        }
        if (length(a_reduce_b_set) > 0) { # 添加检查
            for (i in a_reduce_b_set) {
                matrix_indicator <- (4 * i - 3):(4 * i)
                matrix_big_i <- matrix_big[matrix_indicator, ]
                matrix_o_i <- matrix(0, nrow = 2, ncol = 2)
                matrix_o_i[1, 1] <- matrix_big_i[1, 1]
                matrix_o_i[2, 2] <- matrix_big_i[3, 3]
                matrix_o_i[1, 2] <- matrix_big_i[1, 3]
                matrix_o_i[2, 1] <- matrix_big_i[1, 3]
                matrix_o_i_inv <- solve(matrix_o_i)
                beta_hat_exp_i <- beta_hat_exp[i, 1]
                beta_hat_out_i <- beta_hat_out[i, 1]
                beta_exp_i_current_iter <- beta_exp[i, 1] # 使用当前迭代的 beta_exp

                fenzi_reduce <- beta_hat_exp_i * (matrix_o_i_inv[2, 1] *
                    (beta_hat_exp_i - beta_exp_i_current_iter) +
                    matrix_o_i_inv[2, 2] * beta_hat_out_i) + fenzi_reduce
                fenmu_reduce <- beta_exp_i_current_iter^2 * matrix_o_i_inv[2, 2] + fenmu_reduce
            }
        }

        if ((fenmu_complete + fenmu_reduce) == 0) { # 避免除以零
            warning("更新 alpha 时分母为零，迭代可能无法收敛。")
            alpha_new <- alpha # 保持 alpha 不变或采取其他错误处理
        } else {
            alpha_new_val <- (fenzi_complete + fenzi_reduce) / (fenmu_complete + fenmu_reduce)
            alpha_new <- as.numeric(alpha_new_val) # 确保 alpha_new 是标量
        }

        # 检查收敛性
        if (abs(alpha_new - alpha_old) < tol) {
            converged <- TRUE
            alpha <- alpha_new # 更新 alpha 为最终收敛值
            break # 跳出循环
        }
        alpha <- alpha_new # 更新 alpha 以供下一轮迭代
    }

    if (!converged) {
        warning(paste("算法在", max_iter, "次迭代后未达到收敛阈值", tol))
    }

    return(list(alpha_new = alpha, beta_exp_updated = beta_exp, iterations = iterations, converged = converged))
}

# %% 迭代算法

# # 更新r
# ## 构造异质性统计量
# t <- matrix(0, ncol = 1, nrow = n_snps)
# f <- matrix(0, ncol = 1, nrow = n_snps)
# for (i in 1:n_snps) {
#     matrix_indicator <- 4 * i - 1
#     # 计算关于子代的统计量
#     beta_hat_exp_i <- beta_hat_exp[i, 1]
#     beta_exp_i <- beta_exp[i, 1]
#     beta_out_se <- matrix_big[matrix_indicator, 3]
#     t[i] <- (beta_hat_exp_i - alpha * beta_exp_i)^2 / beta_out_se
#     # 计算关于父母的统计量
#     beta_hat_exp_i <- beta_hat_exp[i, 2]
#     beta_exp_i <- beta_exp[i, 2]
#     beta_out_se <- matrix_big[matrix_indicator + 1, 4]
#     f[i] <- (beta_hat_exp_i - alpha * beta_exp_i)^2 / beta_out_se
# }
# a_snps_set <- arrange(data.frame(snps_id = 1:n_snps, statistic = t), t)
# b_snps_set <- arrange(data.frame(snps_id = 1:n_snps, statistic = f), f)

# a_legal_set <- a_snps_set[1:a_legal, 1]
# b_legal_set <- b_snps_set[1:b_legal, 1]
# complete_set <- intersect(a_legal_set, b_legal_set)
# a_reduce_b_set <- setdiff(a_legal_set, b_legal_set)

# # 更新gamma
# beta_exp_new <- beta_exp
# ## 更新全数据集中的gamma
# for (i in complete_set) {
#     matrix_indicator <- (4 * i - 3):(4 * i)
#     matrix_big_i <- matrix_big[matrix_indicator, ]
#     matrix_big_i_inv <- solve(matrix_big_i)
#     omega_11 <- matrix_big_i_inv[1:2, 1:2]
#     omega_21 <- matrix_big_i_inv[3:4, 1:2]
#     omega_12 <- matrix_big_i_inv[1:2, 3:4]
#     omega_22 <- matrix_big_i_inv[3:4, 3:4]
#     beta_hat_exp_i <- beta_hat_exp[i, ]
#     beta_hat_out_i <- beta_hat_out[i, ]
#     gamma_fenmu <- (omega_11 + alpha * (omega_21 + omega_12) +
#         alpha^2 * omega_22)
#     gamma_fenzi <- ((omega_11 + alpha * omega_21) %*% beta_hat_exp_i +
#         (omega_12 + alpha * omega_22) %*% beta_hat_out_i)

#     beta_exp_new[i, ] <- t(solve(gamma_fenmu) %*% gamma_fenzi)
# }
# ## 更新辅数据集的gamma
# for (i in a_reduce_b_set) {
#     matrix_indicator <- (4 * i - 3):(4 * i)
#     matrix_big_i <- matrix_big[matrix_indicator, ]
#     matrix_o_i <- matrix(0, nrow = 2, ncol = 2)
#     matrix_o_i[1, 1] <- matrix_big_i[1, 1]
#     matrix_o_i[2, 2] <- matrix_big_i[3, 3]
#     matrix_o_i[1, 2] <- matrix_big_i[1, 3]
#     matrix_o_i[2, 1] <- matrix_big_i[1, 3]
#     matrix_o_i_inv <- solve(matrix_o_i)
#     beta_hat_exp_i <- beta_hat_exp[i, 1]
#     beta_hat_out_i <- beta_hat_out[i, 1]
#     beta_exp_new[i, 2] <- beta_hat_exp[i, 2]

#     fenzi <- (matrix_o_i_inv[1, 1] + alpha * matrix_o_i_inv[2, 1]) * beta_hat_exp_i +
#         (matrix_o_i_inv[1, 2] + alpha * matrix_o_i_inv[2, 2]) * beta_hat_out_i
#     fenmu <- matrix_o_i_inv[1, 1] + alpha * matrix_o_i_inv[2, 1] +
#         alpha * (matrix_o_i_inv[1, 2] + alpha * matrix_o_i_inv[2, 2])
#     beta_exp_new[i, 1] <- fenzi / fenmu
# }
# beta_exp <- beta_exp_new

# # 更新alpha
# fenzi_complete <- 0
# fenmu_complete <- 0
# fenzi_reduce <- 0
# fenmu_reduce <- 0
# for (i in complete_set) {
#     matrix_indicator <- (4 * i - 3):(4 * i)
#     matrix_big_i <- matrix_big[matrix_indicator, ]
#     matrix_big_i_inv <- solve(matrix_big_i)
#     omega_11 <- matrix_big_i_inv[1:2, 1:2]
#     omega_21 <- matrix_big_i_inv[3:4, 1:2]
#     omega_12 <- matrix_big_i_inv[1:2, 3:4]
#     omega_22 <- matrix_big_i_inv[3:4, 3:4]
#     beta_hat_exp_i <- beta_hat_exp[i, ]
#     beta_hat_out_i <- beta_hat_out[i, ]
#     beta_exp_i <- beta_exp[i, ]
#     fenzi_complete <- t(beta_hat_exp_i) %*% (omega_21 %*%
#         (beta_hat_exp_i - beta_exp_i) + omega_22 %*% beta_hat_out_i) + fenzi_complete
#     fenmu_complete <- t(beta_exp_i) %*% omega_22 %*% beta_exp_i + fenmu_complete
# }
# for (i in a_reduce_b_set) {
#     matrix_indicator <- (4 * i - 3):(4 * i)
#     matrix_big_i <- matrix_big[matrix_indicator, ]
#     matrix_o_i <- matrix(0, nrow = 2, ncol = 2)
#     matrix_o_i[1, 1] <- matrix_big_i[1, 1]
#     matrix_o_i[2, 2] <- matrix_big_i[3, 3]
#     matrix_o_i[1, 2] <- matrix_big_i[1, 3]
#     matrix_o_i[2, 1] <- matrix_big_i[1, 3]
#     matrix_o_i_inv <- solve(matrix_o_i)
#     beta_hat_exp_i <- beta_hat_exp[i, 1]
#     beta_hat_out_i <- beta_hat_out[i, 1]
#     beta_exp_i <- beta_exp[i, 1]

#     fenzi_reduce <- beta_hat_exp_i * (matrix_o_i_inv[2, 1] *
#         (beta_hat_exp_i - beta_exp_i) +
#         matrix_o_i_inv[2, 2] * beta_hat_out_i) + fenzi_reduce
#     fenmu_reduce <- beta_exp_i^2 * matrix_o_i_inv[2, 2] + fenmu_reduce
# }
# alpha_new <- (fenzi_complete + fenzi_reduce) / (fenmu_complete + fenmu_reduce)
# alpha_new

# %% 迭代算法函数cpp版本

# function_test <- perform_iterative_update_rcpp(
#     n_snps, beta_hat_exp, beta_hat_out,
#     alpha_input = 0, matrix_big, a_legal = 3, b_legal = 3,
#     max_iter = 100, tol = 1e-06
# )

# %% 方差计算函数
# n_snps <- 3
# if (exists("generate_mr_trio_data_ultra_updata") &&
#     is.function(generate_mr_trio_data_ultra_updata)) {
#     data_set_example <- generate_mr_trio_data_ultra_updata(
#         beta_ms_oe_exp = 0.5, beta_os_oe_exp = 0.5,
#         beta_fs_oe_out = 0.5, n_snps = n_snps
#     )
# } else {
#     stop("函数 'generate_mr_trio_data_ultra_updata' 未定义。")
# }

# if (exists("perform_fgwas_analysis") && is.function(perform_fgwas_analysis)) {
#     data_set_fgwas_example <- perform_fgwas_analysis(data_set_example, n_snps = n_snps)
# } else {
#     stop("函数 'perform_fgwas_analysis' 未定义。")
# }

# # 定义输入变量
# beta_sigma_exp <- as.matrix(data_set_fgwas_example$beta_sigma_exp)
# beta_sigma_out <- as.matrix(data_set_fgwas_example$beta_sigma_out)
# beta_hat_exp <- as.matrix(data_set_fgwas_example$beta_hat_exp)
# beta_hat_out <- as.matrix(data_set_fgwas_example$beta_hat_out)
# beta_sigma_t_exp <- invert_beta_sigma_out_matrices(beta_sigma_exp)
# beta_sigma_t_out <- invert_beta_sigma_out_matrices(beta_sigma_out)
# matrix_big <- create_combined_diagonal_matrices(
#     beta_sigma_t_exp,
#     beta_sigma_t_out,
#     off_diagonal_elements = NULL
# )
# a <- 0 # 两种集合中不同的snps个数
# b <- 0 # 两种集合中不同的snps个数
# a_legal <- n_snps - a
# b_legal <- n_snps - b
# function_test <- perform_iterative_update_rcpp(
#     n_snps, beta_hat_exp, beta_hat_out,
#     alpha_input = 0, matrix_big, a_legal = 3, b_legal = 3,
#     max_iter = 100, tol = 1e-06
# )
# complete_set <- function_test$complete_set
# a_reduce_b_set <- function_test$a_reduce_b_set
# alpha <- function_test$alpha_new
# beta_exp <- function_test$beta_exp_updated

# %% 定义方差计算函数
calculate_alpha_se <- function(
    complete_set, a_reduce_b_set,
    matrix_big, beta_hat_exp, beta_hat_out, beta_exp, alpha) {
    # 开始进行算法
    h_dimension <- 2 * nrow(beta_hat_exp)
    h_matrix <- matrix(0, nrow = h_dimension + 1, ncol = h_dimension + 1)
    h_a_a <- matrix(0, nrow = 1, ncol = 1)
    h_a_gamma <- matrix(0, nrow = h_dimension, ncol = 1)
    h_gamma_gamma <- matrix(0, nrow = h_dimension, ncol = h_dimension)

    for (i in complete_set) {
        matrix_indicator <- (4 * i - 3):(4 * i)
        matrix_big_i <- matrix_big[matrix_indicator, ]
        matrix_big_i_inv <- solve(matrix_big_i)
        omega_11 <- matrix_big_i_inv[1:2, 1:2]
        omega_21 <- matrix_big_i_inv[3:4, 1:2]
        omega_12 <- matrix_big_i_inv[1:2, 3:4]
        omega_22 <- matrix_big_i_inv[3:4, 3:4]
        beta_hat_exp_i <- beta_hat_exp[i, ]
        beta_hat_out_i <- beta_hat_out[i, ]
        print(i)
        beta_exp_i <- beta_exp[i, ]
        h_a_a <- t(beta_hat_exp_i) %*% omega_22 %*% beta_hat_exp_i + h_a_a
        h_a_gamma[(2 * i - 1):(2 * i), 1] <- (-omega_12) %*% beta_hat_exp_i +
            (omega_12 + omega_21) %*% beta_exp_i - omega_22 %*% beta_hat_out_i +
            2 * alpha * omega_22 %*% beta_hat_exp_i
        h_gamma_gamma[(2 * i - 1):(2 * i), (2 * i - 1):(2 * i)] <- omega_11 + alpha * (omega_12 + omega_21) +
            alpha^2 * omega_22
    }
    for (i in a_reduce_b_set) {
        matrix_indicator <- (4 * i - 3):(4 * i)
        matrix_big_i <- matrix_big[matrix_indicator, ]
        matrix_o_i <- matrix(0, nrow = 2, ncol = 2)
        matrix_o_i[1, 1] <- matrix_big_i[1, 1]
        matrix_o_i[2, 2] <- matrix_big_i[3, 3]
        matrix_o_i[1, 2] <- matrix_big_i[1, 3]
        matrix_o_i[2, 1] <- matrix_big_i[1, 3]
        matrix_o_i_inv <- solve(matrix_o_i)

        beta_hat_exp_i <- beta_hat_exp[i, 1]
        beta_hat_out_i <- beta_hat_out[i, 1]
        beta_exp_i <- beta_exp[i, 1]
        h_a_a <- beta_hat_exp_i^2 * matrix_o_i_inv[2, 2] + h_a_a
        h_a_gamma[(2 * i - 1):(2 * i), 1] <- c(-(matrix_o_i_inv[2, 1] *
            beta_hat_exp_i -
            matrix_o_i_inv[2, 2] * beta_hat_out_i + 2 * beta_exp_i * (matrix_o_i_inv[2, 1] +
                alpha * matrix_o_i_inv[2, 2])), 0)
        h_gamma_gamma_input <- matrix_o_i_inv[1, 1] +
            alpha * (matrix_o_i_inv[1, 2] + matrix_o_i_inv[2, 1]) +
            alpha^2 * matrix_o_i_inv[2, 2]
        h_gamma_gamma_input_matrix <- matrix(c(h_gamma_gamma_input, 0, 0, 0),
            nrow = 2, ncol = 2
        )
        h_gamma_gamma[(2 * i - 1):(2 * i), (2 * i - 1):(2 * i)] <-
            h_gamma_gamma_input_matrix
    }
    h_matrix[1, 1] <- h_a_a
    h_matrix[-1, 1] <- h_a_gamma
    h_matrix[1, -1] <- t(h_a_gamma)
    h_matrix[-1, -1] <- h_gamma_gamma

    # 检查 h_matrix 中是否存在全为0的行和列
    alpha_se_value <- NA # 初始化 alpha_se_value
    if (nrow(h_matrix) > 0 && ncol(h_matrix) > 0) {
        # 找出全为0的行
        zero_rows <- which(apply(h_matrix, 1, function(row) all(row == 0)))
        # 找出全为0的列
        zero_cols <- which(apply(h_matrix, 2, function(col) all(col == 0)))

        # 找出同时为0的行和列的索引 (确保行列索引一致才移除)
        # 注意：R的索引是从1开始的
        common_zero_indices <- intersect(zero_rows, zero_cols)

        if (length(common_zero_indices) > 0) {
            cat("发现全零行和列的索引:", common_zero_indices, "\n")
            # 创建新的索引以排除这些行和列
            rows_to_keep <- setdiff(1:nrow(h_matrix), common_zero_indices)
            cols_to_keep <- setdiff(1:ncol(h_matrix), common_zero_indices)

            if (length(rows_to_keep) > 0 && length(cols_to_keep) > 0) {
                h_matrix_cleaned <- h_matrix[rows_to_keep,
                    cols_to_keep,
                    drop = FALSE
                ]
                cat("清理后的 h_matrix 维度:", dim(h_matrix_cleaned), "\n")
                # 更新 alpha_se 的计算，如果 h_matrix_cleaned 非空且可逆
                if (nrow(h_matrix_cleaned) > 0 && ncol(h_matrix_cleaned) > 0 && det(h_matrix_cleaned) != 0) {
                    alpha_se_value <- solve(h_matrix_cleaned)[1, 1]
                    cat("使用清理后的 h_matrix 重新计算 alpha_se:", alpha_se_value, "\n")
                } else {
                    warning("清理后的 h_matrix 为空或奇异，无法计算 alpha_se。")
                    # alpha_se_value 保持 NA
                }
            } else {
                warning("清理后没有剩余的行或列，h_matrix 变为空。")
                # alpha_se_value 保持 NA
            }
        } else {
            cat("未发现同时为全零的行和列。\n")
            # 如果没有需要移除的行和列，并且原始h_matrix可逆
            if (nrow(h_matrix) > 0 && ncol(h_matrix) > 0 && det(h_matrix) != 0) {
                alpha_se_value <- solve(h_matrix)[1, 1]
            } else {
                warning("原始 h_matrix 为空或奇异，无法计算 alpha_se。")
                # alpha_se_value 保持 NA
            }
        }
    } else {
        warning("h_matrix 为空，无法进行处理。")
        # alpha_se_value 保持 NA
    }
    return(sqrt(alpha_se_value))
}

# %% 定义计算 bic的函数
calculate_bic <- function(function_test_output,
                          beta_hat_exp,
                          beta_hat_out,
                          matrix_big,
                          min_n,
                          a = 0,
                          b = 0) {
    # --- 1. 从输入列表中解析所需变量 ---
    complete_set <- function_test_output$complete_set
    a_reduce_b_set <- function_test_output$a_reduce_b_set
    alpha <- function_test_output$alpha_new
    beta_exp <- function_test_output$beta_exp_updated
    n_snps <- nrow(beta_hat_exp)

    # --- 2. 根据SNP分组，分别计算各部分的似然贡献 ---

    # 2.1 计算 complete_set 部分
    l_complete <- 0
    if (length(complete_set) > 0) {
        for (i in complete_set) {
            matrix_index_all <- (4 * i - 3):(4 * i)
            matrix_big_i_inv <- solve(matrix_big[matrix_index_all, ])
            beta_hat_exp_i <- beta_hat_exp[i, ]
            beta_hat_out_i <- beta_hat_out[i, ]
            beta_exp_i <- beta_exp[i, ]
            d_diff <- c(
                beta_hat_exp_i - beta_exp_i,
                beta_hat_out_i - alpha * beta_exp_i
            )
            l_complete <- l_complete + t(d_diff) %*% matrix_big_i_inv %*% d_diff
        }
    }

    # 2.2 计算 a_reduce_b_set 部分
    l_reduce <- 0
    if (length(a_reduce_b_set) > 0) {
        for (i in a_reduce_b_set) {
            matrix_index_all <- (4 * i - 3):(4 * i)
            # 注意：保留了您原始代码的逻辑。
            # 此处使用 matrix_index_all 提取子矩阵。
            matrix_big_i_sub <- matrix_big[matrix_index_all, ]

            matrix_o_i <- matrix(0, nrow = 2, ncol = 2)
            matrix_o_i[1, 1] <- matrix_big_i_sub[1, 1]
            matrix_o_i[2, 2] <- matrix_big_i_sub[3, 3]
            matrix_o_i[1, 2] <- matrix_big_i_sub[1, 3]
            matrix_o_i[2, 1] <- matrix_big_i_sub[1, 3]
            matrix_o_i_inv <- solve(matrix_o_i)

            beta_hat_exp_i <- beta_hat_exp[i, ]
            beta_hat_out_i <- beta_hat_out[i, ]
            beta_exp_i <- beta_exp[i, ]
            d_o <- c(beta_hat_exp_i[1] - beta_exp_i[1], beta_hat_out_i[1] - alpha * beta_exp_i[1])

            l_reduce <- l_reduce + t(d_o) %*% matrix_o_i_inv %*% d_o +
                (beta_hat_exp_i[2] - beta_exp_i[2])^2 / matrix_big_i_sub[2, 2]
        }
    }

    # 2.3 计算剩余集合 (a_c_b_c_set) 部分
    a_c_b_c_set <- setdiff(1:n_snps, c(complete_set, a_reduce_b_set))
    l_cc <- 0
    if (length(a_c_b_c_set) > 0) {
        for (i in a_c_b_c_set) {
            matrix_index_exp <- (4 * i - 3):(4 * i - 2)
            # 注意：保留您的原始逻辑，但已修正以确保可运行。
            # 原始代码在维度上可能存在问题，这里提取 2x2 子矩阵并求逆。
            matrix_exp_i_inv <- solve(matrix_big[
                matrix_index_exp, 1:2
            ])

            beta_hat_exp_i <- beta_hat_exp[i, ]
            beta_exp_i <- beta_exp[i, ]

            l_cc <- l_cc + t(beta_hat_exp_i - beta_exp_i) %*%
                matrix_exp_i_inv %*% (beta_hat_exp_i - beta_exp_i)
        }
    }

    # --- 3. 计算最终的BIC值 ---
    bic <- l_complete + l_reduce + l_cc + log(min_n) * (a + b)

    return(bic)
}


# %% bic 计算代码版



# %% 完整的cml_overlap_2
robust_cml_fit <- function(
    n_snps, beta_hat_exp, beta_hat_out,
    matrix_big, a_legal, b_legal,
    max_iter = 100, tol = 1e-6, n_starts = 5, alpha_start_range = c(-0.5, 0.5),
    digits = 4) {
    results_list <- list()

    # 多次运行迭代算法，每次使用不同的随机起点
    for (i in 1:n_starts) {
        # 从指定范围生成一个随机的 alpha 起始值
        random_alpha_start <- runif(1, min = alpha_start_range[1], max = alpha_start_range[2])

        # 使用 try 来捕获迭代过程中可能发生的错误
        result <- try(
            perform_iterative_update_rcpp(
                n_snps = n_snps,
                beta_hat_exp = beta_hat_exp,
                beta_hat_out = beta_hat_out,
                alpha_input = random_alpha_start,
                matrix_big = matrix_big,
                a_legal = a_legal,
                b_legal = b_legal,
                max_iter = max_iter,
                tol = tol
            ),
            silent = TRUE
        )

        # 如果迭代成功，则保存结果
        if (!inherits(result, "try-error") && result$converged) {
            results_list[[length(results_list) + 1]] <- result
        }
    }

    # 检查是否有任何成功的收敛结果
    if (length(results_list) == 0) {
        stop("在所有尝试中，算法均未能收敛。")
    }

    # 提取所有收敛的 alpha 值，并四舍五入到指定小数位数
    final_alphas <- sapply(results_list, function(res) round(res$alpha_new, digits = digits))

    # 找到出现次数最多的 alpha 值（众数）
    unique_alphas <- unique(final_alphas)
    mode_alpha <- unique_alphas[which.max(tabulate(match(final_alphas, unique_alphas)))]

    # 找到与众数 alpha 对应的第一个结果
    best_result_index <- which(final_alphas == mode_alpha)[1]
    best_result <- results_list[[best_result_index]]

    # 返回最稳定的结果
    return(best_result)
}

# %% cml本体
cml_family_overlap_2_init <- function(
    n_snps, beta_hat_exp, beta_hat_out,
    matrix_big,
    max_iter = 100, tol = 1e-6, n_starts = 5,
    alpha_start_range = c(-0.01, 0.01),
    digits = 4) {
    # 生成组合列表
    if (n_snps > 1) {
        # 使用 lapply 和 do.call(rbind, ...) 的高效方法来创建数据框
        legal_pairs_list <- lapply(2:n_snps, function(a) {
            data.frame(a_legal = a, b_legal = 1:(a))
        })
        legal_pairs_table <- do.call(rbind, legal_pairs_list)

        # 打印生成的表格头部以供预览
        cat("\n--- a_legal 和 b_legal 组合表格 ---\n")
        print(head(legal_pairs_table))
        cat("...\n")
        cat("总共生成了", nrow(legal_pairs_table), "个组合。\n")
    } else {
        # 处理 n_snps 不足以生成任何有效组合的情况
        cat("\n'n_snps' 的值 (", n_snps, ") 小于或等于 1，无法生成任何 b_legal < a_legal 的组合。\n")
        legal_pairs_table <- data.frame(a_legal = integer(0), b_legal = integer(0))
    }
    # 默认初始值
    legal_pairs_table$bic <- 0
    legal_pairs_table$alpha <- 0
    legal_pairs_table$alpha_se <- 0

    for (i in 1:nrow(legal_pairs_table)) {
        a_legal_i <- legal_pairs_table[i, 1]
        b_legal_i <- legal_pairs_table[i, 2]
        cml_fit_i <- robust_cml_fit(
            n_snps, beta_hat_exp, beta_hat_out,
            matrix_big, a_legal_i, b_legal_i,
            max_iter = max_iter, tol = tol, n_starts = n_starts,
            alpha_start_range = alpha_start_range,
            digits = digits
        )
        legal_pairs_table$alpha[i] <- cml_fit_i$alpha_new
        # 准备方差计算的输入
        complete_set_i <- cml_fit_i$complete_set
        a_reduce_b_set <- cml_fit_i$a_reduce_b_set
        beta_exp_i <- cml_fit_i$beta_exp_updated
        alpha_i <- cml_fit_i$alpha_new

        legal_pairs_table$alpha_se[i] <- calculate_alpha_se_rcpp(
            complete_set_i,
            a_reduce_b_set,
            matrix_big,
            beta_hat_exp,
            beta_hat_out,
            beta_exp_i,
            alpha_i
        )
        legal_pairs_table$bic[i] <- calculate_bic(cml_fit_i,
            beta_hat_exp,
            beta_hat_out,
            matrix_big,
            min_n = 1000,
            a = n_snps - a_legal_i,
            b = n_snps - b_legal_i
        )
    }
    #
    weights <- exp(-(legal_pairs_table$bic / 2))

    # 归一化权重
    # 在 R 中，is.na() 函数可以同时检测 NA 和 NaN，这通常是更稳健的做法。
    valid_rows_filter <- !is.na(legal_pairs_table$alpha_se)

    # 2. 根据过滤器筛选出有效的数据表和对应的原始权重
    # 后续的所有计算都将基于这个有效子集 (valid_table)
    valid_table <- legal_pairs_table[valid_rows_filter, ]
    valid_weights <- weights[valid_rows_filter]

    # 3. 检查是否存在有效行，以防止因数据为空而出错
    if (nrow(valid_table) > 0 && sum(valid_weights, na.rm = TRUE) > 0) {
        # 4. 对筛选后的有效权重进行归一化
        # 注意：分母现在是有效权重的总和
        normalized_weights <- valid_weights / sum(valid_weights)
        valid_table$weight <- normalized_weights

        # 5. 计算加权 alpha (仅使用有效数据)
        weighted_alpha <- sum(valid_table$weight * valid_table$alpha)

        # 6. 计算加权标准误 (沿用您原有的加权方式)
        # 这实际上是标准误（alpha_se）的加权平均值
        weighted_se <- sum(valid_table$weight * sqrt((valid_table$alpha_se)^2 +
            (valid_table$alpha - weighted_alpha)^2))
        result <- list(
            weighted_alpha = weighted_alpha,
            weighted_se = weighted_se,
            legal_pairs_table = legal_pairs_table
        )
        return(result)
    }
}
# %% cml鲁棒性增强
cml_family_overlap_2_robust <- function(
    n_snps, beta_hat_exp, beta_hat_out,
    matrix_big,
    max_iter = 100, tol = 1e-6, n_starts = 5,
    alpha_start_range = c(-0.01, 0.01),
    digits = 4) {
    # --- 1. 组合生成部分 (您的代码已经很高效) ---
    if (n_snps > 1) {
        legal_pairs_list <- lapply(2:n_snps, function(a) {
            # 修正：b_legal 应该小于 a_legal，所以是 1:(a-1)
            # 如果 b_legal 可以等于 a_legal，那么您的 1:a 是正确的
            # 这里我假设 b 必须严格小于 a，如果不是，请改回 1:a
            if (a > 1) {
                return(data.frame(a_legal = a, b_legal = 1:(a - 1)))
            }
        })
        legal_pairs_table <- do.call(rbind, legal_pairs_list)

        cat("\n--- a_legal 和 b_legal 组合表格 ---\n")
        print(head(legal_pairs_table))
        cat("...\n")
        cat("总共生成了", nrow(legal_pairs_table), "个待计算的组合。\n\n")
    } else {
        cat("\n'n_snps' 的值 (", n_snps, ") 小于等于 1，无法生成任何 b_legal < a_legal 的组合。\n")
        # 直接返回一个空结果，而不是继续执行
        return(list(
            weighted_alpha = NA,
            weighted_se = NA,
            legal_pairs_table = data.frame(a_legal = integer(0), b_legal = integer(0))
        ))
    }

    # --- 2. 结果预分配 ---
    # 预先分配列表来存储结果，这比在循环中反复扩展数据框更高效
    results_list <- vector("list", nrow(legal_pairs_table))

    # --- 3. 核心计算循环 (增加 tryCatch 和日志) ---
    for (i in 1:nrow(legal_pairs_table)) {
        a_legal_i <- legal_pairs_table$a_legal[i]
        b_legal_i <- legal_pairs_table$b_legal[i]

        # 增加进度条反馈
        cat(sprintf("正在处理第 %d / %d 个组合: a=%d, b=%d\n", i, nrow(legal_pairs_table), a_legal_i, b_legal_i))

        # 【核心改动】使用 tryCatch 包裹所有可能出错的计算
        results_list[[i]] <- tryCatch(
            {
                # --- 尝试执行的代码块 ---
                cml_fit_i <- robust_cml_fit(
                    n_snps, beta_hat_exp, beta_hat_out,
                    matrix_big, a_legal_i, b_legal_i,
                    max_iter = max_iter, tol = tol, n_starts = n_starts,
                    alpha_start_range = alpha_start_range,
                    digits = digits
                )

                complete_set_i <- cml_fit_i$complete_set
                a_reduce_b_set <- cml_fit_i$a_reduce_b_set
                beta_exp_i <- cml_fit_i$beta_exp_updated
                alpha_i <- cml_fit_i$alpha_new

                alpha_se_i <- calculate_alpha_se_rcpp(
                    complete_set_i, a_reduce_b_set, matrix_big,
                    beta_hat_exp, beta_hat_out, beta_exp_i, alpha_i
                )

                bic_i <- calculate_bic(
                    cml_fit_i, beta_hat_exp, beta_hat_out, matrix_big,
                    min_n = 1000, a = n_snps - a_legal_i, b = n_snps - b_legal_i
                )

                # 如果成功，返回一个包含所有结果的列表
                list(alpha = alpha_i, alpha_se = alpha_se_i, bic = bic_i)
            },
            error = function(e) {
                # --- 如果发生任何错误，则执行此处的代码 ---
                cat(sprintf("  -> 组合 (a=%d, b=%d) 计算失败。错误信息: %s\n", a_legal_i, b_legal_i, e$message))
                # 返回包含 NA 的列表，以标记失败
                return(list(alpha = NA, alpha_se = NA, bic = NA))
            }
        )
    }

    # --- 4. 整合循环结果 ---
    # 将列表中的结果高效地合并回数据框
    results_df <- do.call(rbind, lapply(results_list, as.data.frame))
    legal_pairs_table <- cbind(legal_pairs_table, results_df)

    # --- 5. 模型平均 (增加 is.finite 检查) ---
    cat("\n--- 所有组合计算完成，开始进行模型平均 ---\n")

    # 【核心改动】使用 is.finite() 来过滤掉 NA, NaN, Inf, -Inf
    valid_rows_filter <- is.finite(legal_pairs_table$alpha) &
        is.finite(legal_pairs_table$alpha_se) &
        is.finite(legal_pairs_table$bic)

    valid_table <- legal_pairs_table[valid_rows_filter, ]
    cat("在", nrow(legal_pairs_table), "个组合中，有", nrow(valid_table), "个得到了有效的数值结果。\n")

    if (nrow(valid_table) > 0) {
        # 计算权重时也要注意BIC可能为Inf或-Inf的情况
        # 通过在有效数据子集上计算来避免此问题
        weights <- exp(-(valid_table$bic / 2))

        # 再次检查权重是否有效
        if (any(!is.finite(weights)) || sum(weights, na.rm = TRUE) == 0) {
            warning("计算出的权重包含无效值(Inf/NaN)或总和为零，无法进行加权平均。")
            return(list(weighted_alpha = NA, weighted_se = NA, legal_pairs_table = legal_pairs_table))
        }

        normalized_weights <- weights / sum(weights)
        valid_table$weight <- normalized_weights

        weighted_alpha <- sum(valid_table$weight * valid_table$alpha)

        # 计算加权SE的公式保持不变
        weighted_se <- sum(valid_table$weight * sqrt((valid_table$alpha_se)^2 +
            (valid_table$alpha - weighted_alpha)^2))

        result <- list(
            weighted_alpha = weighted_alpha,
            weighted_se = weighted_se,
            legal_pairs_table = legal_pairs_table # 返回包含所有尝试（包括失败）的完整表格
        )
        return(result)
    } else {
        warning("没有任何组合能够成功计算出有效结果。")
        return(list(
            weighted_alpha = NA,
            weighted_se = NA,
            legal_pairs_table = legal_pairs_table
        ))
    }
}
