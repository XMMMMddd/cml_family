# %% 加载包
library(Matrix)
library(dplyr)
library(Rcpp) # 需要 Rcpp 包来编译和加载 C++ 代码
# %% 加载需要用到的包用于开发
source("FGWAS 优化版本/FGWAS 函数.R")
source("样本生成函数/三联体家庭结构生成函数.R")
# %% 加载 C++ 代码
cpp_file_path <- file.path("cml家庭函数", "cml家庭函数overlap版本", "perform_iterative_update_rcpp.cpp")

Rcpp::sourceCpp(cpp_file_path)
print(paste("成功加载 Rcpp 函数从:", cpp_file_path))
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
# %% 拟合算法(自适应选择有效的snps位点)
n_snps <- 3
if (exists("generate_mr_trio_data_ultra_updata") &&
    is.function(generate_mr_trio_data_ultra_updata)) {
    data_set_example <- generate_mr_trio_data_ultra_updata(
        beta_ms_oe_exp = 0.5, beta_os_oe_exp = 0.5,
        beta_fs_oe_out = 0.5, n_snps = n_snps
    )
} else {
    stop("函数 'generate_mr_trio_data_ultra_updata' 未定义。")
}

if (exists("perform_fgwas_analysis") && is.function(perform_fgwas_analysis)) {
    data_set_fgwas_example <- perform_fgwas_analysis(data_set_example, n_snps = n_snps)
} else {
    stop("函数 'perform_fgwas_analysis' 未定义。")
}

# 定义输入变量
beta_sigma_exp <- as.matrix(data_set_fgwas_example$beta_sigma_exp)
beta_sigma_out <- as.matrix(data_set_fgwas_example$beta_sigma_out)
beta_hat_exp <- as.matrix(data_set_fgwas_example$beta_hat_exp)
beta_hat_out <- as.matrix(data_set_fgwas_example$beta_hat_out)
beta_sigma_t_exp <- invert_beta_sigma_out_matrices(beta_sigma_exp)
beta_sigma_t_out <- invert_beta_sigma_out_matrices(beta_sigma_out)
matrix_big <- create_combined_diagonal_matrices(
    beta_sigma_t_exp,
    beta_sigma_t_out,
    off_diagonal_elements = NULL
)
a <- 0 # 两种集合中不同的snps个数
b <- 0 # 两种集合中不同的snps个数
a_legal <- n_snps - a
b_legal <- n_snps - b
# 定义初始值
alpha <- 0
beta_exp <- beta_hat_exp
r <- matrix(0, ncol = 2, nrow = n_snps)

# %% 迭代算法

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
beta_exp_new <- beta_exp
## 更新全数据集中的gamma
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
## 更新辅数据集的gamma
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
beta_exp <- beta_exp_new

# 更新alpha
fenzi_complete <- 0
fenmu_complete <- 0
fenzi_reduce <- 0
fenmu_reduce <- 0
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
    beta_exp_i <- beta_exp[i, ]
    fenzi_complete <- t(beta_hat_exp_i) %*% (omega_21 %*%
        (beta_hat_exp_i - beta_exp_i) + omega_22 %*% beta_hat_out_i) + fenzi_complete
    fenmu_complete <- t(beta_exp_i) %*% omega_22 %*% beta_exp_i + fenmu_complete
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

    fenzi_reduce <- beta_hat_exp_i * (matrix_o_i_inv[2, 1] *
        (beta_hat_exp_i - beta_exp_i) +
        matrix_o_i_inv[2, 2] * beta_hat_out_i) + fenzi_reduce
    fenmu_reduce <- beta_exp_i^2 * matrix_o_i_inv[2, 2] + fenmu_reduce
}
alpha_new <- (fenzi_complete + fenzi_reduce) / (fenmu_complete + fenmu_reduce)
alpha_new

# %% 迭代算法函数

function_test <- perform_iterative_update_rcpp(
    n_snps, beta_hat_exp, beta_hat_out,
    alpha_input = 0, matrix_big, a_legal = 3, b_legal = 3,
    max_iter = 100, tol = 1e-06
)

# %% 方差计算函数
n_snps <- 3
if (exists("generate_mr_trio_data_ultra_updata") &&
    is.function(generate_mr_trio_data_ultra_updata)) {
    data_set_example <- generate_mr_trio_data_ultra_updata(
        beta_ms_oe_exp = 0.5, beta_os_oe_exp = 0.5,
        beta_fs_oe_out = 0.5, n_snps = n_snps
    )
} else {
    stop("函数 'generate_mr_trio_data_ultra_updata' 未定义。")
}

if (exists("perform_fgwas_analysis") && is.function(perform_fgwas_analysis)) {
    data_set_fgwas_example <- perform_fgwas_analysis(data_set_example, n_snps = n_snps)
} else {
    stop("函数 'perform_fgwas_analysis' 未定义。")
}

# 定义输入变量
beta_sigma_exp <- as.matrix(data_set_fgwas_example$beta_sigma_exp)
beta_sigma_out <- as.matrix(data_set_fgwas_example$beta_sigma_out)
beta_hat_exp <- as.matrix(data_set_fgwas_example$beta_hat_exp)
beta_hat_out <- as.matrix(data_set_fgwas_example$beta_hat_out)
beta_sigma_t_exp <- invert_beta_sigma_out_matrices(beta_sigma_exp)
beta_sigma_t_out <- invert_beta_sigma_out_matrices(beta_sigma_out)
matrix_big <- create_combined_diagonal_matrices(
    beta_sigma_t_exp,
    beta_sigma_t_out,
    off_diagonal_elements = NULL
)
a <- 0 # 两种集合中不同的snps个数
b <- 0 # 两种集合中不同的snps个数
a_legal <- n_snps - a
b_legal <- n_snps - b
function_test <- perform_iterative_update_rcpp(
    n_snps, beta_hat_exp, beta_hat_out,
    alpha_input = 0, matrix_big, a_legal = 3, b_legal = 2,
    max_iter = 100, tol = 1e-06
)
complete_set <- function_test$complete_set
a_reduce_b_set <- function_test$a_reduce_b_set
alpha <- function_test$alpha_new
beta_exp <- function_test$beta_exp_updated


h_dimension <- length(complete_set) * 2 + length(a_reduce_b_set) * 2
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
    beta_exp_i <- beta_exp[i, ]
    h_a_a <- t(beta_hat_exp_i) %*% omega_22 %*% beta_hat_exp_i + h_a_a
    h_a_gamma[(2 * i - 1):(2 * i), 1] <- -(omega_12) %*% beta_hat_exp_i +
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
    h_a_gamma[(2 * i - 1):(2 * i), 1] <- c(-(matrix_o_i_inv[2, 1] * beta_hat_exp_i -
        matrix_o_i_inv[2, 2] * beta_hat_out_i + 2 * beta_exp_i * (matrix_o_i_inv[2, 1] +
            alpha * matrix_o_i_inv[2, 2])), 0)
    h_gamma_gamma_input <- matrix_o_i_inv[1, 1] +
        alpha * (matrix_o_i_inv[1, 2] + matrix_o_i_inv[2, 1]) + alpha^2 * matrix_o_i_inv[2, 2]
    h_gamma_gamma_input_matrix <- matrix(c(h_gamma_gamma_input, 0, 0, 0), nrow = 2, ncol = 2)
    h_gamma_gamma[(2 * i - 1):(2 * i), (2 * i - 1):(2 * i)] <- h_gamma_gamma_input_matrix
}
h_matrix[1, 1] <- h_a_a
h_matrix[-1, 1] <- h_a_gamma
h_matrix[1, -1] <- t(h_a_gamma)
h_matrix[-1, -1] <- h_gamma_gamma

# 检查 h_matrix 中是否存在全为0的行和列
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
             h_matrix_cleaned <- h_matrix[rows_to_keep, cols_to_keep, drop = FALSE]
             cat("清理后的 h_matrix 维度:", dim(h_matrix_cleaned), "\n")
             # 更新 alpha_se 的计算，如果 h_matrix_cleaned 非空且可逆
             if (nrow(h_matrix_cleaned) > 0 && ncol(h_matrix_cleaned) > 0 && det(h_matrix_cleaned) != 0) {
                 alpha_se <- solve(h_matrix_cleaned)[1, 1]
                 cat("使用清理后的 h_matrix 重新计算 alpha_se:", alpha_se, "\n")
             } else {
                 warning("清理后的 h_matrix 为空或奇异，无法计算 alpha_se。")
                 alpha_se <- NA # 或者其他合适的默认/错误值
             }
        } else {
            warning("清理后没有剩余的行或列，h_matrix 变为空。")
            alpha_se <- NA # 或者其他合适的默认/错误值
        }
    } else {
        cat("未发现同时为全零的行和列。\n")
        # 如果没有需要移除的行和列，并且原始h_matrix可逆
        if (nrow(h_matrix) > 0 && ncol(h_matrix) > 0 && det(h_matrix) != 0) {
            alpha_se <- solve(h_matrix)[1, 1]
        } else {
            warning("原始 h_matrix 为空或奇异，无法计算 alpha_se。")
            alpha_se <- NA
        }
    }
} else {
    warning("h_matrix 为空，无法进行处理。")
    alpha_se <- NA # 或者其他合适的默认/错误值
}
alpha_se
