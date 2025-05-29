# %%
library(RcppArmadillo)
library(Rcpp)
sourceCpp("cml家庭函数/cml家庭函数二代/cml_family_rcpp.cpp")


# %%
cml_family_ver_2 <- function(
    beta_hat_exp, beta_hat_out,
    beta_sigma_exp, beta_sigma_out,
    permissible_error = 0.0000001,
    maximum_number_of_iterations = 1000,
    illegal_snps = 0, initial_alpha = 1) {
    # --- 1. 初始化参数 ---

    snp_number <- nrow(beta_hat_exp)

    # 初始化因果效应 alpha (α)
    alpha <- 0

    r_snp <- matrix(rep(0, 2 * snp_number), ncol = 2)

    # 暴露 beta 的初始估计值，使用传入的观测值。该值将在迭代中更新。
    beta_exp <- beta_hat_exp

    # --- 2. 设置迭代终止条件 ---
    # 存储上一次迭代的 alpha 值，用于检查收敛性。
    # 初始化为一个与当前 alpha 相差较大的值，以确保循环至少执行一次。
    alpha_last <- alpha + 1000

    # 循环终止条件标志。
    termination_conditions <- TRUE
    # 迭代次数计数器。
    number_of_iterations <- 0

    # --- 3. 主迭代循环 ---
    # 当 (上一次 alpha 与当前 alpha 的绝对差 > 容许误差) 且 (迭代次数 < 最大迭代次数) 时，继续循环。
    while (termination_conditions) {
        number_of_iterations <- number_of_iterations + 1
        # %%
        termination_conditions <- (abs(alpha_last - alpha) >
            permissible_error &&
            number_of_iterations < maximum_number_of_iterations)
        alpha_last <- alpha
        # 开始更新参数
        # 更新r
        # 声明snps差的向量
        diff_snp <- rep(0, snp_number)
        for (i in 1:snp_number) {
            current_index <- (2 * i - 1):(2 * i)
            current_beta_exp <- beta_exp[i, ]
            current_beta_hat_out <- beta_hat_out[i, ]
            current_beta_sigma_out <- beta_sigma_out[current_index, ]
            diff_beta <- current_beta_hat_out - as.numeric(alpha) * current_beta_exp
            diff_snp[i] <- t(diff_beta) %*% current_beta_sigma_out %*% diff_beta
            r_snp[i, ] <- diff_beta
        }
        diff_r_dataframe <- data.frame(diff_snp,
            r_snp_exp = r_snp[, 1],
            r_snp_out = r_snp[, 2]
        )
        diff_r_dataframe_sort <- arrange(diff_r_dataframe, desc(diff_snp))
        r_snp_new <- matrix(rep(0, 2 * snp_number), ncol = 2)
        if (illegal_snps == 0) {
            r_snp_new <- matrix(rep(0, 2 * snp_number), ncol = 2)
        } else {
            r_snp_new[1:illegal_snps, 1] <- diff_r_dataframe_sort[1:illegal_snps, 2]
            r_snp_new[1:illegal_snps, 2] <- diff_r_dataframe_sort[1:illegal_snps, 3]
        }
        # 更新r_snp
        r_snp <- r_snp_new

        # 更新 beta_exp
        beta_exp_new <- beta_exp
        for (i in 1:snp_number) {
            current_index <- (2 * i - 1):(2 * i)
            current_beta_exp <- beta_exp[i, ]
            current_beta_hat_exp <- beta_hat_exp[i, ]
            current_beta_hat_out <- beta_hat_out[i, ]

            current_beta_sigma_exp <- beta_sigma_exp[current_index, ]
            current_beta_sigma_out <- beta_sigma_out[current_index, ]

            denominator_portion <- solve((current_beta_sigma_exp +
                as.numeric(alpha^2) * current_beta_sigma_out))

            molecular_portion <- current_beta_sigma_exp %*% current_beta_hat_exp +
                as.numeric(alpha) * (current_beta_sigma_out %*% (current_beta_hat_out - r_snp[i, ]))

            beta_exp_new[i, ] <- denominator_portion %*% molecular_portion
        }
        beta_exp <- beta_exp_new

        # 更新alpha
        denominator_portion <- 0
        molecular_portion <- 0
        for (i in 1:snp_number) {
            current_index <- (2 * i - 1):(2 * i)
            current_beta_exp <- beta_exp[i, ]
            current_beta_hat_exp <- beta_hat_exp[i, ]
            current_beta_hat_out <- beta_hat_out[i, ]
            current_beta_sigma_exp <- beta_sigma_exp[current_index, ]
            current_beta_sigma_out <- beta_sigma_out[current_index, ]
            current_r_snp <- r_snp[i, ]
            molecular_portion <- t(current_beta_hat_out - current_r_snp) %*%
                current_beta_sigma_out %*% current_beta_exp +
                molecular_portion
            denominator_portion <- t(current_beta_exp) %*%
                current_beta_sigma_out %*% current_beta_exp +
                denominator_portion
        }
        alpha_new <- molecular_portion / denominator_portion
        alpha <- alpha_new
    }

    # --- 4. 检查收敛性并返回结果 ---
    # 判断算法是否因为达到最大迭代次数而停止。
    if (number_of_iterations >= maximum_number_of_iterations &&
        abs(alpha_last - alpha) > permissible_error) { # 也检查误差条件
        convergence <- FALSE
    } else {
        convergence <- TRUE
    }

    # 将结果组织成一个列表。
    # 返回的 alpha, beta_exp, r_snp 都是矩阵或向量形式，取决于其内部计算。
    # alpha 通常是一个标量。
    results <- list(
        alpha = alpha,
        beta_exp = beta_exp,
        r_snp = r_snp,
        convergence = convergence,
        iterations = number_of_iterations # 也可返回迭代次数
    )

    return(results)
}

# %%
cml_multi_start_optimization_ver_2 <- function(
    beta_hat_exp, beta_hat_out,
    beta_sigma_exp, beta_sigma_out,
    permissible_error = 1e-7,
    maximum_number_of_iterations = 1000,
    illegal_snps = 0, # R 函数的参数
    number_of_attempts = 5,
    consensus_proportion = 0.5) {
    if (number_of_attempts <= 0) {
        warning("number_of_attempts 必须是正数。")
        return(list(
            selected_model_result = NA,
            num_in_consensus_cluster = 0,
            consensus_cluster_alphas = NA,
            consensus_cluster_initial_alphas = NA,
            all_initial_alphas = numeric(0),
            all_resulting_alphas = NA,
            all_convergence_status = NA
        ))
    }

    initial_alphas_used <- rnorm(number_of_attempts, mean = 0, sd = 0.05)
    model_results <- vector("list", length = number_of_attempts)
    original_indices <- 1:number_of_attempts # 用于追踪原始索引

    for (i in 1:number_of_attempts) {
        # 假设 cml_family_ver_2_rcpp 函数已加载并可用
        # 并且其参数签名如之前讨论：
        # beta_hat_exp, beta_hat_out, beta_sigma_exp_stack, beta_sigma_out_stack,
        # initial_alpha, permissible_error, maximum_number_of_iterations, illegal_snps_arg
        model_results[[i]] <- cml_family_ver_2_rcpp(
            beta_hat_exp = beta_hat_exp,
            beta_hat_out = beta_hat_out,
            beta_sigma_exp = beta_sigma_exp, # 传递给Rcpp的参数
            beta_sigma_out = beta_sigma_out, # 传递给Rcpp的参数
            initial_alpha = initial_alphas_used[i],
            permissible_error = permissible_error,
            maximum_number_of_iterations = maximum_number_of_iterations,
            illegal_snps = illegal_snps # 将R函数的illegal_snps传递给Rcpp
        )
    }

    # 提取所有结果的 alpha 值和收敛状态
    # 使用安全的提取方式，以防某个model_results[[i]]不完整（尽管Rcpp应保证结构）
    all_resulting_alphas <- sapply(model_results, function(res) {
        if (is.list(res) && "alpha" %in% names(res)) res$alpha else NA
    })
    all_convergence_status <- sapply(model_results, function(res) {
        if (is.list(res) && "convergence" %in% names(res)) res$convergence else FALSE
    })

    # 筛选出成功收敛的结果及其原始索引和alpha值
    converged_mask <- all_convergence_status & !is.na(all_resulting_alphas)
    converged_indices <- original_indices[converged_mask]
    converged_alphas <- all_resulting_alphas[converged_mask]

    num_converged <- length(converged_alphas)

    # 默认返回值
    default_return_value <- list(
        selected_model_result = NA,
        num_in_consensus_cluster = 0,
        consensus_cluster_alphas = NA,
        consensus_cluster_initial_alphas = NA,
        all_initial_alphas = initial_alphas_used,
        all_resulting_alphas = all_resulting_alphas,
        all_convergence_status = all_convergence_status
    )

    if (num_converged == 0) {
        warning("没有优化尝试成功收敛。")
        return(default_return_value)
    }

    # 根据 permissible_error 确定四舍五入的精度
    # 例如: permissible_error = 1e-7, precision_digits = 7
    precision_digits <- max(0, -floor(log10(permissible_error)))

    rounded_converged_alphas <- round(converged_alphas, digits = precision_digits)

    # 处理只有一个收敛结果的情况
    if (num_converged == 1) {
        if (1 >= ceiling(number_of_attempts * consensus_proportion)) {
            return(list(
                selected_model_result = model_results[[converged_indices[1]]],
                num_in_consensus_cluster = 1,
                consensus_cluster_alphas = converged_alphas[1],
                consensus_cluster_initial_alphas = initial_alphas_used[converged_indices[1]],
                all_initial_alphas = initial_alphas_used,
                all_resulting_alphas = all_resulting_alphas,
                all_convergence_status = all_convergence_status
            ))
        } else {
            warning(paste0(
                "单个收敛结果 (alpha=", signif(converged_alphas[1], 4), ") 未达到共识比例要求 (需要至少 ",
                ceiling(number_of_attempts * consensus_proportion), " 个)。"
            ))
            return(default_return_value)
        }
    }

    # 找到出现频率最高的四舍五入后的 alpha 值 (众数)
    alpha_counts <- table(rounded_converged_alphas)
    # which.max 在有多个最大值时返回第一个的索引
    modal_rounded_alpha_str <- names(alpha_counts)[which.max(alpha_counts)]

    # 识别所有属于这个众数群组的原始收敛结果
    # indices_in_modal_group_relative 是相对于 converged_alphas/converged_indices 的索引
    indices_in_modal_group_relative <- which(as.character(rounded_converged_alphas) == modal_rounded_alpha_str)

    # 将这些相对索引映射回 model_results 中的原始索引
    original_indices_of_modal_group <- converged_indices[indices_in_modal_group_relative]

    alphas_in_modal_group <- converged_alphas[indices_in_modal_group_relative]
    num_in_modal_group <- length(alphas_in_modal_group)

    min_required_for_consensus <- ceiling(number_of_attempts * consensus_proportion)

    if (num_in_modal_group >= min_required_for_consensus) {
        # 达到共识！从此群组中随机抽取一个结果
        selected_original_index <- NA
        if (length(original_indices_of_modal_group) == 1) {
            selected_original_index <- original_indices_of_modal_group[1]
        } else {
            set.seed(Sys.time()) # 为了可重复性，可以在外部设置种子，或移除此行以获得纯随机
            selected_original_index <- sample(original_indices_of_modal_group, 1)
        }

        selected_model_full_result <- model_results[[selected_original_index]]

        return(list(
            selected_model_result = selected_model_full_result,
            num_in_consensus_cluster = num_in_modal_group,
            consensus_cluster_alphas = alphas_in_modal_group,
            consensus_cluster_initial_alphas = initial_alphas_used[original_indices_of_modal_group],
            all_initial_alphas = initial_alphas_used,
            all_resulting_alphas = all_resulting_alphas,
            all_convergence_status = all_convergence_status
        ))
    } else {
        warning(paste0(
            "基于四舍五入到 ", precision_digits, " 位小数，最密集的alpha集群 (围绕值约 ",
            signif(as.numeric(modal_rounded_alpha_str), 4), ") 包含 ", num_in_modal_group,
            " 个结果，未达到共识比例要求的至少 ", min_required_for_consensus, " 个。"
        ))
        return(default_return_value)
    }
}

# BIC 计算函数
bic_calculate <- function(
    model_cml, beta_hat_exp,
    beta_hat_out,
    beta_sigma_exp,
    beta_sigma_out,
    penalty_factor = 0.5,
    n_eff) {
    n_snps <- nrow(beta_hat_out)
    bic_term1 <- 0
    bic_term2 <- 0
    r_snp <- model_cml$r_snp
    k_snp <- sum((r_snp[, 1] != 0) == TRUE)
    # 计算拟合优度项 (与 -2 * logLikelihood 相关)
    for (i in 1:n_snps) {
        current_index <- (2 * i - 1):(2 * i)
        current_beta_hat_exp <- beta_hat_exp[i, ]
        current_beta_hat_out <- beta_hat_out[i, ]
        current_beta_exp <- model_cml$beta_exp[i, ]
        alpha <- model_cml$alpha
        current_r_snp <- model_cml$r_snp[i, ]
        current_beta_sigma_exp <- beta_sigma_exp[current_index, ]
        current_beta_sigma_out <- beta_sigma_out[current_index, ]

        bic_term1 <- t(current_beta_hat_exp - current_beta_exp) %*%
            current_beta_sigma_exp %*%
            (current_beta_hat_exp - current_beta_exp) +
            bic_term1

        bic_term2 <- t(current_beta_hat_out -
            alpha * current_beta_exp - current_r_snp) %*%
            current_beta_sigma_out %*%
            (current_beta_hat_out -
                alpha * current_beta_exp - current_r_snp) + bic_term2
    }


    penalty <- log(n_eff) * (1 + 2 * n_snps + 2 * k_snp) * penalty_factor

    # 计算总 BIC
    bic <- bic_term1 + bic_term2 + penalty

    return(as.numeric(bic)) # 确保返回标量数值
}
# 方差计算函数
## 输入的是模型
variance_calculation_function <- function(
    model_cml, beta_hat_out,
    beta_sigma_exp, beta_sigma_out) {
    snp_number <- nrow(beta_hat_out)
    k_c_snp_list <- 1:nrow(model_cml$r_snp)
    k_c_snp_list <- setdiff(k_c_snp_list, which(model_cml$r_snp[, 1] != 0))
    snp_number_valid <- length(k_c_snp_list)
    h_alpha_alpha <- matrix(0, ncol = 1, nrow = 1)
    h_alpha_gamma <- matrix(0, ncol = 1, nrow = 2 * snp_number_valid)
    h_gamma_gamma <- matrix(0, ncol = 2 * snp_number_valid, nrow = 2 * snp_number_valid)
    for (j in 1:length(k_c_snp_list)) {
        i <- k_c_snp_list[j]
        current_index_h <- (2 * j - 1):(2 * j)
        current_index <- (2 * i - 1):(2 * i)
        current_beta_hat_out <- beta_hat_out[i, ]
        current_beta_exp <- model_cml$beta_exp[i, ]
        alpha <- model_cml$alpha
        current_r_snp <- model_cml$r_snp[i, ]
        current_beta_sigma_exp <- beta_sigma_exp[current_index, ]
        current_beta_sigma_out <- beta_sigma_out[current_index, ]

        # 左上角
        h_alpha_alpha <- t(current_beta_exp) %*%
            current_beta_sigma_out %*% current_beta_exp + h_alpha_alpha

        # 两条相关的部分
        h_alpha_gamma_i <- current_beta_sigma_out %*%
            (2 * alpha * current_beta_exp -
                (current_beta_hat_out - current_r_snp))
        h_alpha_gamma[current_index_h] <- h_alpha_gamma_i

        # 右下角对角矩阵
        h_gamma_gamma_i <- current_beta_sigma_exp +
            as.numeric(alpha^2) * current_beta_sigma_out
        h_gamma_gamma[current_index_h, current_index_h] <- h_gamma_gamma_i
    }
    h <- matrix(0, ncol = 2 * snp_number_valid + 1, nrow = 2 * snp_number_valid + 1)
    h[1, 1] <- h_alpha_alpha
    h[2:(2 * snp_number_valid + 1), 1] <- h_alpha_gamma
    h[1, 2:(2 * snp_number_valid + 1)] <- t(h_alpha_gamma)
    h[2:(2 * snp_number_valid + 1), 2:(2 * snp_number_valid + 1)] <- h_gamma_gamma
    sd_alpha <- sqrt(solve(h)[1, 1])
    return(sd_alpha)
}

# %%
# 设计一个遍历函数

cml_family_ver_2_cpp <- function(
    beta_hat_exp, beta_hat_out,
    beta_sigma_exp, beta_sigma_out,
    permissible_error = 1e-7, # 更明确的默认容许误差
    maximum_number_of_iterations = 1000,
    number_of_attempts = 5,
    consensus_proportion = 0.5, n_eff = 1000) {
    snp_number <- nrow(beta_hat_exp)
    bic_score <- data.frame(
        k_snp = 0:(snp_number - 1),
        bic = rep(0, snp_number),
        alpha = rep(0, snp_number),
        alpha_sd = rep(0, snp_number)
    )
    for (i in 1:nrow(bic_score)) {
        current_model <- cml_multi_start_optimization_ver_2(
            beta_hat_exp, beta_hat_out,
            beta_sigma_exp, beta_sigma_out,
            permissible_error = permissible_error,
            maximum_number_of_iterations = maximum_number_of_iterations,
            number_of_attempts = number_of_attempts,
            consensus_proportion = consensus_proportion,
            illegal_snps = bic_score$k_snp[i]
        )
        # 加入收敛失败防爆

        if (all(current_model$all_convergence_status) == FALSE) {
            bic_score$bic[i] <- Inf
            bic_score$alpha[i] <- 0
            bic_score$alpha_sd[i] <- 10
            next()
        }
        current_bic <- bic_calculate(
            current_model$selected_model_result,
            beta_hat_exp,
            beta_hat_out,
            beta_sigma_exp,
            beta_sigma_out,
            penalty_factor = 0.5,
            n_eff = n_eff
        )
        alpha_var <- variance_calculation_function(
            current_model$selected_model_result, beta_hat_out,
            beta_sigma_exp, beta_sigma_out
        )
        bic_score$bic[i] <- current_bic
        bic_score$alpha[i] <- current_model$selected_model_result$alpha
        bic_score$alpha_sd[i] <- alpha_var
    }
    best_model_index <- arrange(bic_score, bic)[1, 1]
    best_model <- cml_multi_start_optimization_ver_2(
        beta_hat_exp, beta_hat_out,
        beta_sigma_exp, beta_sigma_out,
        permissible_error = permissible_error,
        maximum_number_of_iterations = maximum_number_of_iterations,
        number_of_attempts = number_of_attempts,
        consensus_proportion = consensus_proportion,
        illegal_snps = best_model_index
    )
    # 1. 计算原始权重
    # 权重公式为 exp(-bic/2)
    raw_weights <- exp(-bic_score$bic / 2)
    # 2. 计算 BIC 加权的 alpha
    bic_weighted_alpha <- weighted.mean(bic_score$alpha, raw_weights, na.rm = TRUE)
    # 3. 计算 BIC 加权的 alpha_sd (或 alpha_ad，取决于你的列名)
    bic_weighted_alpha_sd <- weighted.mean(bic_score$alpha_sd, raw_weights, na.rm = TRUE)
    results <- list(bic_weighted_alpha, bic_weighted_alpha_sd, best_model)
    return(results)
}




