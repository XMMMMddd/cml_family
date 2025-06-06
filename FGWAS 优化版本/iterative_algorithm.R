iterative_algorithm <- function(
    matrix_big, beta_sigma_exp, beta_sigma_out, beta_hat_exp, beta_hat_out, k,
    initial_alpha = 0, 
    initial_beta_exp = NULL, 
    initial_r = NULL,
    max_iter = 100, 
    tolerance = 1e-6
) {
    # 输入验证
    if (nrow(matrix_big) %% 4 != 0) stop("matrix_big must have row count divisible by 4")
    if (nrow(beta_hat_exp) != nrow(beta_hat_out)) stop("beta_hat_exp and beta_hat_out must have same row count")
    if (k < 0 || k > nrow(beta_hat_exp)) stop("k must be between 0 and number of SNPs (inclusive of n_snps for k=n_snps case)")
    if (max_iter <= 0) stop("max_iter must be a positive integer")
    if (tolerance <= 0) stop("tolerance must be a positive number")

    # 初始化变量
    n_snps <- nrow(matrix_big) / 4
    
    alpha <- initial_alpha
    
    if (is.null(initial_beta_exp)) {
        beta_exp <- beta_hat_exp 
    } else {
        if (nrow(initial_beta_exp) != n_snps || ncol(initial_beta_exp) != 2) {
            stop("initial_beta_exp has incorrect dimensions")
        }
        beta_exp <- initial_beta_exp
    }
    
    if (is.null(initial_r)) {
        r <- matrix(0, nrow = n_snps, ncol = 2)
    } else {
        if (nrow(initial_r) != n_snps || ncol(initial_r) != 2) {
            stop("initial_r has incorrect dimensions")
        }
        r <- initial_r
    }
    
    # 预计算大矩阵的逆矩阵 (4x4)
    matrix_big_inverses <- list()
    if (n_snps > 0) { # Avoid loop if n_snps is 0
        for (i in 1:n_snps) {
            matrix_index_all <- (4 * i - 3):(4 * i)
            matrix_big_i <- matrix_big[matrix_index_all, ]
            matrix_big_inverses[[i]] <- tryCatch(
                solve(matrix_big_i),
                error = function(e) {
                    warning(paste("Singular matrix at SNP", i, "- using pseudoinverse. Consider checking matrix_big inputs."))
                    if (!requireNamespace("MASS", quietly = TRUE)) {
                        stop("Package 'MASS' needed for ginv() when matrix is singular. Please install it.", call. = FALSE)
                    }
                    MASS::ginv(matrix_big_i)
                }
            )
        }
    }
    
    iterations_taken <- 0
    converged <- FALSE
    alpha_old <- alpha # Initialize alpha_old before loop for first convergence check

    if (n_snps > 0) { # Only iterate if there are SNPs
        for (iter_count in 1:max_iter) {
            iterations_taken <- iter_count
            alpha_old <- alpha # Store alpha from previous iteration (or initial for first)
            
            # 1. 更新r
            d <- numeric(n_snps)
            r_new_iter <- matrix(0, nrow = n_snps, ncol = 2)
            
            for (i in 1:n_snps) {
                matrix_index_2x2 <- (2 * i - 1):(2 * i)
                
                beta_hat_exp_i <- beta_hat_exp[i, , drop = FALSE]
                beta_hat_out_i <- beta_hat_out[i, , drop = FALSE]
                beta_exp_i <- beta_exp[i, , drop = FALSE]
                beta_sigma_exp_i <- beta_sigma_exp[matrix_index_2x2, , drop = FALSE]
                beta_sigma_out_i <- beta_sigma_out[matrix_index_2x2, , drop = FALSE]
                beta_sigma_rho_i <- matrix_big[(4*i-3):(4*i-2), 3:4, drop = FALSE]
                
                residual <- beta_hat_out_i - alpha * beta_exp_i
                d[i] <- residual %*% beta_sigma_out_i %*% t(residual)
                
                exposure_residual <- beta_hat_exp_i - beta_exp_i
                r_i_calc <- residual - exposure_residual %*% t(beta_sigma_exp_i) %*% beta_sigma_rho_i
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

            # 2. 更新beta_exp
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

            # 3. 更新alpha
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
            
            current_alpha_val <- alpha # Store current alpha before potential update
            if (abs(fenmu_alpha_total) < .Machine$double.eps) {
                warning(paste("Denominator for alpha is close to zero in iteration", iter_count, ". Alpha not updated this iteration."))
            } else {
                alpha <- as.numeric(fenzi_alpha_total / fenmu_alpha_total)
            }
            
            # 检查收敛性 (compare newly computed alpha with alpha_old from start of this iteration)
            if (iter_count > 0 && abs(alpha - alpha_old) < tolerance) { # iter_count > 0 condition is redundant here
                converged <- TRUE
                break
            }
        } # End of main iteration loop
    } else { # n_snps is 0
        # If no SNPs, algorithm doesn't run iterations.
        # Alpha remains initial_alpha. Beta_exp and r remain initial or default empty.
        # Convergence is trivially true if max_iter is 0, or not applicable.
        # Let's set converged to TRUE if no iterations are run.
        if (max_iter == 0) converged <- TRUE 
        # Or, consider what converged means if n_snps = 0.
        # For now, if n_snps = 0, iterations_taken = 0, alpha = initial_alpha.
        # converged can be set to TRUE as there's no change possible.
        converged <- TRUE 
    }
    
    if (iterations_taken == max_iter && !converged) {
        # Simplified warning message construction
        final_alpha_change <- abs(alpha - alpha_old) # alpha_old here is from the start of the last iteration
        warning_message <- paste("Algorithm did not converge within", max_iter, 
                                 "iterations. Last alpha change:", format(final_alpha_change, scientific = FALSE, digits = max(3, getOption("digits") - 3 )))
        warning(warning_message)
    }
    
    # 返回结果
    return(list(
        alpha = alpha,
        beta_exp = beta_exp,
        r = r,
        iterations = iterations_taken,
        converged = converged
    ))
}
