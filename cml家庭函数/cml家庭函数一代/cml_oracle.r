oracle_function <- function(beta_x_hat, Sigma_inv_x, beta_y_hat, Sigma_inv_y) {
    # beta_y_hat
    fenzi <- 0
    fenmu <- 0
    theta_se <- 0
    for (i in 1:(nrow(Sigma_inv_y) / 2)) {
        fenzi <- t(beta_y_hat[i, ]) %*% Sigma_inv_y[(2 * i - 1):(2 * i), ] %*% beta_x_hat[i, ] + fenzi
        fenmu <- t(beta_x_hat[i, ]) %*% Sigma_inv_y[(2 * i - 1):(2 * i), ] %*% beta_x_hat[i, ] + fenmu
        theta_se <- t(beta_x_hat[i, ]) %*% Sigma_inv_y[(2 * i - 1):(2 * i), ] %*% beta_x_hat[i, ] + theta_se
    }

    theta <- fenzi / fenmu
    theta_se <- sqrt(1 / theta_se)

    return(c(theta, theta_se))
}
