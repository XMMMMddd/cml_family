// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <algorithm>
#include <cmath>

struct SortItem {
    double value;
    arma::uword index;

    bool operator<(const SortItem& other) const {
        return value > other.value;
    }
};


// [[Rcpp::export]]
Rcpp::List cMl_para_yuanshi_cpp(
    const arma::mat& beta_y_hat,
    const arma::mat& beta_x_hat,
    const arma::mat& Sigma_inv_x,
    const arma::mat& Sigma_inv_y,
    Rcpp::Nullable<double> theta_init_nullable = R_NilValue,           // Optional initial theta
    const Rcpp::Nullable<arma::mat>& beta_x_init_nullable = R_NilValue, // Optional initial beta_x
    int k_par = 0,
    int k_x = 0,
    int max_iterations = 1000,
    double value_diff = 1e-9) {

    arma::uword n_rows = beta_y_hat.n_rows;
    if (n_rows == 0) {
        Rcpp::stop("Input beta_y_hat has zero rows.");
    }

    if (beta_y_hat.n_cols != 2 || beta_x_hat.n_cols != 2 ||
        beta_x_hat.n_rows != n_rows ||
        Sigma_inv_x.n_rows != 2 * n_rows || Sigma_inv_x.n_cols != 2 ||
        Sigma_inv_y.n_rows != 2 * n_rows || Sigma_inv_y.n_cols != 2) {
        Rcpp::stop("Input matrix dimensions are inconsistent.");
    }

     if (k_par < 0 || static_cast<arma::uword>(k_par) > n_rows || k_x < 0 || static_cast<arma::uword>(k_x) > n_rows) {
         Rcpp::stop("k_par or k_x is out of valid range [0, n_rows].");
    }

    arma::vec r_x_vec(n_rows, arma::fill::zeros);
    arma::vec r_p_vec(n_rows, arma::fill::zeros);
    arma::mat r = arma::join_rows(r_x_vec, r_p_vec);

    double theta;
    if (theta_init_nullable.isNotNull()) {
        theta = Rcpp::as<double>(theta_init_nullable);
    } else {
        theta = 0.0; // Default initial value
    }

    arma::mat beta_x;
    if (beta_x_init_nullable.isNotNull()) {
        arma::mat beta_x_init = Rcpp::as<arma::mat>(beta_x_init_nullable);
        if (beta_x_init.n_rows != n_rows || beta_x_init.n_cols != 2) {
             Rcpp::stop("Provided beta_x_init has incorrect dimensions. Expected %u rows and 2 columns.", n_rows);
        }
        beta_x = beta_x_init;
    } else {
        beta_x = beta_x_hat; // Default initial value
    }


    double theta_new = theta;
    double theta_old = theta_new + value_diff * 10.0; // Ensure initial difference

    arma::mat beta_x_old = beta_x;
    // Initialize beta_x_old slightly differently to ensure the loop runs at least once if initial values are already converged
    if (n_rows > 0 && beta_x.n_elem > 0) {
         beta_x_old(0) = beta_x(0) + value_diff * 10.0;
    }


    int n_times = 0;
    bool converged = false;

    arma::vec r_x_current(n_rows);
    arma::vec r_p_current(n_rows);
    arma::vec A_values(n_rows);
    arma::mat beta_x_cal(n_rows, 2);
    arma::mat Sigma_inv_y_i(2, 2);
    arma::mat Sigma_inv_x_i(2, 2);
    arma::rowvec diff_vec(2);
    std::vector<SortItem> sort_items(n_rows);
    arma::vec r_x_new_final(n_rows, arma::fill::zeros);
    arma::vec r_p_new_final(n_rows, arma::fill::zeros);
    arma::mat beta_x_new(n_rows, 2);
    arma::rowvec fenzi_beta_row(2);
    arma::mat fenmu_beta_mat(2, 2);
    arma::mat inv_fenmu_beta;
    double fenzi_theta_sum = 0.0;
    double fenmu_theta_sum = 0.0;
    arma::rowvec beta_y_hat_i_row(2);
    arma::rowvec r_i_row(2);
    arma::rowvec beta_x_i_row(2);
    arma::mat Sigma_y_i(2, 2);
    bool success_inv;


    while (n_times < max_iterations) {

        double beta_diff = arma::accu(arma::abs(beta_x_old - beta_x));
        double theta_diff = std::abs(theta_old - theta_new);

        // Check convergence at the beginning of the loop (after potential first iteration)
        if (n_times > 0 && theta_diff <= value_diff && beta_diff <= value_diff) {
             converged = true;
             break;
        }
         // Special case: if started exactly at converged state, run one iter anyway? Or exit?
         // Current logic exits after 0 iterations if already converged. Let's ensure at least one calculation if needed.
         // The initialization of theta_old/beta_x_old should prevent immediate exit unless diffs are truly zero.

        if (n_times > 0) { // Only update old values after the first check
             beta_x_old = beta_x;
             theta_old = theta_new; // theta_new is the result from the *previous* iteration
        }
        // Update theta for the current iteration's calculations *after* storing the old value
        theta = theta_new;


        n_times++;


        r_p_current = -(theta * beta_x.col(1)) + beta_y_hat.col(1);
        r_x_current.zeros();

        arma::vec beta_x_cal_x = theta * beta_x.col(0) + r_x_current;
        arma::vec beta_x_cal_p = theta * beta_x.col(1) + r_p_current;
        beta_x_cal = arma::join_rows(beta_x_cal_x, beta_x_cal_p);

        for (arma::uword i = 0; i < n_rows; ++i) {
            Sigma_inv_y_i = Sigma_inv_y.rows(2 * i, 2 * i + 1);
            diff_vec = beta_y_hat.row(i) - beta_x_cal.row(i);
            A_values(i) = arma::as_scalar(diff_vec * Sigma_inv_y_i * diff_vec.t());
            sort_items[i] = {A_values(i), i};
        }

        std::sort(sort_items.begin(), sort_items.end());

        arma::vec r_x_potential_new = -(theta * beta_x.col(0)) + beta_y_hat.col(0);

        r_x_new_final.zeros();
        if (k_x > 0) {
            for (int j = 0; j < k_x; ++j) {
                arma::uword original_index = sort_items[j].index;
                r_x_new_final(original_index) = r_x_potential_new(original_index);
            }
        }


         r_x_current = -(theta * beta_x.col(0)) + beta_y_hat.col(0);
         r_p_current.zeros();

        beta_x_cal_x = theta * beta_x.col(0) + r_x_current;
        beta_x_cal_p = theta * beta_x.col(1) + r_p_current;
        beta_x_cal = arma::join_rows(beta_x_cal_x, beta_x_cal_p);

        for (arma::uword i = 0; i < n_rows; ++i) {
            Sigma_inv_y_i = Sigma_inv_y.rows(2 * i, 2 * i + 1);
            diff_vec = beta_y_hat.row(i) - beta_x_cal.row(i);
            A_values(i) = arma::as_scalar(diff_vec * Sigma_inv_y_i * diff_vec.t());
            sort_items[i] = {A_values(i), i};
        }

        std::sort(sort_items.begin(), sort_items.end());

        arma::vec r_p_potential_new = -(theta * beta_x.col(1)) + beta_y_hat.col(1);

        r_p_new_final.zeros();
        if (k_par > 0) {
            for (int j = 0; j < k_par; ++j) {
                arma::uword original_index = sort_items[j].index;
                r_p_new_final(original_index) = r_p_potential_new(original_index);
            }
        }

        r = arma::join_rows(r_x_new_final, r_p_new_final);


        for (arma::uword i = 0; i < n_rows; ++i) {
            Sigma_inv_x_i = Sigma_inv_x.rows(2 * i, 2 * i + 1);
            Sigma_inv_y_i = Sigma_inv_y.rows(2 * i, 2 * i + 1);

            fenzi_beta_row = beta_x_hat.row(i) * Sigma_inv_x_i +
                             theta * beta_y_hat.row(i) * Sigma_inv_y_i -
                             theta * r.row(i) * Sigma_inv_y_i;

            fenmu_beta_mat = Sigma_inv_x_i + std::pow(theta, 2) * Sigma_inv_y_i;

            success_inv = arma::inv_sympd(inv_fenmu_beta, fenmu_beta_mat);
             if (!success_inv) {
                  success_inv = arma::inv(inv_fenmu_beta, fenmu_beta_mat);
                  if (!success_inv) {
                      Rcpp::warning("Matrix inversion failed for beta_x update at index %u. Using pseudo-inverse.", i + 1);
                      inv_fenmu_beta = arma::pinv(fenmu_beta_mat);
                  }
            }
            beta_x_new.row(i) = fenzi_beta_row * inv_fenmu_beta;
        }
        // Update beta_x *before* theta calculation for the current iteration
        beta_x = beta_x_new;


        fenzi_theta_sum = 0.0;
        fenmu_theta_sum = 0.0;

        for (arma::uword i = 0; i < n_rows; ++i) {
            beta_y_hat_i_row = beta_y_hat.row(i);
            r_i_row = r.row(i);
            Sigma_inv_y_i = Sigma_inv_y.rows(2 * i, 2 * i + 1);
            beta_x_i_row = beta_x.row(i); // Use the just updated beta_x


            success_inv = arma::inv_sympd(Sigma_y_i, Sigma_inv_y_i);
             if (!success_inv) {
                  success_inv = arma::inv(Sigma_y_i, Sigma_inv_y_i);
                   if (!success_inv) {
                        Rcpp::warning("Matrix inversion failed for Sigma_inv_y at index %u in theta update. Skipping contribution.", i + 1);
                        continue;
                   }
            }


            if (r_i_row(0) == 0.0 && r_i_row(1) == 0.0) {
                fenzi_theta_sum += arma::as_scalar((beta_y_hat_i_row - r_i_row) * Sigma_inv_y_i * beta_x_i_row.t());
                fenmu_theta_sum += arma::as_scalar(beta_x_i_row * Sigma_inv_y_i * beta_x_i_row.t());
            } else if (r_i_row(0) != 0.0 && r_i_row(1) == 0.0) {
                double sigma_y_22 = Sigma_y_i(1, 1);
                if (std::abs(sigma_y_22) > 1e-10) {
                    fenzi_theta_sum += beta_y_hat_i_row(1) / sigma_y_22 * beta_x_i_row(1);
                    fenmu_theta_sum += beta_x_i_row(1) / sigma_y_22 * beta_x_i_row(1);
                }
            } else if (r_i_row(0) == 0.0 && r_i_row(1) != 0.0) {
                double sigma_y_11 = Sigma_y_i(0, 0);
                 if (std::abs(sigma_y_11) > 1e-10) {
                    fenzi_theta_sum += beta_y_hat_i_row(0) / sigma_y_11 * beta_x_i_row(0);
                    fenmu_theta_sum += beta_x_i_row(0) / sigma_y_11 * beta_x_i_row(0);
                 }
            }

        }

        // Calculate the candidate for the *next* iteration's theta
        if (std::abs(fenmu_theta_sum) > 1e-10) {
            theta_new = fenzi_theta_sum / fenmu_theta_sum;
        } else {
            Rcpp::warning("Denominator for theta update is close to zero. Theta (%f) remains unchanged for the next iteration.", theta);
            theta_new = theta; // Keep theta the same if denom is zero
        }

    } // End while loop


    if (n_times == max_iterations && !converged) {
        Rcpp::warning("Maximum number of iterations (%d) reached. Convergence may not have been achieved.", max_iterations);
    }


    arma::mat H_1(2, 2, arma::fill::zeros);
    arma::mat H_2(2, 1, arma::fill::zeros);
    double H_3_scalar = 0.0;
    arma::mat Sigma_inv_y_i_orig(2, 2);
    arma::mat Sigma_inv_y_i_eff(2, 2);
    arma::vec beta_x_i_eff_vec(2);

    // Use the final converged theta for SE calculation
    double final_theta = theta_new;


    for (arma::uword i = 0; i < n_rows; ++i) {
        r_i_row = r.row(i); // Use the final r
        beta_x_i_row = beta_x.row(i); // Use the final beta_x
        Sigma_inv_y_i_orig = Sigma_inv_y.rows(2 * i, 2 * i + 1);
        Sigma_inv_x_i = Sigma_inv_x.rows(2 * i, 2 * i + 1);


        Sigma_inv_y_i_eff = Sigma_inv_y_i_orig;
        beta_x_i_eff_vec = beta_x_i_row.t();

        if (r_i_row(0) == 0.0 && r_i_row(1) == 0.0) {
            H_1 += Sigma_inv_x_i + std::pow(final_theta, 2) * Sigma_inv_y_i_eff;
            H_2 += 2.0 * final_theta * Sigma_inv_y_i_eff * beta_x_i_eff_vec;
            H_3_scalar += arma::as_scalar(beta_x_i_eff_vec.t() * Sigma_inv_y_i_eff * beta_x_i_eff_vec);
        } else {
             success_inv = arma::inv_sympd(Sigma_y_i, Sigma_inv_y_i_eff);
             if (!success_inv) {
                  success_inv = arma::inv(Sigma_y_i, Sigma_inv_y_i_eff);
                  if (!success_inv) {
                       Rcpp::warning("Matrix inversion failed for Sigma_inv_y at index %u in SE calculation. Skipping contribution.", i + 1);
                        if (r_i_row(0) != 0.0 && r_i_row(1) != 0.0){
                             H_1 += Sigma_inv_x_i;
                        }
                       continue;
                  }
             }

             if (r_i_row(0) != 0.0 && r_i_row(1) == 0.0) {
                 double sigma_y_22 = Sigma_y_i(1, 1);
                 if (std::abs(sigma_y_22) > 1e-10) {
                     H_3_scalar += std::pow(beta_x_i_eff_vec(1), 2) / sigma_y_22;
                 }
                 H_1 += Sigma_inv_x_i;
             } else if (r_i_row(0) == 0.0 && r_i_row(1) != 0.0) {
                 double sigma_y_11 = Sigma_y_i(0, 0);
                 if (std::abs(sigma_y_11) > 1e-10) {
                     H_3_scalar += std::pow(beta_x_i_eff_vec(0), 2) / sigma_y_11;
                 }
                 H_1 += Sigma_inv_x_i;
             } else if (r_i_row(0) != 0.0 && r_i_row(1) != 0.0) {
                 H_1 += Sigma_inv_x_i;
             }
        }
    }


    arma::mat H_matrix(3, 3);
    H_matrix.submat(0, 0, 1, 1) = H_1;
    H_matrix.submat(0, 2, 1, 2) = H_2;
    H_matrix.submat(2, 0, 2, 1) = H_2.t();
    H_matrix(2, 2) = H_3_scalar;

    arma::mat H_inv;
    double theta_se = NA_REAL;

    success_inv = arma::inv_sympd(H_inv, H_matrix);
    if (!success_inv) {
        success_inv = arma::inv(H_inv, H_matrix);
    }
    if (!success_inv) {
         Rcpp::warning("Hessian matrix inversion failed for SE calculation. SE will be NA.");
    } else {
        if (H_inv(2, 2) < 0) {
             Rcpp::warning("Variance estimate for theta is negative (%g). SE will be NA.", H_inv(2,2));
        } else {
            theta_se = std::sqrt(H_inv(2, 2));
        }
    }



    return Rcpp::List::create(
        Rcpp::Named("theta") = final_theta, // Return the final converged value
        Rcpp::Named("beta_x") = beta_x,
        Rcpp::Named("r") = r,
        Rcpp::Named("theta_se") = theta_se,
        Rcpp::Named("iterations") = n_times,
        Rcpp::Named("converged") = converged
    );
}