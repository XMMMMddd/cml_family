#include <RcppArmadillo.h>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List iterative_algorithm_rcpp(
    const arma::mat& matrix_big,
    const arma::mat& beta_sigma_exp,
    const arma::mat& beta_sigma_out,
    const arma::mat& beta_hat_exp,
    const arma::mat& beta_hat_out,
    int k,
    double initial_alpha = 0.0,
    Rcpp::Nullable<Rcpp::NumericMatrix> initial_beta_exp_r = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> initial_r_r = R_NilValue,
    int max_iter = 100,
    double tolerance = 1e-6) {

    if (matrix_big.n_rows % 4 != 0) {
        Rcpp::stop("参数 'matrix_big' 的行数必须是 4 的倍数。");
    }
    if (beta_hat_exp.n_rows != beta_hat_out.n_rows) {
        Rcpp::stop("参数 'beta_hat_exp' 和 'beta_hat_out' 的行数必须一致。");
    }
    if (k < 0 || k > static_cast<int>(beta_hat_exp.n_rows)) {
        Rcpp::stop("参数 'k' 必须在 0 和 SNP 总数之间 (包含 0 和 SNP 总数)。");
    }
    if (max_iter <= 0) {
        Rcpp::stop("参数 'max_iter' 必须是一个正整数。");
    }
    if (tolerance <= 0) {
        Rcpp::stop("参数 'tolerance' 必须是一个正数。");
    }

    int n_snps = matrix_big.n_rows / 4;
    double alpha = initial_alpha;
    arma::mat beta_exp;
    if (initial_beta_exp_r.isNotNull()) {
        beta_exp = Rcpp::as<arma::mat>(initial_beta_exp_r);
    } else {
        beta_exp = beta_hat_exp;
    }

    arma::mat r_arma;
    if (initial_r_r.isNotNull()) {
        r_arma = Rcpp::as<arma::mat>(initial_r_r);
    } else {
        if (n_snps > 0) {
            r_arma.zeros(n_snps, 2);
        } else {
            r_arma.set_size(0, 2);
        }
    }

    std::vector<arma::mat> matrix_big_inverses(n_snps);
    if (n_snps > 0) {
        for (int i = 0; i < n_snps; ++i) {
            arma::mat matrix_big_i = matrix_big.rows(4 * i, 4 * i + 3);
            arma::mat inv_matrix_big_i;
            bool success = arma::inv(inv_matrix_big_i, matrix_big_i);
            if (!success) {
                Rcpp::warning("Rcpp version: 4x4 matrix at SNP %d is singular - using pseudo-inverse.", i + 1);
                inv_matrix_big_i = arma::pinv(matrix_big_i);
            }
            matrix_big_inverses[i] = inv_matrix_big_i;
        }
    }

    int iterations_taken = 0;
    bool converged = false;
    double alpha_old = alpha;

    if (n_snps > 0) {
        for (int iter_count = 1; iter_count <= max_iter; ++iter_count) {
            iterations_taken = iter_count;
            alpha_old = alpha;
            arma::vec d(n_snps);
            arma::mat r_new_iter(n_snps, 2);

            for (int i = 0; i < n_snps; ++i) {
                arma::rowvec beta_hat_exp_i = beta_hat_exp.row(i);
                arma::rowvec beta_hat_out_i = beta_hat_out.row(i);
                arma::rowvec beta_exp_i = beta_exp.row(i);
                arma::mat beta_sigma_exp_i = beta_sigma_exp.rows(2 * i, 2 * i + 1);
                arma::mat beta_sigma_out_i = beta_sigma_out.rows(2 * i, 2 * i + 1);
                arma::mat beta_sigma_rho_i = matrix_big.submat(4 * i, 2, 4 * i + 1, 3);
                
                arma::rowvec residual_d = beta_hat_out_i - alpha * beta_exp_i;
                d(i) = arma::as_scalar(residual_d * beta_sigma_out_i * residual_d.t());
                
                arma::rowvec exposure_residual_r = beta_hat_exp_i - beta_exp_i;
                arma::rowvec r_i_calc = residual_d - exposure_residual_r * beta_sigma_exp_i.t() * beta_sigma_rho_i;
                r_new_iter.row(i) = r_i_calc;
            }

            arma::uvec d_order_indices = arma::sort_index(d, "descend");
            arma::mat r_new_sorted_temp = r_new_iter.rows(d_order_indices);

            if (k < n_snps) {
                for (int j = k; j < n_snps; ++j) {
                    r_new_sorted_temp.row(j).zeros();
                }
            }
            
            r_arma.zeros(n_snps, 2);
            for (unsigned int j = 0; j < (unsigned int)n_snps; ++j) {
                 r_arma.row(d_order_indices(j)) = r_new_sorted_temp.row(j);
            }


            arma::mat beta_exp_new_iter(n_snps, 2);
            for (int i = 0; i < n_snps; ++i) {
                arma::mat omega_11 = matrix_big_inverses[i].submat(0, 0, 1, 1);
                arma::mat omega_12 = matrix_big_inverses[i].submat(0, 2, 1, 3);
                arma::mat omega_21 = matrix_big_inverses[i].submat(2, 0, 3, 1);
                arma::mat omega_22 = matrix_big_inverses[i].submat(2, 2, 3, 3);
                
                arma::rowvec r_i_current = r_arma.row(i);
                
                arma::mat fenmu_beta = omega_11 + alpha * (omega_12 + omega_21) + alpha * alpha * omega_22;
                arma::mat fenzi_beta_term1 = (omega_11 + alpha * omega_21) * beta_hat_exp.row(i).t();
                arma::mat fenzi_beta_term2 = (omega_12 + alpha * omega_22) * (beta_hat_out.row(i) - r_i_current).t();
                arma::mat fenzi_beta = fenzi_beta_term1 + fenzi_beta_term2;
                
                beta_exp_new_iter.row(i) = arma::solve(fenmu_beta, fenzi_beta).t();
            }
            beta_exp = beta_exp_new_iter;

            double fenzi_alpha_total = 0;
            double fenmu_alpha_total = 0;

            for (int i = 0; i < n_snps; ++i) {
                arma::mat omega_12 = matrix_big_inverses[i].submat(0, 2, 1, 3);
                arma::mat omega_22 = matrix_big_inverses[i].submat(2, 2, 3, 3);
                
                arma::rowvec beta_hat_exp_i = beta_hat_exp.row(i);
                arma::rowvec beta_exp_i_current = beta_exp.row(i);
                arma::rowvec beta_hat_out_i = beta_hat_out.row(i);
                arma::rowvec r_i_current = r_arma.row(i);
                
                double fenzi_term_alpha = arma::as_scalar(
                    (beta_hat_exp_i - beta_exp_i_current) * omega_12 * beta_exp_i_current.t() +
                    beta_exp_i_current * omega_22 * (beta_hat_out_i - r_i_current).t()
                );
                double fenmu_term_alpha = arma::as_scalar(
                    beta_exp_i_current * omega_22 * beta_exp_i_current.t()
                );
                
                fenzi_alpha_total += fenzi_term_alpha;
                fenmu_alpha_total += fenmu_term_alpha;
            }

            if (std::abs(fenmu_alpha_total) < std::numeric_limits<double>::epsilon() * 100) { // Increased tolerance slightly from machine epsilon
                 Rcpp::warning("Rcpp version: Alpha denominator near zero in iteration %d.", iter_count);
            } else {
                alpha = fenzi_alpha_total / fenmu_alpha_total;
            }

            if (std::abs(alpha - alpha_old) < tolerance) {
                converged = true;
                break;
            }
        }
    } else {
        converged = true;
        iterations_taken = 0;
    }

    if (iterations_taken >= max_iter && !converged) {
        double final_alpha_change = std::abs(alpha - alpha_old);
         Rcpp::warning("Rcpp version: Algorithm did not converge within %d iterations. Final alpha change: %f", max_iter, final_alpha_change);
    }
    
    return Rcpp::List::create(
        Rcpp::Named("alpha") = alpha,
        Rcpp::Named("beta_exp") = beta_exp,
        Rcpp::Named("r") = r_arma,
        Rcpp::Named("iterations") = iterations_taken,
        Rcpp::Named("converged") = converged
    );
}