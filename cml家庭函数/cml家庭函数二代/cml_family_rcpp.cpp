#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List cml_family_ver_2_rcpp(
    arma::mat beta_hat_exp,
    arma::mat beta_hat_out,
    arma::mat beta_sigma_exp_stack,
    arma::mat beta_sigma_out_stack,
    double initial_alpha = 0.0, // 新增：用户自定义的alpha起始点
    double permissible_error = 0.0000001,
    int maximum_number_of_iterations = 1000,
    int illegal_snps_arg = 0) {

    if (beta_hat_exp.n_cols != 2) {
        Rcpp::stop("beta_hat_exp 必须有 2 列。");
    }
    if (beta_hat_out.n_cols != 2) {
        Rcpp::stop("beta_hat_out 必须有 2 列。");
    }
    int snp_number = beta_hat_exp.n_rows;
    if (beta_hat_out.n_rows != static_cast<arma::uword>(snp_number)) {
        Rcpp::stop("beta_hat_exp 和 beta_hat_out 必须有相同的行数 (SNP 数量)。");
    }

    if (beta_sigma_exp_stack.n_cols != 2) {
        Rcpp::stop("beta_sigma_exp_stack 必须有 2 列。");
    }
    if (beta_sigma_out_stack.n_cols != 2) {
        Rcpp::stop("beta_sigma_out_stack 必须有 2 列。");
    }
    if (beta_sigma_exp_stack.n_rows != static_cast<arma::uword>(2 * snp_number)) {
        Rcpp::stop("beta_sigma_exp_stack 必须有 2*snp_number 行。");
    }
    if (beta_sigma_out_stack.n_rows != static_cast<arma::uword>(2 * snp_number)) {
        Rcpp::stop("beta_sigma_out_stack 必须有 2*snp_number 行。");
    }

    double alpha = initial_alpha; // 使用用户提供的或默认的 initial_alpha
    arma::mat r_snp( (snp_number > 0 ? snp_number : 0), 2, arma::fill::zeros);
    arma::mat beta_exp = beta_hat_exp;

    double alpha_last = alpha + 1000.0;
    bool R_style_loop_condition = true;
    int number_of_iterations = 0;

    while (R_style_loop_condition) {
        number_of_iterations++;

        bool check_for_next_iter;
        if (number_of_iterations == 1 && (alpha_last == alpha + 1000.0)) { // 此处的alpha是initial_alpha
             check_for_next_iter = (std::abs( (alpha + 1000.0) - alpha) > permissible_error &&
                                    number_of_iterations < maximum_number_of_iterations);
        } else {
             check_for_next_iter = (std::abs(alpha_last - alpha) > permissible_error &&
                                    number_of_iterations < maximum_number_of_iterations);
        }

        alpha_last = alpha;

        if (snp_number > 0) {
            r_snp.zeros();
        }

        arma::mat beta_exp_new = beta_exp;
        if (snp_number > 0) {
            for (int i = 0; i < snp_number; ++i) {
                arma::rowvec current_beta_hat_exp_row = beta_hat_exp.row(i);
                arma::rowvec current_beta_hat_out_row_val = beta_hat_out.row(i);
                arma::mat current_beta_sigma_exp_mat = beta_sigma_exp_stack.rows(2 * i, 2 * i + 1);
                arma::mat current_beta_sigma_out_mat = beta_sigma_out_stack.rows(2 * i, 2 * i + 1);
                arma::rowvec r_snp_i_zero = r_snp.row(i);

                arma::mat inv_denom_beta_exp = arma::inv(current_beta_sigma_exp_mat + std::pow(alpha, 2) * current_beta_sigma_out_mat);

                arma::colvec term1_num_beta_exp = current_beta_sigma_exp_mat * current_beta_hat_exp_row.t();
                arma::colvec term2_num_beta_exp_factor = current_beta_sigma_out_mat * (current_beta_hat_out_row_val - r_snp_i_zero).t();
                arma::colvec num_beta_exp = term1_num_beta_exp + alpha * term2_num_beta_exp_factor;

                beta_exp_new.row(i) = (inv_denom_beta_exp * num_beta_exp).t();
            }
            beta_exp = beta_exp_new;
        }

        double num_alpha_update = 0.0;
        double den_alpha_update = 0.0;
        if (snp_number > 0) {
            for (int i = 0; i < snp_number; ++i) {
                arma::rowvec current_beta_exp_row = beta_exp.row(i);
                arma::rowvec current_beta_hat_out_row_val = beta_hat_out.row(i);
                arma::mat current_beta_sigma_out_mat = beta_sigma_out_stack.rows(2 * i, 2 * i + 1);
                arma::rowvec current_r_snp_zero_row = r_snp.row(i);

                num_alpha_update += arma::as_scalar((current_beta_hat_out_row_val - current_r_snp_zero_row) *
                                             current_beta_sigma_out_mat * current_beta_exp_row.t());
                den_alpha_update += arma::as_scalar(current_beta_exp_row *
                                             current_beta_sigma_out_mat * current_beta_exp_row.t());
            }
        }

        if (snp_number == 0) {
             alpha = initial_alpha; // 如果没有SNP，alpha保持初始值或根据迭代逻辑变为0
             if (number_of_iterations > 0 && den_alpha_update == 0.0) alpha = 0.0; // 避免NaN，如果迭代过且分母为0
        } else if (std::abs(den_alpha_update) < 1e-15) {
            Rcpp::warning("alpha 更新的分母在迭代 %d 时接近于零。alpha 被设为 0。结果可能不稳定。", number_of_iterations);
            alpha = 0.0;
        } else {
            alpha = num_alpha_update / den_alpha_update;
        }

        R_style_loop_condition = check_for_next_iter;
    }

    bool convergence_met;
    if (number_of_iterations >= maximum_number_of_iterations &&
        std::abs(alpha_last - alpha) > permissible_error) {
        convergence_met = false;
    } else {
        convergence_met = true;
    }

    if (snp_number == 0) {
        if (number_of_iterations == 0 && maximum_number_of_iterations == 0) {
             convergence_met = (std::abs(alpha_last - alpha) <= permissible_error);
        } else if (number_of_iterations == 0 && maximum_number_of_iterations > 0) {
             convergence_met = false;
        }
    }

    return Rcpp::List::create(
        Rcpp::Named("alpha") = alpha,
        Rcpp::Named("beta_exp") = beta_exp,
        Rcpp::Named("r_snp") = r_snp,
        Rcpp::Named("convergence") = convergence_met,
        Rcpp::Named("iterations") = number_of_iterations
    );
}