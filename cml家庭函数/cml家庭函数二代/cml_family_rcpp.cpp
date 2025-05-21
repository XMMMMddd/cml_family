// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::export]]
Rcpp::List cml_family_ver_2_rcpp(
    const arma::mat& beta_hat_exp,
    const arma::mat& beta_hat_out,
    const arma::mat& beta_sigma_exp,
    const arma::mat& beta_sigma_out,
    double permissible_error = 0.0000001,
    int maximum_number_of_iterations = 1000,
    int illegal_snps = 0,
    double initial_alpha_param = 1.0) {

    int snp_number = beta_hat_exp.n_rows;
    double alpha = 0.0;
    arma::mat r_snp = arma::zeros<arma::mat>(snp_number, 2);
    arma::mat beta_exp = beta_hat_exp;
    double alpha_last = alpha + 1000.0;
    bool termination_conditions_met_for_loop = true;
    int number_of_iterations = 0;

    while (termination_conditions_met_for_loop) {
        number_of_iterations++;
        bool should_continue_next_iteration;
        if (number_of_iterations < maximum_number_of_iterations) {
            if (std::abs(alpha_last - alpha) > permissible_error) {
                should_continue_next_iteration = true;
            } else {
                should_continue_next_iteration = false;
            }
        } else {
            should_continue_next_iteration = false;
        }

        alpha_last = alpha;

        arma::vec diff_snp_scalar_vec(snp_number, arma::fill::zeros);
        arma::mat r_snp_intermediate_diff_betas(snp_number, 2, arma::fill::zeros);

        if (snp_number > 0) {
            for (int i = 0; i < snp_number; ++i) {
                arma::rowvec current_beta_exp_i = beta_exp.row(i);
                arma::rowvec current_beta_hat_out_i = beta_hat_out.row(i);
                arma::mat current_beta_sigma_out_i = beta_sigma_out.submat(2 * i, 0, 2 * i + 1, 1);
                arma::rowvec diff_beta_i = current_beta_hat_out_i - alpha * current_beta_exp_i;
                diff_snp_scalar_vec(i) = arma::as_scalar(diff_beta_i * current_beta_sigma_out_i * diff_beta_i.t());
                r_snp_intermediate_diff_betas.row(i) = diff_beta_i;
            }
        }

        arma::mat r_snp_new = arma::zeros<arma::mat>(snp_number, 2);
        if (illegal_snps > 0 && snp_number > 0) {
            arma::uvec sorted_indices = arma::sort_index(diff_snp_scalar_vec, "descend");
            for (int j = 0; j < std::min(illegal_snps, snp_number); ++j) {
                r_snp_new.row(j) = r_snp_intermediate_diff_betas.row(sorted_indices(j));
            }
        }
        r_snp = r_snp_new;

        arma::mat beta_exp_new = beta_exp;
        if (snp_number > 0) {
            for (int i = 0; i < snp_number; ++i) {
                arma::colvec beta_hat_exp_i_col = beta_hat_exp.row(i).t();
                arma::rowvec beta_hat_out_i_row = beta_hat_out.row(i);
                arma::rowvec r_snp_i_row = r_snp.row(i);
                arma::colvec term_in_parentheses_col = (beta_hat_out_i_row - r_snp_i_row).t();
                arma::mat sigma_exp_i = beta_sigma_exp.submat(2 * i, 0, 2 * i + 1, 1);
                arma::mat sigma_out_i = beta_sigma_out.submat(2 * i, 0, 2 * i + 1, 1);
                arma::mat A_for_solve = sigma_exp_i + alpha * alpha * sigma_out_i;
                arma::colvec B_for_solve = sigma_exp_i * beta_hat_exp_i_col + alpha * (sigma_out_i * term_in_parentheses_col);
                arma::colvec beta_exp_new_i_col;
                bool solve_success = arma::solve(beta_exp_new_i_col, A_for_solve, B_for_solve);

                if (!solve_success) {
                    Rcpp::stop("arma::solve failed for beta_exp update at SNP index %d (0-based), iteration %d. Matrix A might be singular.", i + 1, number_of_iterations);
                }
                beta_exp_new.row(i) = beta_exp_new_i_col.t();
            }
        }
        beta_exp = beta_exp_new;

        double molecular_sum = 0.0;
        double denominator_sum = 0.0;

        if (snp_number > 0) {
            for (int i = 0; i < snp_number; ++i) {
                arma::rowvec beta_exp_i_row = beta_exp.row(i);
                arma::rowvec beta_hat_out_i_row = beta_hat_out.row(i);
                arma::mat sigma_out_i = beta_sigma_out.submat(2 * i, 0, 2 * i + 1, 1);
                arma::rowvec r_snp_i_row = r_snp.row(i);
                arma::rowvec term_for_molecular = beta_hat_out_i_row - r_snp_i_row;
                molecular_sum += arma::as_scalar(term_for_molecular * sigma_out_i * beta_exp_i_row.t());
                denominator_sum += arma::as_scalar(beta_exp_i_row * sigma_out_i * beta_exp_i_row.t());
            }
        }
        
        if (denominator_sum != 0.0) {
            alpha = molecular_sum / denominator_sum;
        } else {
            if (molecular_sum == 0.0) {
                 alpha = arma::datum::nan;
            } else {
                 alpha = (molecular_sum > 0 ? arma::datum::inf : -arma::datum::inf);
            }
        }
        
        termination_conditions_met_for_loop = should_continue_next_iteration;
    }

    bool convergence_flag;
    if (number_of_iterations >= maximum_number_of_iterations) {
        if (std::abs(alpha_last - alpha) > permissible_error) {
            convergence_flag = false;
        } else {
            convergence_flag = true;
        }
    } else {
        convergence_flag = true;
    }
    
    return Rcpp::List::create(
        Rcpp::Named("alpha") = alpha,
        Rcpp::Named("beta_exp") = beta_exp,
        Rcpp::Named("r_snp") = r_snp,
        Rcpp::Named("convergence") = convergence_flag,
        Rcpp::Named("iterations") = number_of_iterations
    );
}