#include <RcppArmadillo.h>
#include <vector>
#include <algorithm> // For std::sort, std::set_difference (though not directly used here for setdiff)
#include <cmath>     // For std::log, std::abs
#include <stdexcept> // For std::logic_error, std::runtime_error

// [[Rcpp::depends(RcppArmadillo)]]

//' Calculate BIC (Bayesian Information Criterion) - Rcpp version
//'
//' This function calculates BIC based on model outputs and input data.
//' It is an Rcpp implementation of the original R function.
//'
//' @param function_test_output List, containing 'complete_set', 'a_reduce_b_set',
//'                             'alpha_new', and 'beta_exp_updated'.
//' @param beta_hat_exp_r NumericMatrix, beta_hat estimates for exposure (n_snps x 2).
//' @param beta_hat_out_r NumericMatrix, beta_hat estimates for outcome (n_snps x 2).
//' @param matrix_big_r NumericMatrix, the large covariance matrix, expected to be (4*n_snps) x (4*n_snps).
//' @param min_n double, the minimum sample size.
//' @param a double, penalty parameter (default 0.0).
//' @param b double, penalty parameter (default 0.0).
//' @return double, the calculated BIC value, or NaN if calculation fails.
// [[Rcpp::export]]
double calculate_bic_rcpp(
    Rcpp::List function_test_output,
    Rcpp::NumericMatrix beta_hat_exp_r,
    Rcpp::NumericMatrix beta_hat_out_r,
    Rcpp::NumericMatrix matrix_big_r,
    double min_n,
    double a = 0.0,
    double b = 0.0) {

    // --- 1. Parse inputs from R objects to Armadillo objects ---
    Rcpp::IntegerVector complete_set_r_vec;
    Rcpp::IntegerVector a_reduce_b_set_r_vec;
    double alpha_val;
    Rcpp::NumericMatrix beta_exp_r_mat;

    try {
        complete_set_r_vec = Rcpp::as<Rcpp::IntegerVector>(function_test_output["complete_set"]);
        a_reduce_b_set_r_vec = Rcpp::as<Rcpp::IntegerVector>(function_test_output["a_reduce_b_set"]);
        alpha_val = Rcpp::as<double>(function_test_output["alpha_new"]);
        beta_exp_r_mat = Rcpp::as<Rcpp::NumericMatrix>(function_test_output["beta_exp_updated"]);
    } catch (Rcpp::not_compatible& e) {
        Rcpp::stop("Error accessing elements from 'function_test_output' list: " + std::string(e.what()));
    } catch (...) {
        Rcpp::stop("Unknown error accessing elements from 'function_test_output' list.");
    }

    arma::mat beta_exp_arma = Rcpp::as<arma::mat>(beta_exp_r_mat);
    arma::mat beta_hat_exp_arma = Rcpp::as<arma::mat>(beta_hat_exp_r);
    arma::mat beta_hat_out_arma = Rcpp::as<arma::mat>(beta_hat_out_r);
    arma::mat matrix_big_arma = Rcpp::as<arma::mat>(matrix_big_r);

    int n_snps = beta_hat_exp_arma.n_rows;
    if (n_snps == 0) {
        Rcpp::warning("Number of SNPs (n_snps) is 0. BIC calculation might be trivial or invalid.");
    }
    if (beta_exp_arma.n_rows != (arma::uword)n_snps || beta_hat_out_arma.n_rows != (arma::uword)n_snps) {
        Rcpp::stop("Row count mismatch between beta_hat_exp, beta_hat_out, and beta_exp_updated.");
    }
    if (matrix_big_arma.n_rows != (arma::uword)(4*n_snps) || (n_snps > 0 && matrix_big_arma.n_cols != 4)) {
        if (n_snps > 0) { // Only warn if n_snps > 0
             Rcpp::warning("matrix_big_r dimensions are not (4*n_snps) x 4. Expected " +
                        std::to_string(4*n_snps) + "x4" +
                        ", got " + std::to_string(matrix_big_arma.n_rows) + "x" + std::to_string(matrix_big_arma.n_cols));
        }
    }


    // Convert 1-based R indices from IntegerVector to 0-based arma::uvec
    arma::uvec complete_set_arma_0based(complete_set_r_vec.size());
    for (int i = 0; i < complete_set_r_vec.size(); ++i) {
        if (complete_set_r_vec[i] <= 0 || complete_set_r_vec[i] > n_snps) {
            Rcpp::stop("Invalid index in complete_set: " + std::to_string(complete_set_r_vec[i]));
        }
        complete_set_arma_0based(i) = complete_set_r_vec[i] - 1;
    }

    arma::uvec a_reduce_b_set_arma_0based(a_reduce_b_set_r_vec.size());
    for (int i = 0; i < a_reduce_b_set_r_vec.size(); ++i) {
         if (a_reduce_b_set_r_vec[i] <= 0 || a_reduce_b_set_r_vec[i] > n_snps) {
            Rcpp::stop("Invalid index in a_reduce_b_set: " + std::to_string(a_reduce_b_set_r_vec[i]));
        }
        a_reduce_b_set_arma_0based(i) = a_reduce_b_set_r_vec[i] - 1;
    }

    double l_complete = 0.0;
    double l_reduce = 0.0;
    double l_cc = 0.0;

    // --- 2.1 Calculate likelihood contribution from 'complete_set' ---
    if (complete_set_arma_0based.n_elem > 0) {
        for (arma::uword k = 0; k < complete_set_arma_0based.n_elem; ++k) {
            arma::uword cpp_idx = complete_set_arma_0based(k);
            
            arma::uword block_start = 4 * cpp_idx;
            arma::uword block_end = block_start + 3;
            
            arma::mat matrix_big_i_4x4; // The 4x4 block for SNP i
            try {
                 matrix_big_i_4x4 = matrix_big_arma.submat(block_start, 0, block_end, 3); // Select all 4 columns
            } catch (const std::logic_error& e) {
                Rcpp::stop("Submatrix extraction failed for matrix_big_i_4x4 (complete_set) at cpp_idx " + 
                           std::to_string(cpp_idx) + ". Error: " + e.what());
            }

            arma::mat matrix_big_i_inv;
            try {
                matrix_big_i_inv = arma::inv_sympd(matrix_big_i_4x4);
            } catch (const std::runtime_error& e) {
                Rcpp::warning("Matrix inversion failed for matrix_big_i_4x4 (complete_set, cpp_idx=" + 
                              std::to_string(cpp_idx) + "). Error: " + std::string(e.what()));
                return R_NaN;
            }

            arma::rowvec beta_hat_exp_i_row = beta_hat_exp_arma.row(cpp_idx); // Expected 1x2
            arma::rowvec beta_hat_out_i_row = beta_hat_out_arma.row(cpp_idx); // Expected 1x2
            arma::rowvec beta_exp_i_row = beta_exp_arma.row(cpp_idx);         // Expected 1x2

            arma::vec d_diff(4);
            d_diff.subvec(0,1) = (beta_hat_exp_i_row - beta_exp_i_row).t(); 
            d_diff.subvec(2,3) = (beta_hat_out_i_row - alpha_val * beta_exp_i_row).t();
            
            l_complete += arma::as_scalar(d_diff.t() * matrix_big_i_inv * d_diff);
        }
    }

    // --- 2.2 Calculate likelihood contribution from 'a_reduce_b_set' ---
    if (a_reduce_b_set_arma_0based.n_elem > 0) {
        for (arma::uword k = 0; k < a_reduce_b_set_arma_0based.n_elem; ++k) {
            arma::uword cpp_idx = a_reduce_b_set_arma_0based(k);

            arma::uword block_start = 4 * cpp_idx;
            arma::uword block_end = block_start + 3;

            arma::mat matrix_big_i_sub_4x4; // 4x4 block for SNP i
            try {
                matrix_big_i_sub_4x4 = matrix_big_arma.submat(block_start, 0, block_end, 3); // Select all 4 columns
            } catch (const std::logic_error& e) {
                Rcpp::stop("Submatrix extraction failed for matrix_big_i_sub_4x4 (a_reduce_b_set) at cpp_idx " + 
                           std::to_string(cpp_idx) + ". Error: " + e.what());
            }

            arma::mat matrix_o_i(2, 2);
            matrix_o_i(0, 0) = matrix_big_i_sub_4x4(0, 0); // R: matrix_big_i_sub[1,1]
            matrix_o_i(1, 1) = matrix_big_i_sub_4x4(2, 2); // R: matrix_big_i_sub[3,3]
            matrix_o_i(0, 1) = matrix_big_i_sub_4x4(0, 2); // R: matrix_big_i_sub[1,3]
            matrix_o_i(1, 0) = matrix_big_i_sub_4x4(0, 2); // R: matrix_big_i_sub[1,3] for [2,1]

            arma::mat matrix_o_i_inv;
            try {
                matrix_o_i_inv = arma::inv_sympd(matrix_o_i);
            } catch (const std::runtime_error& e) {
                Rcpp::warning("Matrix inversion failed for matrix_o_i (a_reduce_b_set, cpp_idx=" + 
                              std::to_string(cpp_idx) + "). Error: " + std::string(e.what()));
                return R_NaN;
            }
            
            arma::vec d_o(2);
            d_o(0) = beta_hat_exp_arma(cpp_idx, 0) - beta_exp_arma(cpp_idx, 0);
            d_o(1) = beta_hat_out_arma(cpp_idx, 0) - alpha_val * beta_exp_arma(cpp_idx, 0);
            
            l_reduce += arma::as_scalar(d_o.t() * matrix_o_i_inv * d_o);
            
            double term2_num_sq = beta_hat_exp_arma(cpp_idx, 1) - beta_exp_arma(cpp_idx, 1);
            term2_num_sq *= term2_num_sq; // Square it
            double term2_den = matrix_big_i_sub_4x4(1, 1); // R: matrix_big_i_sub[2,2] -> C++ (1,1) of 4x4 block

            if (std::abs(term2_den) < 1e-12) { 
                 Rcpp::warning("Denominator matrix_big_i_sub_4x4(1,1) is near zero for a_reduce_b_set, cpp_idx=" + 
                               std::to_string(cpp_idx) + ". Adding Inf to l_reduce.");
                 l_reduce += R_PosInf; 
            } else {
                l_reduce += term2_num_sq / term2_den;
            }
        }
    }
    
    // --- 2.3 Calculate likelihood contribution from remaining set (a_c_b_c_set) ---
    std::vector<bool> is_used_snp(n_snps, false);
    for(arma::uword idx : complete_set_arma_0based) { if(idx < (arma::uword)n_snps) is_used_snp[idx] = true; }
    for(arma::uword idx : a_reduce_b_set_arma_0based) { if(idx < (arma::uword)n_snps) is_used_snp[idx] = true; }

    arma::uvec a_c_b_c_set_arma_0based; 
    for (int j = 0; j < n_snps; ++j) {
        if (!is_used_snp[j]) {
            a_c_b_c_set_arma_0based.insert_rows(a_c_b_c_set_arma_0based.n_elem, arma::uvec{ (arma::uword)j });
        }
    }

    if (a_c_b_c_set_arma_0based.n_elem > 0) {
        for (arma::uword k = 0; k < a_c_b_c_set_arma_0based.n_elem; ++k) {
            arma::uword cpp_idx = a_c_b_c_set_arma_0based(k);

            // R: matrix_index_exp <- (4*i-3):(4*i-2) -> C++ 0-based: rows 4*cpp_idx and 4*cpp_idx+1
            arma::uword exp_block_start = 4 * cpp_idx;
            arma::uword exp_block_end = exp_block_start + 1; 

            arma::mat matrix_exp_i_2x2; // The 2x2 block
            try {
                matrix_exp_i_2x2 = matrix_big_arma.submat(exp_block_start, 0, exp_block_end, 1); // Select first 2 columns
            } catch (const std::logic_error& e) {
                 Rcpp::stop("Submatrix extraction failed for matrix_exp_i_2x2 (a_c_b_c_set) at cpp_idx " + 
                            std::to_string(cpp_idx) + ". Error: " + e.what());
            }

            arma::mat matrix_exp_i_inv;
            try {
                matrix_exp_i_inv = arma::inv_sympd(matrix_exp_i_2x2);
            } catch (const std::runtime_error& e) {
                Rcpp::warning("Matrix inversion failed for matrix_exp_i_2x2 (a_c_b_c_set, cpp_idx=" + 
                              std::to_string(cpp_idx) + "). Error: " + std::string(e.what()));
                return R_NaN;
            }

            arma::rowvec beta_hat_exp_i_row = beta_hat_exp_arma.row(cpp_idx); // Expected 1x2
            arma::rowvec beta_exp_i_row = beta_exp_arma.row(cpp_idx);         // Expected 1x2
            
            arma::vec diff_exp_2x1 = (beta_hat_exp_i_row - beta_exp_i_row).t(); 
            
            l_cc += arma::as_scalar(diff_exp_2x1.t() * matrix_exp_i_inv * diff_exp_2x1);
        }
    }

    // --- 3. Calculate final BIC value ---
    if (min_n <= 0) {
        Rcpp::warning("min_n is non-positive, log(min_n) will be NaN or -Inf.");
    }
    double log_min_n = std::log(min_n);
    if (R_IsNaN(log_min_n) || !R_finite(log_min_n)) {
         Rcpp::warning("log(min_n) resulted in NaN or Inf/ -Inf. BIC will be affected.");
    }
    
    double bic = l_complete + l_reduce + l_cc + log_min_n * (a + b);
    
    if (R_IsNaN(bic) || !R_finite(bic)) {
        Rcpp::warning("Final BIC value is NaN or Inf. Check intermediate calculations and inputs.");
    }

    return bic;
}