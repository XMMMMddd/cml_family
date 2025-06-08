#include <RcppArmadillo.h>
#include <algorithm> // For std::sort, std::set_intersection
#include <vector>    // For std::vector

// [[Rcpp::depends(RcppArmadillo)]]

// Helper function to find indices of all-zero rows in a matrix
arma::uvec find_zero_rows(const arma::mat& M) {
    arma::uvec zero_rows_indices;
    if (M.n_rows == 0) {
        return zero_rows_indices;
    }
    for (arma::uword i = 0; i < M.n_rows; ++i) {
        if (arma::all(M.row(i) == 0.0)) {
            zero_rows_indices.insert_rows(zero_rows_indices.n_elem, arma::uvec{i});
        }
    }
    return zero_rows_indices;
}

// Helper function to find indices of all-zero columns in a matrix
arma::uvec find_zero_cols(const arma::mat& M) {
    arma::uvec zero_cols_indices;
    if (M.n_cols == 0) {
        return zero_cols_indices;
    }
    for (arma::uword i = 0; i < M.n_cols; ++i) {
        if (arma::all(M.col(i) == 0.0)) {
            zero_cols_indices.insert_rows(zero_cols_indices.n_elem, arma::uvec{i});
        }
    }
    return zero_cols_indices;
}

// Helper function to find the intersection of two unsigned integer vectors (arma::uvec)
arma::uvec intersect_uvec(const arma::uvec& v1, const arma::uvec& v2) {
    std::vector<arma::uword> vec1_std(v1.begin(), v1.end());
    std::vector<arma::uword> vec2_std(v2.begin(), v2.end());
    std::sort(vec1_std.begin(), vec1_std.end());
    std::sort(vec2_std.begin(), vec2_std.end());
    std::vector<arma::uword> common_elements_std;
    std::set_intersection(vec1_std.begin(), vec1_std.end(),
                          vec2_std.begin(), vec2_std.end(),
                          std::back_inserter(common_elements_std));
    return arma::conv_to<arma::uvec>::from(common_elements_std);
}


//' Calculate Alpha Standard Error (Rcpp version)
//'
//' This function calculates the standard error of alpha using matrix operations.
//' It is an Rcpp implementation of the original R function.
//'
//' @param complete_set IntegerVector, 1-based indices for the complete set.
//' @param a_reduce_b_set IntegerVector, 1-based indices for the a_reduce_b set.
//' @param matrix_big_r NumericMatrix, the large input matrix.
//' @param beta_hat_exp_r NumericMatrix, beta_hat estimates for exposure.
//' @param beta_hat_out_r NumericMatrix, beta_hat estimates for outcome.
//' @param beta_exp_r NumericMatrix, true beta for exposure.
//' @param alpha double, the alpha parameter.
//' @return double, the calculated standard error of alpha, or NaN if calculation fails.
// [[Rcpp::export]]
double calculate_alpha_se_rcpp(
    Rcpp::IntegerVector complete_set,
    Rcpp::IntegerVector a_reduce_b_set,
    Rcpp::NumericMatrix matrix_big_r,
    Rcpp::NumericMatrix beta_hat_exp_r,
    Rcpp::NumericMatrix beta_hat_out_r,
    Rcpp::NumericMatrix beta_exp_r,
    double alpha) {

    // Convert Rcpp objects to Armadillo objects
    arma::ivec complete_set_arma = Rcpp::as<arma::ivec>(complete_set);
    arma::ivec a_reduce_b_set_arma = Rcpp::as<arma::ivec>(a_reduce_b_set);
    arma::mat matrix_big_arma = Rcpp::as<arma::mat>(matrix_big_r);
    arma::mat beta_hat_exp_arma = Rcpp::as<arma::mat>(beta_hat_exp_r);
    arma::mat beta_hat_out_arma = Rcpp::as<arma::mat>(beta_hat_out_r);
    arma::mat beta_exp_arma = Rcpp::as<arma::mat>(beta_exp_r);

    // --- Input Validation (Recommended) ---
    if (complete_set_arma.n_elem > 0) {
        if (complete_set_arma.max() > matrix_big_arma.n_rows / 4 || 
            complete_set_arma.max() > beta_hat_exp_arma.n_rows) {
            Rcpp::stop("An index in 'complete_set' is out of bounds for the provided data matrices.");
        }
    }
    if (a_reduce_b_set_arma.n_elem > 0) {
        if (a_reduce_b_set_arma.max() > matrix_big_arma.n_rows / 4 ||
            a_reduce_b_set_arma.max() > beta_hat_exp_arma.n_rows) {
            Rcpp::stop("An index in 'a_reduce_b_set' is out of bounds for the provided data matrices.");
        }
    }

    // Determine the dimension of the gamma part of the H matrix
    int h_dimension = complete_set_arma.n_elem * 2 + a_reduce_b_set_arma.n_elem * 2;
    
    // Initialize components of the H matrix
    arma::mat h_matrix_arma(h_dimension + 1, h_dimension + 1, arma::fill::zeros);
    arma::mat h_a_a_arma(1, 1, arma::fill::zeros);
    arma::mat h_a_gamma_arma(h_dimension, 1, arma::fill::zeros);
    arma::mat h_gamma_gamma_arma(h_dimension, h_dimension, arma::fill::zeros);

    // --- Loop over the 'complete_set' ---
    for (arma::uword k = 0; k < complete_set_arma.n_elem; ++k) {
        int r_index = complete_set_arma(k);
        int cpp_index = r_index - 1;

        arma::uword row_start_big = 4 * cpp_index;
        arma::mat matrix_big_i = matrix_big_arma.rows(row_start_big, row_start_big + 3);
        
        arma::mat matrix_big_i_inv;
        try {
            matrix_big_i_inv = arma::inv_sympd(matrix_big_i);
        } catch (const std::runtime_error& e) {
            Rcpp::warning("Matrix inversion failed (complete_set, r_index=" + std::to_string(r_index) + "). Error: " + std::string(e.what()));
            return R_NaN;
        }

        arma::mat omega_11 = matrix_big_i_inv.submat(0, 0, 1, 1);
        arma::mat omega_21 = matrix_big_i_inv.submat(2, 0, 3, 1);
        arma::mat omega_12 = matrix_big_i_inv.submat(0, 2, 1, 3);
        arma::mat omega_22 = matrix_big_i_inv.submat(2, 2, 3, 3);

        arma::mat beta_hat_exp_i_colvec = beta_hat_exp_arma.row(cpp_index).t();
        arma::mat beta_hat_out_i_colvec = beta_hat_out_arma.row(cpp_index).t();
        arma::mat beta_exp_i_colvec = beta_exp_arma.row(cpp_index).t();

        h_a_a_arma += beta_hat_exp_i_colvec.t() * omega_22 * beta_hat_exp_i_colvec;
        
        // =================================================================
        // *** KEY FIX #1: Use loop counter 'k' for H-matrix indexing ***
        arma::uword row_start_gamma_block = 2 * k; 
        // =================================================================
        arma::uword row_end_gamma_block = row_start_gamma_block + 1;

        h_a_gamma_arma.submat(row_start_gamma_block, 0, row_end_gamma_block, 0) = 
            -(omega_12 * beta_hat_exp_i_colvec) +
            (omega_12 + omega_21) * beta_exp_i_colvec - 
            omega_22 * beta_hat_out_i_colvec +
            2 * alpha * omega_22 * beta_hat_exp_i_colvec;

        h_gamma_gamma_arma.submat(row_start_gamma_block, row_start_gamma_block, 
                                  row_end_gamma_block, row_end_gamma_block) = 
            omega_11 + alpha * (omega_12 + omega_21) +
            alpha * alpha * omega_22;
    }

    // --- Loop over the 'a_reduce_b_set' ---
    for (arma::uword k = 0; k < a_reduce_b_set_arma.n_elem; ++k) {
        int r_index = a_reduce_b_set_arma(k);
        int cpp_index = r_index - 1;

        arma::mat matrix_big_i_subset = matrix_big_arma.rows(4 * cpp_index, 4 * cpp_index + 3);

        arma::mat matrix_o_i(2, 2);
        matrix_o_i(0, 0) = matrix_big_i_subset(0, 0);
        matrix_o_i(1, 1) = matrix_big_i_subset(2, 2);
        matrix_o_i(0, 1) = matrix_big_i_subset(0, 2);
        matrix_o_i(1, 0) = matrix_big_i_subset(0, 2);

        arma::mat matrix_o_i_inv;
        try {
            matrix_o_i_inv = arma::inv_sympd(matrix_o_i);
        } catch (const std::runtime_error& e) {
            Rcpp::warning("Matrix inversion failed (a_reduce_b_set, r_index=" + std::to_string(r_index) + "). Error: " + std::string(e.what()));
            return R_NaN;
        }

        double beta_hat_exp_i_val = beta_hat_exp_arma(cpp_index, 0);
        double beta_hat_out_i_val = beta_hat_out_arma(cpp_index, 0);
        double beta_exp_i_val = beta_exp_arma(cpp_index, 0);

        h_a_a_arma(0,0) += beta_hat_exp_i_val * beta_hat_exp_i_val * matrix_o_i_inv(1, 1);

        // =========================================================================================
        // *** KEY FIX #2: Use loop counter 'k' plus offset for H-matrix indexing ***
        arma::uword row_start_gamma_block = complete_set_arma.n_elem * 2 + 2 * k;
        // =========================================================================================
        
        h_a_gamma_arma(row_start_gamma_block, 0) = 
            -(matrix_o_i_inv(1, 0) * beta_hat_exp_i_val -
              matrix_o_i_inv(1, 1) * beta_hat_out_i_val +
              2 * beta_exp_i_val * (matrix_o_i_inv(1, 0) + alpha * matrix_o_i_inv(1, 1)));
        h_a_gamma_arma(row_start_gamma_block + 1, 0) = 0.0;

        double h_gamma_gamma_input_val = matrix_o_i_inv(0, 0) +
                                       alpha * (matrix_o_i_inv(0, 1) + matrix_o_i_inv(1, 0)) +
                                       alpha * alpha * matrix_o_i_inv(1, 1);
        
        arma::mat h_gamma_gamma_input_matrix(2, 2, arma::fill::zeros);
        h_gamma_gamma_input_matrix(0, 0) = h_gamma_gamma_input_val;

        h_gamma_gamma_arma.submat(row_start_gamma_block, row_start_gamma_block, 
                                  row_start_gamma_block + 1, row_start_gamma_block + 1) = h_gamma_gamma_input_matrix;
    }

    // --- Assemble the full H matrix ---
    if (h_dimension == 0) {
        if (h_a_a_arma.n_elem > 0) h_matrix_arma(0, 0) = h_a_a_arma(0,0);
    } else {
        h_matrix_arma(0, 0) = h_a_a_arma(0,0);
        h_matrix_arma.submat(1, 0, h_dimension, 0) = h_a_gamma_arma; 
        h_matrix_arma.submat(0, 1, 0, h_dimension) = h_a_gamma_arma.t();
        h_matrix_arma.submat(1, 1, h_dimension, h_dimension) = h_gamma_gamma_arma;
    }

    // --- Process h_matrix and calculate SE ---
    double alpha_se_value = R_NaN;

    if (h_matrix_arma.n_rows == 0) {
        Rcpp::warning("h_matrix is empty, cannot process.");
        return R_NaN;
    }
    
    arma::uvec zero_row_indices = find_zero_rows(h_matrix_arma);
    arma::uvec zero_col_indices = find_zero_cols(h_matrix_arma);
    arma::uvec common_zero_indices = intersect_uvec(zero_row_indices, zero_col_indices);

    arma::mat h_matrix_final = h_matrix_arma;
    if (common_zero_indices.n_elem > 0) {
        h_matrix_final.shed_cols(common_zero_indices); 
        h_matrix_final.shed_rows(common_zero_indices); 
    }

    if (h_matrix_final.n_rows > 0 && h_matrix_final.is_square()) {
        if (std::abs(arma::det(h_matrix_final)) > 1e-12) {
            try {
                arma::mat inv_h_matrix = arma::inv_sympd(h_matrix_final);
                alpha_se_value = inv_h_matrix(0, 0);
            } catch (const std::runtime_error& e) {
                 Rcpp::warning("Inversion of the final h_matrix failed. Error: " + std::string(e.what()));
            }
        } else {
            Rcpp::warning("Final h_matrix is singular (det ~ 0), cannot calculate alpha_se.");
        }
    } else {
        Rcpp::warning("Final h_matrix is empty or not square after cleaning.");
    }
    
    if (R_IsNaN(alpha_se_value) || alpha_se_value < 0) {
        if (!R_IsNaN(alpha_se_value) && alpha_se_value < 0) {
            Rcpp::warning("alpha_se_value is negative before sqrt, result will be NaN.");
        }
        return R_NaN;
    }

    return std::sqrt(alpha_se_value);
}