#include <RcppArmadillo.h>
#include <algorithm> // For std::sort, std::set_intersection
#include <vector>    // For std::vector

// [[Rcpp::depends(RcppArmadillo)]]

// Helper function to find indices of all-zero rows in a matrix
// Input: M - an Armadillo matrix
// Output: arma::uvec containing 0-based indices of rows where all elements are zero
arma::uvec find_zero_rows(const arma::mat& M) {
    arma::uvec zero_rows_indices;
    if (M.n_rows == 0) {
        return zero_rows_indices; // Return empty if matrix has no rows
    }
    for (arma::uword i = 0; i < M.n_rows; ++i) {
        if (arma::all(M.row(i) == 0.0)) { // Check if all elements in the row are zero
            // Add index i to the list of zero rows
            zero_rows_indices.insert_rows(zero_rows_indices.n_elem, arma::uvec{i});
        }
    }
    return zero_rows_indices;
}

// Helper function to find indices of all-zero columns in a matrix
// Input: M - an Armadillo matrix
// Output: arma::uvec containing 0-based indices of columns where all elements are zero
arma::uvec find_zero_cols(const arma::mat& M) {
    arma::uvec zero_cols_indices;
    if (M.n_cols == 0) {
        return zero_cols_indices; // Return empty if matrix has no columns
    }
    for (arma::uword i = 0; i < M.n_cols; ++i) {
        if (arma::all(M.col(i) == 0.0)) { // Check if all elements in the column are zero
            // Add index i to the list of zero columns
            zero_cols_indices.insert_rows(zero_cols_indices.n_elem, arma::uvec{i});
        }
    }
    return zero_cols_indices;
}

// Helper function to find the intersection of two unsigned integer vectors (arma::uvec)
// Inputs: v1, v2 - two Armadillo uvecs
// Output: arma::uvec containing common elements, sorted
arma::uvec intersect_uvec(const arma::uvec& v1, const arma::uvec& v2) {
    // Convert arma::uvec to std::vector for using STL algorithms
    std::vector<arma::uword> vec1_std(v1.begin(), v1.end());
    std::vector<arma::uword> vec2_std(v2.begin(), v2.end());

    // Sort vectors before finding intersection
    std::sort(vec1_std.begin(), vec1_std.end());
    std::sort(vec2_std.begin(), vec2_std.end());

    std::vector<arma::uword> common_elements_std;
    // Find common elements and store them in common_elements_std
    std::set_intersection(vec1_std.begin(), vec1_std.end(),
                          vec2_std.begin(), vec2_std.end(),
                          std::back_inserter(common_elements_std));
    
    // Convert the result back to arma::uvec
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

    // Convert Rcpp objects to Armadillo objects for efficient matrix operations
    arma::ivec complete_set_arma = Rcpp::as<arma::ivec>(complete_set);
    arma::ivec a_reduce_b_set_arma = Rcpp::as<arma::ivec>(a_reduce_b_set);
    arma::mat matrix_big_arma = Rcpp::as<arma::mat>(matrix_big_r);
    arma::mat beta_hat_exp_arma = Rcpp::as<arma::mat>(beta_hat_exp_r);
    arma::mat beta_hat_out_arma = Rcpp::as<arma::mat>(beta_hat_out_r);
    arma::mat beta_exp_arma = Rcpp::as<arma::mat>(beta_exp_r);

    // Determine the dimension of the gamma part of the H matrix
    int h_dimension = complete_set_arma.n_elem * 2 + a_reduce_b_set_arma.n_elem * 2;
    
    // Initialize components of the H matrix
    // h_matrix_arma will be (h_dimension + 1) x (h_dimension + 1)
    // h_a_a_arma is a scalar (1x1 matrix) representing the (1,1) element related to alpha
    // h_a_gamma_arma is a (h_dimension x 1) column vector for interactions between alpha and gamma
    // h_gamma_gamma_arma is a (h_dimension x h_dimension) matrix for gamma-gamma interactions
    arma::mat h_matrix_arma(h_dimension + 1, h_dimension + 1, arma::fill::zeros);
    arma::mat h_a_a_arma(1, 1, arma::fill::zeros); // Scalar component
    arma::mat h_a_gamma_arma(h_dimension, 1, arma::fill::zeros); // Alpha-gamma interaction part
    arma::mat h_gamma_gamma_arma(h_dimension, h_dimension, arma::fill::zeros); // Gamma-gamma interaction part

    // --- Loop over the 'complete_set' ---
    // Note: Indices from R (complete_set_arma, a_reduce_b_set_arma) are 1-based.
    // Armadillo/C++ uses 0-based indexing. So, conversion is needed.
    for (arma::uword k = 0; k < complete_set_arma.n_elem; ++k) {
        int r_index = complete_set_arma(k); // Get 1-based R index
        int cpp_index = r_index - 1;      // Convert to 0-based C++ index

        // In R: matrix_indicator <- (4 * r_index - 3):(4 * r_index)
        // This corresponds to 0-based rows: (4*r_index - 4) to (4*r_index - 1)
        // Or, using cpp_index: (4*cpp_index) to (4*cpp_index + 3)
        arma::uword row_start_big = 4 * cpp_index;
        arma::uword row_end_big = row_start_big + 3;
        
        // Extract the relevant 4xN submatrix from matrix_big_arma
        arma::mat matrix_big_i = matrix_big_arma.rows(row_start_big, row_end_big);
        arma::mat matrix_big_i_inv;
        try {
            // Attempt to invert matrix_big_i
            matrix_big_i_inv = arma::inv_sympd(matrix_big_i); // Assuming symmetric positive definite, use inv() if not
        } catch (const std::runtime_error& e) {
            Rcpp::warning("Matrix inversion failed for matrix_big_i (complete_set loop, r_index=" + std::to_string(r_index) + "). Error: " + std::string(e.what()));
            return R_NaN; // Return NaN on failure
        }

        // Extract omega submatrices from the inverted matrix_big_i_inv
        // These are 2x2 blocks
        arma::mat omega_11 = matrix_big_i_inv.submat(0, 0, 1, 1);
        arma::mat omega_21 = matrix_big_i_inv.submat(2, 0, 3, 1);
        arma::mat omega_12 = matrix_big_i_inv.submat(0, 2, 1, 3);
        arma::mat omega_22 = matrix_big_i_inv.submat(2, 2, 3, 3);

        // Extract beta values for the current index i
        // These are expected to be 1x2 row vectors from the input matrices
        arma::rowvec beta_hat_exp_i_rowvec = beta_hat_exp_arma.row(cpp_index);
        arma::rowvec beta_hat_out_i_rowvec = beta_hat_out_arma.row(cpp_index);
        arma::rowvec beta_exp_i_rowvec = beta_exp_arma.row(cpp_index);
        
        // Transpose to 2x1 column vectors for matrix multiplication
        arma::mat beta_hat_exp_i_colvec = beta_hat_exp_i_rowvec.t();
        arma::mat beta_hat_out_i_colvec = beta_hat_out_i_rowvec.t();
        arma::mat beta_exp_i_colvec = beta_exp_i_rowvec.t();

        // Update h_a_a component
        h_a_a_arma += beta_hat_exp_i_colvec.t() * omega_22 * beta_hat_exp_i_colvec;
        
        // Define row indices for h_a_gamma_arma and h_gamma_gamma_arma
        // In R: (2 * r_index - 1):(2 * r_index)
        // 0-based: (2 * cpp_index) to (2 * cpp_index + 1)
        arma::uword row_start_gamma_block = 2 * cpp_index;
        arma::uword row_end_gamma_block = row_start_gamma_block + 1;

        // Update h_a_gamma_arma component (a 2x1 block)
        h_a_gamma_arma.submat(row_start_gamma_block, 0, row_end_gamma_block, 0) = 
            -(omega_12 * beta_hat_exp_i_colvec) +
            (omega_12 + omega_21) * beta_exp_i_colvec - 
            omega_22 * beta_hat_out_i_colvec +
            2 * alpha * omega_22 * beta_hat_exp_i_colvec;

        // Update h_gamma_gamma_arma component (a 2x2 diagonal block)
        h_gamma_gamma_arma.submat(row_start_gamma_block, row_start_gamma_block, 
                                  row_end_gamma_block, row_end_gamma_block) = 
            omega_11 + alpha * (omega_12 + omega_21) +
            alpha * alpha * omega_22;
    }

    // --- Loop over the 'a_reduce_b_set' ---
    for (arma::uword k = 0; k < a_reduce_b_set_arma.n_elem; ++k) {
        int r_index = a_reduce_b_set_arma(k); // Get 1-based R index
        int cpp_index = r_index - 1;        // Convert to 0-based C++ index

        // Extract the 4 rows for current r_index from matrix_big_arma
        arma::uword row_start_big = 4 * cpp_index;
        arma::mat matrix_big_i_subset = matrix_big_arma.rows(row_start_big, row_start_big + 3);

        // Construct matrix_o_i (2x2) as per R logic
        arma::mat matrix_o_i(2, 2);
        matrix_o_i(0, 0) = matrix_big_i_subset(0, 0); // R: matrix_big_i[1,1]
        matrix_o_i(1, 1) = matrix_big_i_subset(2, 2); // R: matrix_big_i[3,3]
        matrix_o_i(0, 1) = matrix_big_i_subset(0, 2); // R: matrix_big_i[1,3]
        matrix_o_i(1, 0) = matrix_big_i_subset(0, 2); // R: matrix_big_i[1,3] (used for [2,1] as well)

        arma::mat matrix_o_i_inv;
        try {
            // Attempt to invert matrix_o_i
            matrix_o_i_inv = arma::inv_sympd(matrix_o_i); // Assuming symmetric positive definite
        } catch (const std::runtime_error& e) {
            Rcpp::warning("Matrix inversion failed for matrix_o_i (a_reduce_b_set loop, r_index=" + std::to_string(r_index) + "). Error: " + std::string(e.what()));
            return R_NaN; // Return NaN on failure
        }

        // Extract scalar beta values (first element of the row)
        double beta_hat_exp_i_val = beta_hat_exp_arma(cpp_index, 0); // R: beta_hat_exp[i,1]
        double beta_hat_out_i_val = beta_hat_out_arma(cpp_index, 0); // R: beta_hat_out[i,1]
        double beta_exp_i_val = beta_exp_arma(cpp_index, 0);       // R: beta_exp[i,1]

        // Update h_a_a_arma component
        h_a_a_arma(0,0) += beta_hat_exp_i_val * beta_hat_exp_i_val * matrix_o_i_inv(1, 1);

        // Define row indices for h_a_gamma_arma and h_gamma_gamma_arma
        arma::uword row_start_gamma_block = 2 * cpp_index;

        // Update h_a_gamma_arma component (a 2x1 block, second element is 0)
        h_a_gamma_arma(row_start_gamma_block, 0) = 
            -(matrix_o_i_inv(1, 0) * beta_hat_exp_i_val -
              matrix_o_i_inv(1, 1) * beta_hat_out_i_val +
              2 * beta_exp_i_val * (matrix_o_i_inv(1, 0) + alpha * matrix_o_i_inv(1, 1)));
        h_a_gamma_arma(row_start_gamma_block + 1, 0) = 0.0; // Second element is explicitly zero

        // Calculate the scalar input for the h_gamma_gamma_arma block
        double h_gamma_gamma_input_val = matrix_o_i_inv(0, 0) +
                                       alpha * (matrix_o_i_inv(0, 1) + matrix_o_i_inv(1, 0)) +
                                       alpha * alpha * matrix_o_i_inv(1, 1);
        
        // Construct the 2x2 h_gamma_gamma_input_matrix (only (0,0) element is non-zero)
        arma::mat h_gamma_gamma_input_matrix(2, 2, arma::fill::zeros);
        h_gamma_gamma_input_matrix(0, 0) = h_gamma_gamma_input_val;

        // Update h_gamma_gamma_arma component (a 2x2 diagonal block)
        h_gamma_gamma_arma.submat(row_start_gamma_block, row_start_gamma_block, 
                                  row_start_gamma_block + 1, row_start_gamma_block + 1) = h_gamma_gamma_input_matrix;
    }

    // --- Assemble the full H matrix (h_matrix_arma) ---
    if (h_dimension == 0) { // If no gamma components, H is just h_a_a
        h_matrix_arma(0, 0) = h_a_a_arma(0,0);
    } else { // If gamma components exist
        h_matrix_arma(0, 0) = h_a_a_arma(0,0); // Top-left element
        // Off-diagonal block (alpha-gamma interaction)
        h_matrix_arma.submat(1, 0, h_dimension, 0) = h_a_gamma_arma; 
        // Off-diagonal block (gamma-alpha interaction, transpose of above)
        h_matrix_arma.submat(0, 1, 0, h_dimension) = h_a_gamma_arma.t();
        // Bottom-right block (gamma-gamma interaction)
        h_matrix_arma.submat(1, 1, h_dimension, h_dimension) = h_gamma_gamma_arma;
    }

    // --- Process h_matrix_arma: check for all-zero rows/columns and calculate SE ---
    double alpha_se_value = R_NaN; // Initialize result to NaN (R's NA)

    if (h_matrix_arma.n_rows == 0 || h_matrix_arma.n_cols == 0) {
        Rcpp::warning("h_matrix is empty, cannot process.");
        return R_NaN;
    }
    
    // Find indices of rows and columns that are entirely zero
    arma::uvec zero_row_indices = find_zero_rows(h_matrix_arma);
    arma::uvec zero_col_indices = find_zero_cols(h_matrix_arma);
    
    // Find indices that are common to both all-zero rows and all-zero columns
    arma::uvec common_zero_indices = intersect_uvec(zero_row_indices, zero_col_indices);

    if (common_zero_indices.n_elem > 0) {
        // If common zero rows/columns exist, create a cleaned matrix by removing them
        // Rcpp::Rcout << "Found all-zero rows and columns at (0-based) indices: " << common_zero_indices.t() << std::endl;
        
        arma::mat h_matrix_cleaned = h_matrix_arma;
        // Shed (remove) columns first, then rows, using the common indices
        // Note: shed_cols/shed_rows expect sorted unique indices, which intersect_uvec provides.
        h_matrix_cleaned.shed_cols(common_zero_indices); 
        h_matrix_cleaned.shed_rows(common_zero_indices); 
        
        // Rcpp::Rcout << "Cleaned h_matrix dimensions: " << h_matrix_cleaned.n_rows << "x" << h_matrix_cleaned.n_cols << std::endl;

        if (h_matrix_cleaned.n_rows > 0 && h_matrix_cleaned.is_square()) {
            double det_cleaned = arma::det(h_matrix_cleaned);
            if (std::abs(det_cleaned) > 1e-12) { // Check for non-singularity (tolerance for float comparison)
                try {
                    arma::mat inv_h_matrix_cleaned = arma::inv_sympd(h_matrix_cleaned);
                    alpha_se_value = inv_h_matrix_cleaned(0, 0);
                    // Rcpp::Rcout << "Recalculated alpha_se using cleaned h_matrix: " << alpha_se_value << std::endl;
                } catch (const std::runtime_error& e) {
                     Rcpp::warning("Inversion of cleaned h_matrix failed. Error: " + std::string(e.what()));
                     // alpha_se_value remains R_NaN
                }
            } else {
                Rcpp::warning("Cleaned h_matrix is singular (det ~ 0), cannot calculate alpha_se.");
            }
        } else {
            Rcpp::warning("Cleaned h_matrix is empty or not square after removing zero rows/cols.");
        }
    } else {
        // No common all-zero rows/columns, use the original h_matrix_arma
        // Rcpp::Rcout << "No common all-zero rows and columns found." << std::endl;
        if (h_matrix_arma.n_rows > 0 && h_matrix_arma.is_square()) {
             double det_original = arma::det(h_matrix_arma);
             if (std::abs(det_original) > 1e-12) { // Check for non-singularity
                try {
                    arma::mat inv_h_matrix_arma = arma::inv_sympd(h_matrix_arma);
                    alpha_se_value = inv_h_matrix_arma(0, 0);
                } catch (const std::runtime_error& e) {
                    Rcpp::warning("Inversion of original h_matrix failed. Error: " + std::string(e.what()));
                    // alpha_se_value remains R_NaN
                }
             } else {
                Rcpp::warning("Original h_matrix is singular (det ~ 0), cannot calculate alpha_se.");
             }
        } else {
             Rcpp::warning("Original h_matrix is empty or not square.");
        }
    }
    
    // Final step: take square root, ensuring value is non-negative and not NaN
    if (R_IsNaN(alpha_se_value) || alpha_se_value < 0) {
        if (!R_IsNaN(alpha_se_value) && alpha_se_value < 0) { // Only warn if it was a valid negative number
            Rcpp::warning("alpha_se_value is negative before sqrt, result will be NaN.");
        }
        return R_NaN; // Return NaN if value is already NaN or negative
    }

    return std::sqrt(alpha_se_value);
}