// File: iterative_update.cpp
// Save this as a .cpp file and compile using the R code below

#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>
#include <set>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Perform iterative update algorithm
//' 
//' @param n_snps Number of SNPs
//' @param beta_hat_exp Exposure beta estimates matrix
//' @param beta_hat_out Outcome beta estimates matrix  
//' @param alpha_input Initial alpha value
//' @param matrix_big Big covariance matrix
//' @param a_legal Number of legal SNPs for set a
//' @param b_legal Number of legal SNPs for set b
//' @param max_iter Maximum iterations (default 100)
//' @param tol Convergence tolerance (default 1e-6)
//' @return List containing results
//' @export
// [[Rcpp::export]]
List perform_iterative_update_rcpp(int n_snps, 
                                   arma::mat beta_hat_exp, 
                                   arma::mat beta_hat_out,
                                   double alpha_input, 
                                   arma::mat matrix_big, 
                                   int a_legal, 
                                   int b_legal,
                                   int max_iter = 100, 
                                   double tol = 1e-6) {
    
    // Initialize internal variables using input parameters
    arma::mat beta_exp = beta_hat_exp;
    double alpha = alpha_input;
    
    int iterations = 0;
    bool converged = false;
    
    // Variables to store final sets
    std::vector<int> final_complete_set;
    std::vector<int> final_a_reduce_b_set;
    
    for (int iter_count = 1; iter_count <= max_iter; iter_count++) {
        iterations = iter_count;
        double alpha_old = alpha; // Save alpha from previous iteration
        
        // Update r - Construct heterogeneity statistics
        arma::vec t(n_snps, fill::zeros);
        arma::vec f(n_snps, fill::zeros);
        
        for (int i = 0; i < n_snps; i++) {
            int matrix_indicator = 4 * (i + 1) - 1 - 1; // Convert to 0-based indexing
            
            // Calculate statistics for offspring
            double beta_hat_exp_i = beta_hat_exp(i, 0);
            double beta_exp_i = beta_exp(i, 0);
            double beta_out_se = matrix_big(matrix_indicator, 2); // 0-based indexing
            t(i) = pow(beta_hat_exp_i - alpha * beta_exp_i, 2) / beta_out_se;
            
            // Calculate statistics for parents
            beta_hat_exp_i = beta_hat_exp(i, 1);
            beta_exp_i = beta_exp(i, 1);
            beta_out_se = matrix_big(matrix_indicator + 1, 3); // 0-based indexing
            f(i) = pow(beta_hat_exp_i - alpha * beta_exp_i, 2) / beta_out_se;
        }
        
        // Create sorted indices for SNPs based on statistics
        arma::uvec a_indices = sort_index(t);
        arma::uvec b_indices = sort_index(f);
        
        // Get legal sets (first a_legal and b_legal indices)
        std::vector<int> a_legal_set, b_legal_set;
        for (int i = 0; i < std::min(a_legal, n_snps); i++) {
            a_legal_set.push_back(a_indices(i));
        }
        for (int i = 0; i < std::min(b_legal, n_snps); i++) {
            b_legal_set.push_back(b_indices(i));
        }
        
        // Find intersection and difference sets
        std::vector<int> complete_set;
        std::vector<int> a_reduce_b_set;
        
        // Find intersection
        std::sort(a_legal_set.begin(), a_legal_set.end());
        std::sort(b_legal_set.begin(), b_legal_set.end());
        std::set_intersection(a_legal_set.begin(), a_legal_set.end(),
                             b_legal_set.begin(), b_legal_set.end(),
                             std::back_inserter(complete_set));
        
        // Find difference (a_legal_set - b_legal_set)
        std::set_difference(a_legal_set.begin(), a_legal_set.end(),
                           b_legal_set.begin(), b_legal_set.end(),
                           std::back_inserter(a_reduce_b_set));
        
        // Store sets from current iteration
        final_complete_set = complete_set;
        final_a_reduce_b_set = a_reduce_b_set;
        
        // Update gamma (beta_exp)
        arma::mat beta_exp_new = beta_exp;
        
        // Update gamma for complete dataset
        if (!complete_set.empty()) {
            for (int idx : complete_set) {
                int i = idx; // SNP index (0-based)
                int matrix_indicator = 4 * (i + 1) - 3 - 1; // Convert to 0-based indexing
                
                try {
                    arma::mat matrix_big_i = matrix_big.rows(matrix_indicator, matrix_indicator + 3);
                    arma::mat matrix_big_i_inv = inv(matrix_big_i);
                    
                    arma::mat omega_11 = matrix_big_i_inv.submat(0, 0, 1, 1);
                    arma::mat omega_21 = matrix_big_i_inv.submat(2, 0, 3, 1);
                    arma::mat omega_12 = matrix_big_i_inv.submat(0, 2, 1, 3);
                    arma::mat omega_22 = matrix_big_i_inv.submat(2, 2, 3, 3);
                    
                    arma::rowvec beta_hat_exp_i = beta_hat_exp.row(i);
                    arma::rowvec beta_hat_out_i = beta_hat_out.row(i);
                    
                    arma::mat gamma_fenmu = omega_11 + alpha * (omega_21 + omega_12) + alpha * alpha * omega_22;
                    arma::mat gamma_fenzi = (omega_11 + alpha * omega_21) * beta_hat_exp_i.t() + 
                                           (omega_12 + alpha * omega_22) * beta_hat_out_i.t();
                    
                    arma::mat result = solve(gamma_fenmu, gamma_fenzi);
                    beta_exp_new.row(i) = result.t();
                } catch (const std::exception& e) {
                    Rcpp::warning("Matrix inversion failed for SNP " + std::to_string(i+1) + ": " + e.what());
                }
            }
        }
        
        // Update gamma for auxiliary dataset
        if (!a_reduce_b_set.empty()) {
            for (int idx : a_reduce_b_set) {
                int i = idx; // SNP index (0-based)
                int matrix_indicator = 4 * (i + 1) - 3 - 1; // Convert to 0-based indexing
                
                try {
                    arma::mat matrix_big_i = matrix_big.rows(matrix_indicator, matrix_indicator + 3);
                    arma::mat matrix_o_i(2, 2, fill::zeros);
                    
                    matrix_o_i(0, 0) = matrix_big_i(0, 0);
                    matrix_o_i(1, 1) = matrix_big_i(2, 2);
                    matrix_o_i(0, 1) = matrix_big_i(0, 2);
                    matrix_o_i(1, 0) = matrix_big_i(0, 2);
                    
                    arma::mat matrix_o_i_inv = inv(matrix_o_i);
                    double beta_hat_exp_i = beta_hat_exp(i, 0);
                    double beta_hat_out_i = beta_hat_out(i, 0);
                    
                    beta_exp_new(i, 1) = beta_hat_exp(i, 1);
                    
                    double fenzi = (matrix_o_i_inv(0, 0) + alpha * matrix_o_i_inv(1, 0)) * beta_hat_exp_i +
                                  (matrix_o_i_inv(0, 1) + alpha * matrix_o_i_inv(1, 1)) * beta_hat_out_i;
                    double fenmu = matrix_o_i_inv(0, 0) + alpha * matrix_o_i_inv(1, 0) +
                                  alpha * (matrix_o_i_inv(0, 1) + alpha * matrix_o_i_inv(1, 1));
                    
                    if (std::abs(fenmu) > 1e-10) {
                        beta_exp_new(i, 0) = fenzi / fenmu;
                    }
                } catch (const std::exception& e) {
                    Rcpp::warning("Matrix inversion failed for SNP " + std::to_string(i+1) + ": " + e.what());
                }
            }
        }
        
        beta_exp = beta_exp_new; // Update beta_exp for next iteration or final return
        
        // Update alpha
        double fenzi_complete = 0.0;
        double fenmu_complete = 0.0;
        double fenzi_reduce = 0.0;
        double fenmu_reduce = 0.0;
        
        if (!complete_set.empty()) {
            for (int idx : complete_set) {
                int i = idx; // SNP index (0-based)
                int matrix_indicator = 4 * (i + 1) - 3 - 1; // Convert to 0-based indexing
                
                try {
                    arma::mat matrix_big_i = matrix_big.rows(matrix_indicator, matrix_indicator + 3);
                    arma::mat matrix_big_i_inv = inv(matrix_big_i);
                    
                    arma::mat omega_11 = matrix_big_i_inv.submat(0, 0, 1, 1);
                    arma::mat omega_21 = matrix_big_i_inv.submat(2, 0, 3, 1);
                    arma::mat omega_12 = matrix_big_i_inv.submat(0, 2, 1, 3);
                    arma::mat omega_22 = matrix_big_i_inv.submat(2, 2, 3, 3);
                    
                    arma::rowvec beta_hat_exp_i = beta_hat_exp.row(i);
                    arma::rowvec beta_hat_out_i = beta_hat_out.row(i);
                    arma::rowvec beta_exp_i_current_iter = beta_exp.row(i);
                    
                    arma::mat temp1 = beta_hat_exp_i * (omega_21 * (beta_hat_exp_i.t() - beta_exp_i_current_iter.t()) + 
                                                       omega_22 * beta_hat_out_i.t());
                    arma::mat temp2 = beta_exp_i_current_iter * omega_22 * beta_exp_i_current_iter.t();
                    
                    fenzi_complete += temp1(0, 0);
                    fenmu_complete += temp2(0, 0);
                } catch (const std::exception& e) {
                    Rcpp::warning("Matrix operation failed for SNP " + std::to_string(i+1) + ": " + e.what());
                }
            }
        }
        
        if (!a_reduce_b_set.empty()) {
            for (int idx : a_reduce_b_set) {
                int i = idx; // SNP index (0-based)
                int matrix_indicator = 4 * (i + 1) - 3 - 1; // Convert to 0-based indexing
                
                try {
                    arma::mat matrix_big_i = matrix_big.rows(matrix_indicator, matrix_indicator + 3);
                    arma::mat matrix_o_i(2, 2, fill::zeros);
                    
                    matrix_o_i(0, 0) = matrix_big_i(0, 0);
                    matrix_o_i(1, 1) = matrix_big_i(2, 2);
                    matrix_o_i(0, 1) = matrix_big_i(0, 2);
                    matrix_o_i(1, 0) = matrix_big_i(0, 2);
                    
                    arma::mat matrix_o_i_inv = inv(matrix_o_i);
                    double beta_hat_exp_i = beta_hat_exp(i, 0);
                    double beta_hat_out_i = beta_hat_out(i, 0);
                    double beta_exp_i_current_iter = beta_exp(i, 0);
                    
                    fenzi_reduce += beta_hat_exp_i * (matrix_o_i_inv(1, 0) * 
                                                     (beta_hat_exp_i - beta_exp_i_current_iter) +
                                                     matrix_o_i_inv(1, 1) * beta_hat_out_i);
                    fenmu_reduce += beta_exp_i_current_iter * beta_exp_i_current_iter * matrix_o_i_inv(1, 1);
                } catch (const std::exception& e) {
                    Rcpp::warning("Matrix operation failed for SNP " + std::to_string(i+1) + ": " + e.what());
                }
            }
        }
        
        double alpha_new;
        if (std::abs(fenmu_complete + fenmu_reduce) < 1e-10) { // Avoid division by zero
            Rcpp::warning("Denominator is zero when updating alpha, iteration may not converge.");
            alpha_new = alpha; // Keep alpha unchanged or apply other error handling
        } else {
            alpha_new = (fenzi_complete + fenzi_reduce) / (fenmu_complete + fenmu_reduce);
        }
        
        // Check convergence
        if (std::abs(alpha_new - alpha_old) < tol) {
            converged = true;
            alpha = alpha_new; // Update alpha to final converged value
            break; // Exit loop
        }
        alpha = alpha_new; // Update alpha for next iteration
    }
    
    if (!converged) {
        Rcpp::warning("Algorithm did not reach convergence threshold " + std::to_string(tol) + 
                     " after " + std::to_string(max_iter) + " iterations");
    }
    
    // Convert final sets to R vectors (1-based indexing for R)
    IntegerVector r_complete_set(final_complete_set.size());
    IntegerVector r_a_reduce_b_set(final_a_reduce_b_set.size());
    
    for (size_t i = 0; i < final_complete_set.size(); i++) {
        r_complete_set[i] = final_complete_set[i] + 1; // Convert to 1-based indexing
    }
    for (size_t i = 0; i < final_a_reduce_b_set.size(); i++) {
        r_a_reduce_b_set[i] = final_a_reduce_b_set[i] + 1; // Convert to 1-based indexing
    }
    
    return List::create(
        Named("alpha_new") = alpha,
        Named("beta_exp_updated") = beta_exp,
        Named("iterations") = iterations,
        Named("converged") = converged,
        Named("complete_set") = r_complete_set,
        Named("a_reduce_b_set") = r_a_reduce_b_set
    );
}