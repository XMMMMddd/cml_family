#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// 辅助函数：计算集合交集
arma::uvec intersect_uvec(const arma::uvec& a, const arma::uvec& b) {
    std::vector<arma::uword> result;
    for (arma::uword i = 0; i < a.n_elem; i++) {
        for (arma::uword j = 0; j < b.n_elem; j++) {
            if (a(i) == b(j)) {
                result.push_back(a(i));
                break;
            }
        }
    }
    return arma::conv_to<arma::uvec>::from(result);
}

// 辅助函数：计算集合差集
arma::uvec setdiff_uvec(const arma::uvec& a, const arma::uvec& b) {
    std::vector<arma::uword> result;
    for (arma::uword i = 0; i < a.n_elem; i++) {
        bool found = false;
        for (arma::uword j = 0; j < b.n_elem; j++) {
            if (a(i) == b(j)) {
                found = true;
                break;
            }
        }
        if (!found) {
            result.push_back(a(i));
        }
    }
    return arma::conv_to<arma::uvec>::from(result);
}

// 辅助函数：计算两个集合的并集
arma::uvec union_uvec(const arma::uvec& a, const arma::uvec& b) {
    std::set<arma::uword> result_set;
    for (arma::uword i = 0; i < a.n_elem; i++) {
        result_set.insert(a(i));
    }
    for (arma::uword i = 0; i < b.n_elem; i++) {
        result_set.insert(b(i));
    }
    
    std::vector<arma::uword> result_vec(result_set.begin(), result_set.end());
    return arma::conv_to<arma::uvec>::from(result_vec);
}

//' Core iterative MLE algorithm implemented in C++ with variable selection history
//' 
//' @param gamma_hat Matrix of gamma estimates (2 x num_snps)
//' @param beta_hat Matrix of beta estimates (2 x num_snps)
//' @param Sigma_gamma_list List of gamma covariance matrices
//' @param Sigma_beta_list List of beta covariance matrices  
//' @param Sigma_gamma_inv_list List of inverse gamma covariance matrices
//' @param Sigma_beta_inv_list List of inverse beta covariance matrices
//' @param a Size of set A
//' @param b Size of set B
//' @param max_iter Maximum number of iterations
//' @param tol Convergence tolerance
//' @param alpha_init Initial value for alpha
//' @return List containing results of the iterative algorithm with variable selection history
// [[Rcpp::export]]
List iterative_mle_core_with_history(
    arma::mat gamma_hat,
    arma::mat beta_hat,
    List Sigma_gamma_list,
    List Sigma_beta_list,
    List Sigma_gamma_inv_list,
    List Sigma_beta_inv_list,
    int a, int b,
    int max_iter = 100,
    double tol = 1e-6,
    double alpha_init = 0.0
) {
    int num_snps = gamma_hat.n_cols;
    
    // 初始化参数
    double alpha_current = alpha_init;
    arma::mat gamma_current = gamma_hat;
    
    bool converged = false;
    arma::vec alpha_history(max_iter, fill::zeros);
    arma::uvec set_A_indices, set_B_indices;
    
    // 新增：用于保存每次迭代的变量集历史
    std::vector<arma::uvec> set_A_history;
    std::vector<arma::uvec> set_B_history;
    std::vector<arma::uvec> set_A_intersect_B_history;
    std::vector<arma::uvec> set_A_diff_B_history;
    std::vector<arma::uvec> set_B_diff_A_history;
    std::vector<arma::vec> t_m_history;
    std::vector<arma::vec> f_m_history;
    
    int iter;
    for (iter = 0; iter < max_iter; iter++) {
        double alpha_prev = alpha_current;
        
        // ==========================================================
        // 步骤1: 选择SNPs集合A与B
        // ==========================================================
        arma::vec var_beta_o(num_snps);
        arma::vec var_beta_p(num_snps);
        
        for (int m = 0; m < num_snps; m++) {
            arma::mat Sigma_beta_m = as<arma::mat>(Sigma_beta_list[m]);
            var_beta_o(m) = Sigma_beta_m(0, 0);
            var_beta_p(m) = Sigma_beta_m(1, 1);
        }
        
        arma::vec t_m(num_snps);
        arma::vec f_m(num_snps);
        
        for (int m = 0; m < num_snps; m++) {
            t_m(m) = pow(beta_hat(0, m) - alpha_current * gamma_current(0, m), 2) / var_beta_o(m);
            f_m(m) = pow(beta_hat(1, m) - alpha_current * gamma_current(1, m), 2) / var_beta_p(m);
        }
        
        // 保存当前迭代的t_m和f_m值
        t_m_history.push_back(t_m);
        f_m_history.push_back(f_m);
        
        // 获取排序后的索引
        arma::uvec t_m_sorted = sort_index(t_m);
        arma::uvec f_m_sorted = sort_index(f_m);
        
        int a_actual = std::min(a, num_snps);
        int b_actual = std::min(b, num_snps);
        
        set_A_indices = t_m_sorted.head(a_actual);
        set_B_indices = f_m_sorted.head(b_actual);
        
        // 计算集合交集和差集
        arma::uvec A_intersect_B = intersect_uvec(set_A_indices, set_B_indices);
        arma::uvec A_diff_B = setdiff_uvec(set_A_indices, set_B_indices);
        arma::uvec B_diff_A = setdiff_uvec(set_B_indices, set_A_indices);
        
        // 保存当前迭代的变量集
        set_A_history.push_back(set_A_indices);
        set_B_history.push_back(set_B_indices);
        set_A_intersect_B_history.push_back(A_intersect_B);
        set_A_diff_B_history.push_back(A_diff_B);
        set_B_diff_A_history.push_back(B_diff_A);
        
        // ==========================================================
        // 步骤2: 更新gamma
        // ==========================================================
        arma::mat gamma_new = arma::zeros(2, num_snps);
        
        // 处理A∩B
        for (arma::uword i = 0; i < A_intersect_B.n_elem; i++) {
            int m = A_intersect_B(i);
            arma::mat term1_inv = as<arma::mat>(Sigma_gamma_inv_list[m]);
            arma::mat term2_inv = pow(alpha_current, 2) * as<arma::mat>(Sigma_beta_inv_list[m]);
            arma::mat combined_inv = term1_inv + term2_inv;
            arma::vec term1_val = term1_inv * gamma_hat.col(m);
            arma::vec term2_val = alpha_current * as<arma::mat>(Sigma_beta_inv_list[m]) * beta_hat.col(m);
            gamma_new.col(m) = solve(combined_inv, term1_val + term2_val);
        }
        
        // 处理A-B
        for (arma::uword i = 0; i < A_diff_B.n_elem; i++) {
            int m = A_diff_B(i);
            arma::mat Sigma_gamma_m = as<arma::mat>(Sigma_gamma_list[m]);
            
            double num_go = gamma_hat(0, m) * var_beta_o(m) - alpha_current * Sigma_gamma_m(0, 0) * beta_hat(0, m);
            double den_go = var_beta_o(m) - pow(alpha_current, 2) * Sigma_gamma_m(0, 0);
            
            // 检查分母是否接近零
            if (abs(den_go) < 1e-12) {
                Rcpp::warning("Denominator for gamma_o is near zero in A-B set");
                gamma_new.col(m) = gamma_hat.col(m);
                continue;
            }
            
            double gamma_o_new = num_go / den_go;
            
            double term_gp = Sigma_gamma_m(0, 1) * (-alpha_current / var_beta_o(m)) * 
                           (beta_hat(0, m) - alpha_current * gamma_o_new);
            double gamma_p_new = gamma_hat(1, m) + term_gp;
            
            gamma_new(0, m) = gamma_o_new;
            gamma_new(1, m) = gamma_p_new;
        }
        
        // 处理B-A
        for (arma::uword i = 0; i < B_diff_A.n_elem; i++) {
            int m = B_diff_A(i);
            arma::mat Sigma_gamma_m = as<arma::mat>(Sigma_gamma_list[m]);
            
            double num_gp = gamma_hat(1, m) * var_beta_p(m) - alpha_current * Sigma_gamma_m(1, 1) * beta_hat(1, m);
            double den_gp = var_beta_p(m) - pow(alpha_current, 2) * Sigma_gamma_m(1, 1);
            
            // 检查分母是否接近零
            if (abs(den_gp) < 1e-12) {
                Rcpp::warning("Denominator for gamma_p is near zero in B-A set");
                gamma_new.col(m) = gamma_hat.col(m);
                continue;
            }
            
            double gamma_p_new = num_gp / den_gp;
            
            double term_go = Sigma_gamma_m(0, 1) * (-alpha_current / var_beta_p(m)) * 
                           (beta_hat(1, m) - alpha_current * gamma_p_new);
            double gamma_o_new = gamma_hat(0, m) + term_go;
            
            gamma_new(0, m) = gamma_o_new;
            gamma_new(1, m) = gamma_p_new;
        }
        
        // 其他SNPs保持不变
        arma::uvec all_indices = linspace<arma::uvec>(0, num_snps-1, num_snps);
        arma::uvec selected_indices = union_uvec(set_A_indices, set_B_indices);
        arma::uvec others = setdiff_uvec(all_indices, selected_indices);
        
        for (arma::uword i = 0; i < others.n_elem; i++) {
            int m = others(i);
            gamma_new.col(m) = gamma_hat.col(m);
        }
        
        gamma_current = gamma_new;
        
        // ==========================================================
        // 步骤3: 更新alpha
        // ==========================================================
        double num_alpha = 0.0;
        double den_alpha = 0.0;
        
        // A∩B的贡献
        for (arma::uword i = 0; i < A_intersect_B.n_elem; i++) {
            int m = A_intersect_B(i);
            arma::mat Sigma_beta_inv_m = as<arma::mat>(Sigma_beta_inv_list[m]);
            num_alpha += as_scalar(gamma_current.col(m).t() * Sigma_beta_inv_m * beta_hat.col(m));
            den_alpha += as_scalar(gamma_current.col(m).t() * Sigma_beta_inv_m * gamma_current.col(m));
        }
        
        // A-B的贡献
        for (arma::uword i = 0; i < A_diff_B.n_elem; i++) {
            int m = A_diff_B(i);
            num_alpha += gamma_current(0, m) * beta_hat(0, m) / var_beta_o(m);
            den_alpha += pow(gamma_current(0, m), 2) / var_beta_o(m);
        }
        
        // B-A的贡献
        for (arma::uword i = 0; i < B_diff_A.n_elem; i++) {
            int m = B_diff_A(i);
            num_alpha += gamma_current(1, m) * beta_hat(1, m) / var_beta_p(m);
            den_alpha += pow(gamma_current(1, m), 2) / var_beta_p(m);
        }
        
        if (abs(den_alpha) > 1e-9) {
            alpha_current = num_alpha / den_alpha;
        } else {
            Rcpp::warning("Denominator for alpha is near zero. Halting.");
            break;
        }
        
        alpha_history(iter) = alpha_current;
        
        // 检查收敛
        if (abs(alpha_current - alpha_prev) < tol) {
            converged = true;
            iter++; // 调整迭代次数以反映实际完成的迭代数
            break;
        }
    }
    
    // 清理alpha_history，只保留实际使用的部分
    arma::vec alpha_history_clean = alpha_history.head(iter);
    
    // 将变量集历史转换为R格式（索引从1开始）
    List set_A_history_R(iter);
    List set_B_history_R(iter);
    List set_A_intersect_B_history_R(iter);
    List set_A_diff_B_history_R(iter);
    List set_B_diff_A_history_R(iter);
    List t_m_history_R(iter);
    List f_m_history_R(iter);
    
    for (int i = 0; i < iter; i++) {
        set_A_history_R[i] = set_A_history[i] + 1;
        set_B_history_R[i] = set_B_history[i] + 1;
        set_A_intersect_B_history_R[i] = set_A_intersect_B_history[i] + 1;
        set_A_diff_B_history_R[i] = set_A_diff_B_history[i] + 1;
        set_B_diff_A_history_R[i] = set_B_diff_A_history[i] + 1;
        t_m_history_R[i] = t_m_history[i];
        f_m_history_R[i] = f_m_history[i];
    }
    
    return List::create(
        Named("alpha_final") = alpha_current,
        Named("gamma_final") = gamma_current.t(),
        Named("iterations") = iter,
        Named("converged") = converged,
        Named("set_A_indices") = set_A_indices + 1, // 最终的集合A
        Named("set_B_indices") = set_B_indices + 1, // 最终的集合B
        Named("alpha_history") = alpha_history_clean,
        // 新增的历史记录
        Named("set_A_history") = set_A_history_R,
        Named("set_B_history") = set_B_history_R,
        Named("set_A_intersect_B_history") = set_A_intersect_B_history_R,
        Named("set_A_diff_B_history") = set_A_diff_B_history_R,
        Named("set_B_diff_A_history") = set_B_diff_A_history_R,
        Named("t_m_history") = t_m_history_R,
        Named("f_m_history") = f_m_history_R
    );
}