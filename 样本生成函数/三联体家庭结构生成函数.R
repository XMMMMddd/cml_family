library(dplyr)
library(MASS)
# 辅助函数：创建相关矩阵
.create_correlation_matrix_r <- function(n_snps, correlation_type, correlation_param) {
  if (n_snps <= 0) {
    return(matrix(nrow = 0, ncol = 0))
  }
  if (n_snps == 1) {
    return(matrix(1, nrow = 1, ncol = 1))
  }

  if (correlation_type == "independent") {
    # 不相关的 SNPs，生成单位矩阵
    corr_matrix <- diag(n_snps)
  } else if (correlation_type == "equicorrelated") {
    # 所有 SNP 之间具有相同的相关系数 rho
    rho <- correlation_param
    if (n_snps > 1 && (rho < -1 / (n_snps - 1) - 1e-9 || rho > 1.0 + 1e-9)) { # 允许小的浮点误差
      warning(paste0(
        "对于 ", n_snps, " 个 SNPs 的等相关矩阵, rho (相关系数) = ", rho,
        " 可能导致矩阵非正定。",
        " 建议 rho 在 [", round(-1 / (n_snps - 1), 3), ", 1.0] 范围内。"
      ))
    }
    corr_matrix <- matrix(rho, nrow = n_snps, ncol = n_snps)
    diag(corr_matrix) <- 1.0
  } else if (correlation_type == "ar1") {
    # 一阶自回归 (AR1) 相关结构
    rho <- correlation_param
    if (abs(rho) > 1.0) {
      stop("对于 AR1 相关性, rho (相关系数) 必须在 -1.0 和 1.0 之间。")
    }
    indices <- 1:n_snps
    # 创建指数的绝对差矩阵
    abs_diff_indices <- abs(outer(indices, indices, "-"))
    corr_matrix <- rho^abs_diff_indices
  } else {
    stop(paste0(
      "未知的 correlation_type: ", correlation_type,
      "。有效类型为 'independent', 'equicorrelated', 'ar1'。"
    ))
  }
  return(corr_matrix)
}

# 创建 SNPs 的函数
#' @param correlation_type 可以取 "independent, equicorrelated ar1"
simulate_snps_r <- function(n_samples, n_snps, mafs = 0.1,
                            correlation_type = "independent", correlation_param = 0.5,
                            custom_corr_matrix = NULL) {
  # 处理 n_snps 为 0 的情况
  if (n_snps == 0) {
    return(matrix(NA, nrow = n_samples, ncol = 0))
  }

  # 1. 确定 MAFs
  actual_mafs <- NULL
  actual_mafs <- mafs

  # 2. 定义潜在变量的相关矩阵
  corr_matrix <- NULL
  if (n_snps == 1 && correlation_type != "custom") {
    corr_matrix <- matrix(1.0, nrow = 1, ncol = 1)
  } else if (correlation_type == "custom") {
    if (is.null(custom_corr_matrix)) {
      stop("如果 correlation_type 为 'custom', 则必须提供 custom_corr_matrix。")
    }
    if (!is.matrix(custom_corr_matrix) || nrow(custom_corr_matrix) != n_snps || ncol(custom_corr_matrix) != n_snps) {
      stop(paste0("custom_corr_matrix 必须是一个 (", n_snps, "x", n_snps, ") 的矩阵。"))
    }
    if (!isSymmetric(custom_corr_matrix, tol = 1e-6)) warning("custom_corr_matrix 不是对称的。")
    if (any(abs(diag(custom_corr_matrix) - 1) > 1e-9)) warning("custom_corr_matrix 的对角线元素不全为 1。")
    corr_matrix <- custom_corr_matrix
  } else {
    corr_matrix <- .create_correlation_matrix_r(n_snps, correlation_type, correlation_param)
  }

  # 3. 从多元正态分布中模拟潜在连续变量
  mean_vector <- rep(0, n_snps)
  latent_continuous_values <- matrix(NA, nrow = n_samples, ncol = n_snps) # 初始化

  if (n_snps > 0) {
    tryCatch(
      {
        latent_continuous_values <- MASS::mvrnorm(n = n_samples, mu = mean_vector, Sigma = corr_matrix, tol = 1e-6)
      },
      error = function(e) {
        stop(paste("从多元正态分布采样失败。相关矩阵可能不是正定的，或 mvrnorm 存在其他问题。原始错误:", e$message))
      }
    )
    # 如果 n_samples 是 1, mvrnorm 返回一个向量，需转换为矩阵
    if (n_samples == 1 && n_snps > 0 && is.vector(latent_continuous_values)) {
      latent_continuous_values <- matrix(latent_continuous_values, nrow = 1)
    }
  }


  # 4. 将连续值转换为基因型 (0, 1, 2)
  genotypes <- matrix(0, nrow = n_samples, ncol = n_snps)

  if (n_snps > 0) {
    for (j in 1:n_snps) {
      qj <- actual_mafs # 当前 SNP 的 MAF

      # 根据 HWE 计算基因型频率
      p_00 <- (1 - qj)^2 # 纯合主要等位基因 (基因型 0)
      p_01 <- 2 * qj * (1 - qj) # 杂合子 (基因型 1)
      # p_11 = qj^2 (纯合次要等位基因, 基因型 2)

      if (qj == 0) { # 所有个体都是纯合主要等位基因 (基因型 0)
        genotypes[, j] <- 0
        next # 跳到下一个 SNP
      }
      if (qj == 1) { # 所有个体都是纯合次要等位基因 (基因型 2)
        genotypes[, j] <- 2
        next # 跳到下一个 SNP
      }

      # 从标准正态分布 N(0,1) 计算阈值
      thresh1 <- qnorm(p_00)
      thresh2 <- qnorm(p_00 + p_01) # 等同于 qnorm(1 - qj^2)

      snp_latent_col <- latent_continuous_values[, j]
      genotypes[snp_latent_col <= thresh1, j] <- 0
      genotypes[snp_latent_col > thresh1 & snp_latent_col <= thresh2, j] <- 1
      genotypes[snp_latent_col > thresh2, j] <- 2
    }
  }
  return(genotypes)
}

# 一个样本的生成函数
#' @param compatibility_selection_geno 可以取 "independent, equicorrelated ar1"
generate_mr_trio_data_2sample <- function(
    N_exp = 1000, N_out = 1000, overlap_prop = 0,
    p_f = 0.3, p_m = 0.3,
    # Exposure Effects
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3,
    # Outcome Effects (Direct Pleiotropy / Dynastic)
    beta_FStoOE_out = 0, beta_MStoOE_out = 0,
    beta_OStoOE_out = 0,
    # Causal Effect
    beta_exp_to_out = 0.4,
    # Confounding Effects
    beta_confounding_exp = 0.2, beta_confounding_out = 0.2,
    # Other parameters
    correlation = 0.2, seed = NULL,
    # 选型婚配（基因）
    compatibility_selection_geno = "independent",
    correlation_param = 0.5,
    # 选型婚配（环境）
    compatibility_selection_factor_exp = 0,
    compatibility_selection_factor_out = 0,
    sample.outcome.betas.from.range = FALSE) { # 默认关闭范围采样

  # --- 0. 设置随机种子 ---
  set.seed(seed)
  N_total <- N_exp + N_out

  # 生成总个体ID或索引
  all_indices <- 1:N_total


  #' @title title snp生成函数
  #' @param N 样本量
  #' @param p 等位基因频率 description
  internal_generate_hwe_snps <- function(N, p) {
    if (!is.numeric(p) || length(p) != 1 || p < 0 || p > 1) {
      stop("内部错误: internal_generate_hwe_snps 的 'p' 必须是0到1之间的单个数值")
    }
    # 基因型频率: AA (2), Aa (1), aa (0)
    genotype_freqs <- c(p^2, 2 * p * (1 - p), (1 - p)^2)
    # 使用 R 的 sample 函数生成基因型
    # 注意： R的sample默认取整数1:x，所以显式提供基因型值
    genotypes <- sample(c(2, 1, 0), size = N, replace = TRUE, prob = genotype_freqs)
    return(genotypes)
  }

  #' @title (内部) 模拟从单个亲本传递的等位基因
  #' @param parent_snps 亲本基因型向量 (0, 1, 或 2)
  #' @return 长度相同的向量，包含传递的等位基因 (0 或 1)
  internal_get_transmitted_allele <- function(parent_snps) {
    transmitted_alleles <- numeric(length(parent_snps))
    for (i in seq_along(parent_snps)) {
      genotype <- parent_snps[i]
      if (genotype == 0) { # aa
        transmitted_alleles[i] <- 0
      } else if (genotype == 2) { # AA
        transmitted_alleles[i] <- 1
      } else if (genotype == 1) { # Aa
        # 0 或 1 各有50%概率
        transmitted_alleles[i] <- rbinom(1, 1, 0.5)
      } else {
        warning(paste("发现无效的亲本基因型:", genotype, "在索引", i, "- 将传递NA"))
        transmitted_alleles[i] <- NA # 或者可以停止执行 stop()
      }
    }
    return(transmitted_alleles)
  }

  # --- 2. 为 *所有 N_total* 个体生成 SNP 数据 ---
  Grand_Father_SNPs <- simulate_snps_r(
    n_samples = N_total, n_snps = 2, mafs = p_f,
    correlation_type = compatibility_selection_geno,
    correlation_param = correlation_param
  )
  Grandfather_Father_SNPs_all <- Grand_Father_SNPs[, 1]
  Grandmother_Father_SNPs_all <- Grand_Father_SNPs[, 2]
  Grand_Mother_SNPs <- simulate_snps_r(
    n_samples = N_total, n_snps = 2, mafs = p_m,
    correlation_type = compatibility_selection_geno,
    correlation_param = correlation_param
  )
  Grandfather_Mother_SNPs_all <- Grand_Mother_SNPs[, 1]
  Grandmother_Mother_SNPs_all <- Grand_Mother_SNPs[, 2]

  allele_from_GF_F_all <- internal_get_transmitted_allele(Grandfather_Father_SNPs_all)
  allele_from_GM_F_all <- internal_get_transmitted_allele(Grandmother_Father_SNPs_all)
  Father_SNPs_all <- allele_from_GF_F_all + allele_from_GM_F_all

  allele_from_GF_M_all <- internal_get_transmitted_allele(Grandfather_Mother_SNPs_all)
  allele_from_GM_M_all <- internal_get_transmitted_allele(Grandmother_Mother_SNPs_all)
  Mother_SNPs_all <- allele_from_GF_M_all + allele_from_GM_M_all

  Offspring_SNPs_Father_all <- internal_get_transmitted_allele(Father_SNPs_all)
  Offspring_SNPs_Mother_all <- internal_get_transmitted_allele(Mother_SNPs_all)
  Offspring_SNPs_all <- Offspring_SNPs_Father_all + Offspring_SNPs_Mother_all


  # --- 3. 为 *所有 N_total* 个体生成混杂因素和相关性因素 ---
  Confounder_exp_all <- rnorm(N_total, mean = 0, sd = 1)
  Confounder_out_all <- rnorm(N_total, mean = 0, sd = 1)
  correlation_safe <- max(0, min(1, correlation))
  if (correlation != correlation_safe) {
    warning(paste0("Correlation ", correlation, "修正为 ", correlation_safe))
  }

  # 家庭环境相关性
  correlation_factor_exp_all <- rnorm(N_total, mean = 0, sd = sqrt(correlation_safe))
  correlation_factor_out_all <- rnorm(N_total, mean = 0, sd = sqrt(correlation_safe))

  # 选型婚配相关性
  compatibility_selection_factor_exp <- rnorm(N_total, mean = 0, sd = sqrt(compatibility_selection_factor_exp))
  compatibility_selection_factor_out <- rnorm(N_total, mean = 0, sd = sqrt(compatibility_selection_factor_out))

  # --- 4. 为 *所有 N_total* 个体生成暴露数据 (连续) ---

  #' @title 暴露生成函数
  #' @param beta_StoE_exp SNP对暴露的影响 description
  #' @param beta_FStoOE_exp 父亲的 SNP 对暴露的影响
  #' @param beta_MStoOE_exp 母亲的 SNP 对暴露的影响
  #' @param SNPs_all 自己的 SNP
  #' @param SNPs_Father 父亲的 SNP
  #' @param SNPs_Mother 母亲的 SNP
  #' @param beta_confounding_exp 混杂的效应
  #' @param Confounder_all 混杂
  #' @param correlation_factor_all 相关性因子
  genrate_expose_function <- function(beta_StoE_exp, beta_FStoOE_exp, beta_MStoOE_exp,
                                      SNPs_all, SNPs_Father, SNPs_Mother,
                                      beta_confounding_exp, Confounder_all, correlation_factor_all, compatibility_selection_factor) {
    results <- beta_StoE_exp * SNPs_all + # 自己的效应
      beta_FStoOE_exp * SNPs_Father + # 来自父亲的效应
      beta_MStoOE_exp * SNPs_Mother + # 来自母亲的效应
      beta_confounding_exp * Confounder_all + correlation_factor_all +
      compatibility_selection_factor
    var_results <- var(results, na.rm = TRUE)
    exp_err <- sqrt(max(0, 1 - var_results))
    expose_all <- results + rnorm(N_total, mean = 0, sd = exp_err)
    return(expose_all)
  }

  Father_expose_all <- genrate_expose_function(
    beta_OStoOE_exp, beta_FStoOE_exp, beta_MStoOE_exp,
    Father_SNPs_all, Grandfather_Father_SNPs_all, Grandmother_Father_SNPs_all,
    beta_confounding_exp, Confounder_exp_all,
    correlation_factor_exp_all, compatibility_selection_factor_exp
  )
  Mother_expose_all <- genrate_expose_function(
    beta_OStoOE_exp, beta_FStoOE_exp, beta_MStoOE_exp,
    Mother_SNPs_all, Grandfather_Mother_SNPs_all, Grandmother_Mother_SNPs_all,
    beta_confounding_exp, Confounder_exp_all,
    correlation_factor_exp_all, compatibility_selection_factor_exp
  )
  Offspring_expose_all <- genrate_expose_function(
    beta_OStoOE_exp, beta_FStoOE_exp, beta_MStoOE_exp,
    Offspring_SNPs_all, Father_expose_all, Mother_expose_all,
    beta_confounding_exp, Confounder_exp_all, correlation_factor_exp_all, 0
  )

  # --- 5. 为 *所有 N_total* 个体生成结局数据 (连续) ---
  genrate_outcome_function <- function(beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
                                       expose_all, SNPs_all, SNPs_Father, SNPs_Mother,
                                       beta_confounding_out, Confounder_all,
                                       correlation_factor_all, compatibility_selection_factor) {
    outcome_deterministic_all <- beta_exp_to_out * expose_all +
      beta_OStoOE_out * SNPs_all +
      beta_FStoOE_out * SNPs_Father + beta_MStoOE_out * SNPs_Mother +
      beta_confounding_out * Confounder_all + correlation_factor_all + compatibility_selection_factor

    var_results <- var(outcome_deterministic_all, na.rm = TRUE)
    exp_err <- sqrt(max(0, 1 - var_results))
    out_all <- outcome_deterministic_all + rnorm(N_total, mean = 0, sd = exp_err)
    return(out_all)
  }

  Father_outcome_all <- genrate_outcome_function(
    beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
    Father_expose_all, Father_SNPs_all, Grandfather_Father_SNPs_all, Grandmother_Father_SNPs_all,
    beta_confounding_out, Confounder_out_all, correlation_factor_out_all, compatibility_selection_factor_out
  )
  Mother_outcome_all <- genrate_outcome_function(
    beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
    Mother_expose_all, Mother_SNPs_all, Grandfather_Mother_SNPs_all, Grandmother_Mother_SNPs_all,
    beta_confounding_out, Confounder_out_all, correlation_factor_out_all, compatibility_selection_factor_out
  )
  Offspring_outcome_all <- genrate_outcome_function(
    beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
    Offspring_expose_all, Offspring_SNPs_all, Father_SNPs_all, Mother_SNPs_all,
    beta_confounding_out, Confounder_out_all, correlation_factor_out_all, 0
  )


  # --- 6. 创建完整数据集 (包含所有信息，用于后续抽样) ---
  data_all <- data.frame(
    ID = all_indices, # 添加ID方便追踪
    Father_SNPs = Father_SNPs_all,
    Mother_SNPs = Mother_SNPs_all,
    Offspring_SNPs = Offspring_SNPs_all,
    Father_expose = Father_expose_all,
    Mother_expose = Mother_expose_all,
    Offspring_expose = Offspring_expose_all,
    Father_outcome = Father_outcome_all,
    Mother_outcome = Mother_outcome_all,
    Offspring_outcome = Offspring_outcome_all
  )

  # --- 7. 根据重叠比例抽取样本索引 ---
  # 重新计算重叠数，以防N_total_valid变小
  # 使用 floor 或 round 都可以，这里用 floor 确保不超过比例
  N_overlap_actual <- floor(min(N_exp, N_out) * overlap_prop)

  N_exp_only_final <- N_exp - N_overlap_actual

  N_out_only_final <- N_out - N_overlap_actual


  # 从有效索引中抽样
  indices_pool <- all_indices

  # 抽取重叠样本
  indices_overlap <- sample(indices_pool, size = N_overlap_actual, replace = FALSE)
  indices_pool <- setdiff(indices_pool, indices_overlap) # 从池中移除

  # 抽取暴露组独有样本
  indices_exp_only <- sample(indices_pool, size = N_exp_only_final, replace = FALSE)
  indices_pool <- setdiff(indices_pool, indices_exp_only)

  # 抽取结局组独有样本
  indices_out_only <- sample(indices_pool, size = N_out_only_final, replace = FALSE)

  # 组合最终索引
  indices_exp <- c(indices_overlap, indices_exp_only)
  indices_out <- c(indices_overlap, indices_out_only)

  # --- 8. 创建最终的暴露组和结局组数据框 ---
  exposure_data <- data_all %>%
    filter(ID %in% indices_exp) %>%
    dplyr::select(
      ID, # <-- 添加 ID 列
      Father_SNPs, Mother_SNPs, Offspring_SNPs,
      Father_expose, Mother_expose, Offspring_expose
    )

  outcome_data <- data_all %>%
    filter(ID %in% indices_out) %>%
    dplyr::select(
      ID, # <-- 添加 ID 列
      Father_SNPs, Mother_SNPs, Offspring_SNPs,
      Father_outcome, Mother_outcome, Offspring_outcome
    )

  # --- 9. 返回列表 ---
  return(list(exposure_data = as.data.frame(exposure_data), outcome_data = as.data.frame(outcome_data)))
}



# 3 种数据集类型生成函数 ----
#' @param compatibility_selection_geno 可以取 "independent, equicorrelated ar1"
generate_multiple_datasets_v3 <- function(
    n = 10, # 要生成的总数据集数量
    num_pleiotropic = 1, # 指定多少个数据集具有 *潜在的* 非零结局效应
    N_exp = 1000, N_out = 1000, overlap_prop = 0,
    p_f = 0.3, p_m = 0.3,
    # --- 暴露效应 (当非零时的大小) ---
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3,
    # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
    # 子代基因型 -> 结局
    mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0.05, # 新增: 均值和标准差
    mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0.05, # 新增: 均值和标准差
    # 父母基因型 -> 结局
    mean_beta_OStoOE_out = 0.1, sd_beta_OStoOE_out = 0.05, # 新增: 均值和标准差
    prop_negative_pleiotropy = 0.5, # 新增: 在多效性 SNP 中，效应为负的比例 (0 到 1)
    # 选型婚配
    compatibility_selection_prop = 0, # 选型婚配在数据集中的比例
    compatibility_selection_geno = "independent",
    correlation_param = 0.5,
    compatibility_selection_factor_exp = 0,
    compatibility_selection_factor_out = 0,
    # 人群分层
    ## 定义人群分层的差异(次等位基因频率差异)
    crowd_stratification_differences = 0, # 两个人群
    # --- 其他效应 ---
    beta_exp_to_out = 0.4,
    beta_confounding_exp = 0.2,
    beta_confounding_out = 0.2,
    correlation = 0.2,
    seed = NULL) {
  # --- 1. 输入验证与随机种子设置 ---
  if (!is.null(seed)) {
    set.seed(seed)
  }

  n_pool <- 1:n
  pleiotropic_indices <- sample(n_pool, num_pleiotropic)
  not_pleiotropic_indices <- setdiff(n_pool, pleiotropic_indices)


  # --- 3. 循环生成 n 个数据集 ---
  all_datasets <- vector("list", n)

  for (i in 1:n) {
    # --- 3a. 判断当前数据集 i 的特性 ---
    is_pleiotropic_candidate <- i %in% pleiotropic_indices # 是否是潜在的多效性 SNP

    # --- 3b. 设置当前循环使用的 beta 参数值 ---

    # --- 暴露效应 ---
    current_beta_FStoOE_exp <- beta_FStoOE_exp
    current_beta_MStoOE_exp <- beta_MStoOE_exp
    current_beta_OStoOE_exp <- beta_OStoOE_exp


    # --- 结局效应 / 水平多效性 (核心修改) ---
    if (is_pleiotropic_candidate) {
      # 决定这个多效性 SNP 的效应方向 (正或负)
      # rbinom(1, 1, prob) 返回 0 或 1，1 代表负向
      direction_multiplier <- if (runif(1) < prop_negative_pleiotropy) -1 else 1

      # 为这个特定的 SNP i 抽样生成多效性效应大小
      # 从以 (方向 * 均值) 为中心的正态分布中抽样
      current_beta_FStoOE_out <- rnorm(1, mean = direction_multiplier * mean_beta_FStoOE_out, sd = sd_beta_FStoOE_out)
      current_beta_MStoOE_out <- rnorm(1, mean = direction_multiplier * mean_beta_MStoOE_out, sd = sd_beta_MStoOE_out)
      current_beta_OStoOE_out <- rnorm(1, mean = direction_multiplier * mean_beta_OStoOE_out, sd = sd_beta_OStoOE_out)

      snp_label <- "pleiotropic_variable" # 标记为可变的多效性
    } else {
      # 如果不是多效性 SNP，效应为 0
      current_beta_FStoOE_out <- 0
      current_beta_MStoOE_out <- 0
      current_beta_OStoOE_out <- 0

      # 确定非多效性 SNP 的标签
      snp_label <- "valid_instrument" # 有效工具变量
    }
    # 选型婚配比例
    is_compatibility_selection <- if (runif(1) < compatibility_selection_prop) 1 else 0
    if (is_compatibility_selection == 1) {
      current_compatibility_selection_geno <- compatibility_selection_geno
      current_correlation_param <- correlation_param
      current_compatibility_selection_factor_exp <- compatibility_selection_factor_exp
      current_compatibility_selection_factor_out <- compatibility_selection_factor_out
    } else {
      current_compatibility_selection_geno <- "independent"
      current_correlation_param <- 0
      current_compatibility_selection_factor_exp <- 0
      current_compatibility_selection_factor_out <- 0
    }

    # 人群分层
    current_p_f <- sample(c(p_f, p_f + crowd_stratification_differences), 1)
    current_p_m <- sample(c(p_m, p_m + crowd_stratification_differences), 1)

    # --- 3c. 调用底层函数生成单个数据集 ---
    data_i <- generate_mr_trio_data_2sample(
      N_exp = N_exp, N_out = N_out,
      overlap_prop = overlap_prop,
      p_f = current_p_f, p_m = current_p_m,
      beta_FStoOE_exp = current_beta_FStoOE_exp,
      beta_MStoOE_exp = current_beta_MStoOE_exp,
      beta_OStoOE_exp = current_beta_OStoOE_exp,
      # 传递当前循环抽样生成的或设为 0 的结局效应值
      beta_FStoOE_out = current_beta_FStoOE_out,
      beta_MStoOE_out = current_beta_MStoOE_out,
      beta_OStoOE_out = current_beta_OStoOE_out,
      beta_exp_to_out = beta_exp_to_out,
      beta_confounding_exp = beta_confounding_exp,
      beta_confounding_out = beta_confounding_out,
      correlation = correlation, seed = NULL,
      # 选型婚配参数
      compatibility_selection_geno = current_compatibility_selection_geno,
      correlation_param <- current_correlation_param,
      compatibility_selection_factor_exp =
        current_compatibility_selection_factor_exp,
      compatibility_selection_factor_out =
        current_compatibility_selection_factor_out
    )

    # --- 3d. 为生成的数据添加标识和类型标签 (使用上面确定的 snp_label) ---
    if (length(data_i) >= 2 && is.data.frame(data_i[[1]]) && is.data.frame(data_i[[2]])) {
      data_i[[1]]$dataset_id <- i
      data_i[[2]]$dataset_id <- i
      # 添加类型标签
      data_i[[1]]$snp_type <- snp_label
      data_i[[2]]$snp_type <- snp_label
    } else {
      warning(paste("数据集", i, "的生成结果格式不符合预期，无法添加 ID 和标签。"))
    }

    all_datasets[[i]] <- data_i
  } # 结束 for 循环

  # --- 4. 合并所有数据集 (与 v2 相同) ---
  all_exposure_dfs <- lapply(all_datasets, function(element) {
    if (length(element) >= 1 && is.data.frame(element[[1]])) {
      return(element[[1]])
    } else {
      return(NULL)
    }
  })
  all_outcome_dfs <- lapply(all_datasets, function(element) {
    if (length(element) >= 2 && is.data.frame(element[[2]])) {
      return(element[[2]])
    } else {
      return(NULL)
    }
  })
  combined_exposure_data <- do.call(rbind, all_exposure_dfs)
  combined_outcome_data <- do.call(rbind, all_outcome_dfs)


  # --- [可选] 重命名结局数据框列名 (与 v2 相同) ---
  if (!is.null(combined_outcome_data) && nrow(combined_outcome_data) > 0) {
    current_names <- names(combined_outcome_data)
    new_names <- sub("_outcome$", "_expose", current_names) # 注意: 这里可能需要调整以匹配你的底层函数输出
    # 检查列名是否存在再重命名，防止错误
    cols_to_rename <- intersect(current_names, grep("_outcome$", current_names, value = TRUE))
    names(combined_outcome_data)[match(cols_to_rename, names(combined_outcome_data))] <- sub("_outcome$", "_expose", cols_to_rename)
  }


  # --- 5. 返回列表 ---
  return(list(combined_exposure_data, combined_outcome_data))
}
