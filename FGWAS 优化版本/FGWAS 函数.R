# %% 需要用到的包
library(lme4)
library(dplyr)
library(MASS)
library(nlme) # <-- 已添加
library(parallel)
library(pbapply)
# %% FGWAS 多个数据

# lmm但是只有一个变量
FMR_trio_optimized_glsfunction_single <- function(data_test) {
  n <- nrow(data_test) # 获取样本量

  # --- 1. 计算协方差矩阵 Omega ---
  # 使用 colMeans/cov/var 提高效率和简洁性
  exposures <- data_test[, c("Father_expose", "Mother_expose", "Offspring_expose")]
  father_exp <- exposures[[1]]
  mother_exp <- exposures[[2]]
  offspring_exp <- exposures[[3]]


  # --- 2. 准备基因型矩阵 X 和 X_par ---
  # 提取 SNP 数据
  snps <- data_test[, c("Father_SNPs", "Mother_SNPs", "Offspring_SNPs")]
  father_snps <- snps[[1]]
  mother_snps <- snps[[2]]
  offspring_snps <- snps[[3]]


  data_test_father <- data_test %>% dplyr::select(starts_with("Father"))
  data_test_father$family_id <- 1:nrow(data_test_father)
  data_test_mother <- data_test %>% dplyr::select(starts_with("Mother"))
  data_test_mother$family_id <- 1:nrow(data_test_mother)
  data_test_offspring <- data_test %>% dplyr::select(starts_with("Offspring"))
  data_test_offspring$family_id <- 1:nrow(data_test_offspring)
  names(data_test_father) <- c("snps", "expose", "family_id")
  names(data_test_mother) <- c("snps", "expose", "family_id")
  names(data_test_offspring) <- c("snps", "expose", "family_id")
  data_long <- bind_rows(data_test_father, data_test_mother, data_test_offspring)

  gls_model_gen <- gls(expose ~ snps,
    data = data_long,
    correlation = corSymm(form = ~ 1 | family_id),
    method = "ML"
  )
  informatrix <- solve(vcov(gls_model_gen)[2, 2])


  # --- 5. 返回结果 ---
  # 返回与原函数相同的结构
  A <- list(
    beta = matrix(gls_model_gen$coefficients[2]), # 1x2 矩阵
    informatrix = informatrix,
    beta_se = sqrt(1 / informatrix)
  )

  return(A)
}

# lmm两个变量
FMR_trio_optimized_glsfunction <- function(data_test) {
  n <- nrow(data_test) # 获取样本量

  # --- 1. 计算协方差矩阵 Omega ---
  # 使用 colMeans/cov/var 提高效率和简洁性
  exposures <- data_test[, c("Father_expose", "Mother_expose", "Offspring_expose")]
  father_exp <- exposures[[1]]
  mother_exp <- exposures[[2]]
  offspring_exp <- exposures[[3]]



  # --- 2. 准备基因型矩阵 X 和 X_par ---
  # 提取 SNP 数据
  snps <- data_test[, c("Father_SNPs", "Mother_SNPs", "Offspring_SNPs")]
  father_snps <- snps[[1]]
  mother_snps <- snps[[2]]
  offspring_snps <- snps[[3]]

  # 计算等位基因频率 p
  # 使用 mean 直接计算所有 SNP 的均值，然后除以 2
  ## 父母只用一个防止人群分层
  p <- mean(c(mother_snps, offspring_snps)) / 2

  # 计算 G1, G2, G3 (向量化)
  # 注意: 原始 G1/G2 公式的含义可能需要确认，这里按原样实现
  G1 <- father_snps + 2 * p
  G2 <- mother_snps + 2 * p
  G3 <- father_snps + mother_snps


  data_test_father <- data_test %>% dplyr::select(starts_with("Father"))
  data_test_father$family_id <- 1:nrow(data_test_father)
  data_test_father$g_snps <- G1
  data_test_mother <- data_test %>% dplyr::select(starts_with("Mother"))
  data_test_mother$family_id <- 1:nrow(data_test_mother)
  data_test_mother$g_snps <- G2
  data_test_offspring <- data_test %>% dplyr::select(starts_with("Offspring"))
  data_test_offspring$family_id <- 1:nrow(data_test_offspring)
  data_test_offspring$g_snps <- G3
  names(data_test_father) <- c("snps", "expose", "family_id", "g_snps")
  names(data_test_mother) <- c("snps", "expose", "family_id", "g_snps")
  names(data_test_offspring) <- c("snps", "expose", "family_id", "g_snps")
  data_long <- arrange(bind_rows(
    data_test_father,
    data_test_mother, data_test_offspring
  ), family_id)

  gls_model_gen <- gls(expose ~ snps + g_snps,
    data = data_long,
    correlation = corSymm(form = ~ 1 | family_id),
    method = "ML"
  )
  informatrix <- solve(vcov(gls_model_gen)[2:3, 2:3])


  # --- 5. 返回结果 ---
  # 返回与原函数相同的结构
  A <- list(
    beta = gls_model_gen$coefficients[2:3], # 1x2 矩阵
    informatrix = informatrix
  )

  return(A)
}

# lmm两个变量(FGWAS)
FMR_trio_optimized <- function(data_test) {
  n <- nrow(data_test) # 获取样本量

  # --- 1. 计算协方差矩阵 Omega ---
  # 使用 colMeans/cov/var 提高效率和简洁性
  exposures <- data_test[, c("Father_expose", "Mother_expose", "Offspring_expose")]
  father_exp <- exposures[[1]]
  mother_exp <- exposures[[2]]
  offspring_exp <- exposures[[3]]

  # 计算协方差/方差 （注意：cov计算的是样本协方差，除以 n-1）
  # 如果原始代码的意图是总体协方差（除以 n），需要调整
  sa2 <- (cov(father_exp, offspring_exp) + cov(mother_exp, offspring_exp)) / 2
  ss2 <- cov(father_exp, mother_exp)
  se2_F <- var(father_exp)
  se2_M <- var(mother_exp)
  se2_O <- var(offspring_exp)

  Sigma <- matrix(c(
    se2_F, ss2,   sa2,
    ss2,   se2_M, sa2,
    sa2,   sa2,   se2_O
  ), nrow = 3, byrow = TRUE)

  # 计算 Omega 的逆
  # 使用 tryCatch 处理可能的奇异矩阵
  Omega_inv <- tryCatch(
    {
      chol2inv(chol(Sigma))
    },
    error = function(e) {
      warning("协方差矩阵 Sigma 求逆失败 (可能奇异): ", e$message)
      # 返回一个 NaN 矩阵或根据需要停止
      matrix(NaN, 3, 3)
    }
  )
  if (any(is.nan(Omega_inv))) {
    return(list(beta = matrix(NaN, 1, 2), informatrix = matrix(NaN, 2, 2), beta_se = NaN))
  }


  # --- 2. 准备基因型矩阵 X 和 X_par ---
  # 提取 SNP 数据
  snps <- data_test[, c("Father_SNPs", "Mother_SNPs", "Offspring_SNPs")]
  father_snps <- snps[[1]]
  mother_snps <- snps[[2]]
  offspring_snps <- snps[[3]]

  # 计算等位基因频率 p
  # 使用 mean 直接计算所有 SNP 的均值，然后除以 2
  ## 父母只用一个防止人群分层
  p <- mean(c(mother_snps, offspring_snps)) / 2

  # 计算 G1, G2, G3 (向量化)
  # 注意: 原始 G1/G2 公式的含义可能需要确认，这里按原样实现
  G1 <- ifelse(father_snps == 0, 2 * p,
    ifelse(father_snps == 1, 1 + 2 * p, 2 * (1 + p))
  )
  G2 <- ifelse(mother_snps == 0, 2 * p,
    ifelse(mother_snps == 1, 1 + 2 * p, 2 * (1 + p))
  )
  G3 <- father_snps + mother_snps

  # 构建基因型矩阵
  X <- as.matrix(snps) # n x 3
  X_par <- cbind(G1, G2, G3) # n x 3
  Y <- as.matrix(exposures) # n x 3


  # --- 3. 向量化计算 fenzi 和 informatrix (原 fenmu) ---
  # 预计算常用乘积
  Y_OmInv <- Y %*% Omega_inv # n x 3
  X_OmInv <- X %*% Omega_inv # n x 3
  Xp_OmInv <- X_par %*% Omega_inv # n x 3

  # 计算 fenzi (1x2 矩阵) 的两个元素
  # sum(Y_OmInv * X) 等价于 sum(diag(t(Y) %*% Omega_inv %*% X))
  fenzi_1 <- sum(Y_OmInv * X)
  fenzi_2 <- sum(Y_OmInv * X_par)
  fenzi <- matrix(c(fenzi_1, fenzi_2), nrow = 1)

  # 计算 informatrix (2x2 矩阵) 的元素 (原 fenmu)
  # sum(X_OmInv * X) 等价于 sum(diag(t(X) %*% Omega_inv %*% X))
  info_11 <- sum(X_OmInv * X) # 原 beta_se 累加值
  info_12 <- sum(X_OmInv * X_par) # 原 beta_cov 累加值 (info_21 = info_12)
  info_22 <- sum(Xp_OmInv * X_par) # 原 beta_par_se 累加值

  informatrix <- matrix(c(
    info_11, info_12,
    info_12, info_22
  ), nrow = 2, byrow = TRUE)


  # --- 4. 计算 Beta 和标准误 ---
  # 计算 informatrix 的逆
  informatrix_inv <- tryCatch(
    {
      chol2inv(chol(informatrix))
    },
    error = function(e) {
      warning("信息矩阵 informatrix 求逆失败 (可能奇异): ", e$message)
      matrix(NaN, 2, 2)
    }
  )

  if (any(is.nan(informatrix_inv))) {
    beta <- matrix(NaN, 1, 2)
    beta_se1 <- NaN
  } else {
    # 计算 beta (1x2 矩阵)
    beta <- fenzi %*% informatrix_inv

    # 计算第一个 beta 的标准误
    # 标准误是方差的平方根，方差是协方差矩阵（即 informatrix_inv）的对角元素
    beta_se1 <- sqrt(pmax(informatrix_inv[1, 1], 0)) # 使用 pmax 避免对负数开方
  }

  # --- 5. 返回结果 ---
  # 返回与原函数相同的结构
  A <- list(
    beta = beta, # 1x2 矩阵
    informatrix = informatrix, # 2x2 信息矩阵
    beta_se = beta_se1 # 第一个 beta 的标准误 (标量)
    # 如果需要第二个beta的SE: sqrt(pmax(informatrix_inv[2, 2], 0))
    # 如果需要完整的beta协方差矩阵: informatrix_inv
  )

  return(A)
}


# 只有孩子的暴露被使用(加上父母当协变量)
FMR_trio_optimized_onlyc <- function(data_test) {
  n <- nrow(data_test) # 获取样本量

  # --- 1. 计算协方差矩阵 Omega ---
  # 使用 colMeans/cov/var 提高效率和简洁性
  exposures <- data_test[, c("Father_expose", "Mother_expose", "Offspring_expose")]

  data_test$g_snps <- data_test$Father_SNPs + data_test$Mother_SNPs
  model_ml <- summary(lm(Offspring_expose ~ Offspring_SNPs +
    g_snps, data = data_test))
  coefficients_table <- model_ml$coefficients
  # --- 5. 返回结果 ---
  # 返回修改后的结构
  A <- list(
    beta = matrix(coefficients_table[, "Estimate"][2],
      nrow = 1, ncol = 1
    ), # 1x1 矩阵
    informatrix = matrix(1 / model_ml$coefficients[, "Std. Error"][2],
      nrow = 1, ncol = 1
    ), # 1x1 信息矩阵
    beta_se = model_ml$coefficients[, "Std. Error"][2] # beta 的标准误 (标量)
  )

  return(A)
}

# 只有孩子的暴露被使用(父母的变量当做协变量)
FMR_trio_optimized_onlyc2 <- function(data_test) {
  required_cols <- c(
    "Father_expose", "Mother_expose", "Offspring_expose",
    "Father_SNPs", "Mother_SNPs", "Offspring_SNPs"
  )

  missing_cols <- setdiff(required_cols, names(data_test))
  if (length(missing_cols) > 0) {
    stop("缺少必要的列: ", paste(missing_cols, collapse = ", "))
  }

  # 创建g_snps变量
  data_test$g_snps <- data_test$Father_SNPs + data_test$Mother_SNPs

  # 拟合线性回归模型 - 保留原始lm对象
  lm_model <- lm(Offspring_expose ~ Offspring_SNPs + g_snps, data = data_test)
  model_summary <- summary(lm_model) # 分离summary对象

  # 获取系数表
  coefficients_table <- model_summary$coefficients

  # 正确计算信息矩阵 - 使用原始lm对象
  vcov_matrix <- vcov(lm_model) # 从lm对象获取方差-协方差矩阵
  information_matrix_full <- solve(vcov_matrix) # 完整的信息矩阵

  # 提取非截距项的信息矩阵 (2x2矩阵，对应Offspring_SNPs和g_snps)
  information_matrix <- information_matrix_full[2:3, 2:3]

  # 提取系数估计值 (去除截距项)
  beta_estimates <- coefficients_table[2:3, "Estimate"]

  # 提取标准误 (两个参数的标准误)
  beta_se <- coefficients_table[2:3, "Std. Error"]

  # 返回结果
  result <- list(
    beta = matrix(beta_estimates, nrow = 1, ncol = 2), # 1x2 矩阵
    informatrix = information_matrix, # 2x2 信息矩阵 (不是1x1!)
    beta_se = beta_se # 长度为2的向量 (两个参数的标准误)
  )
  return(result)
}

# 只有孩子的暴露被使用(且没有其他当做协变量)
FMR_trio_optimized_onlyc3 <- function(data_test) {
  n <- nrow(data_test) # 获取样本量



  model_ml <- summary(lm(Offspring_expose ~ Offspring_SNPs, data = data_test))
  coefficients_table <- model_ml$coefficients
  # --- 5. 返回结果 ---
  # 返回修改后的结构
  A <- list(
    beta = matrix(coefficients_table[, "Estimate"][2],
      nrow = 1, ncol = 1
    ), # 1x1 矩阵
    informatrix = matrix(1 / model_ml$coefficients[, "Std. Error"][2],
      nrow = 1, ncol = 1
    ), # 1x1 信息矩阵
    beta_se = model_ml$coefficients[, "Std. Error"][2] # beta 的标准误 (标量)
  )

  return(A)
}

# 对每一个数据集进行 fgwas
fgwas_for_data_optimized <- function(data_df,
                                     grouping_col = "dataset_id",
                                     processing_func = FMR_trio_optimized) { # 假设 FMR_trio 存在

  # 获取唯一的数据集 ID
  unique_ids <- sort(unique(data_df[[grouping_col]])) # 排序确保结果顺序一致
  n_datasets <- length(unique_ids)


  # --- 初始化结果矩阵 ---
  # 使用 NA 初始化，以便于后续检查是否所有都成功计算
  beta_hat_matrix <- matrix(NA_real_, nrow = n_datasets, ncol = 2)
  # 使用数据集 ID 作为行名，方便追踪
  rownames(beta_hat_matrix) <- unique_ids
  colnames(beta_hat_matrix) <- c("beta_col1", "beta_col2") # 或更具体的名称

  Sigma_inv_matrix <- matrix(NA_real_, nrow = 2 * n_datasets, ncol = 2)

  cat(paste("正在处理", n_datasets, "个数据集...\n"))

  # --- 遍历数据集并应用处理函数 ---
  for (i in 1:n_datasets) {
    current_id <- unique_ids[i]

    # 筛选当前数据集的数据 (使用 base R，通常在循环中略快)
    subset_data <- data_df[data_df[[grouping_col]] == current_id, , drop = FALSE]


    processed_result <- tryCatch(
      {
        processing_func(subset_data)
      },
      error = function(e) {
        warning(paste0("处理数据集 '", current_id, "' 时出错: ", e$message))
        return(NULL) # 返回 NULL 表示失败
      }
    )

    # 存储 beta 估计值
    beta_hat_matrix[i, ] <- processed_result[[1]]
    # 存储 2x2 的 Sigma_inv 块
    row_indices <- (2 * i - 1):(2 * i)
    Sigma_inv_matrix[row_indices, ] <- processed_result[[2]]

    # beta_hat_matrix[i,] 和 Sigma_inv_matrix[row_indices,] 保持 NA
  } # 结束 for 循环

  cat("数据集处理完成。\n")

  # --- 返回结果 ---
  # 检查是否有任何数据集处理成功

  final_result <- list(beta_hat = beta_hat_matrix, Sigma_inv = Sigma_inv_matrix)
  return(final_result)
}

fgwas_for_data_optimized_single <- function(data_df,
                                            grouping_col = "dataset_id",
                                            processing_func = FMR_trio_optimized_onlyc3) { # 假设 FMR_trio 存在

  # 获取唯一的数据集 ID
  unique_ids <- sort(unique(data_df[[grouping_col]])) # 排序确保结果顺序一致
  n_datasets <- length(unique_ids)


  # --- 初始化结果矩阵 ---
  # 使用 NA 初始化，以便于后续检查是否所有都成功计算
  beta_hat <- matrix(NA_real_, nrow = n_datasets, ncol = 1)
  beta_se <- matrix(NA_real_, nrow = n_datasets, ncol = 1)

  cat(paste("正在处理", n_datasets, "个数据集...\n"))

  # --- 遍历数据集并应用处理函数 ---
  for (i in 1:n_datasets) {
    current_id <- unique_ids[i]

    # 筛选当前数据集的数据 (使用 base R，通常在循环中略快)
    subset_data <- data_df[data_df[[grouping_col]] == current_id, , drop = FALSE]


    processed_result <- tryCatch(
      {
        processing_func(subset_data)
      },
      error = function(e) {
        warning(paste0("处理数据集 '", current_id, "' 时出错: ", e$message))
        return(NULL) # 返回 NULL 表示失败
      }
    )

    # 存储 beta 估计值
    beta_hat[i] <- processed_result[[1]]

    beta_se[i] <- processed_result[[3]]

    # beta_hat_matrix[i,] 和 Sigma_inv_matrix[row_indices,] 保持 NA
  } # 结束 for 循环

  cat("数据集处理完成。\n")

  # --- 返回结果 ---
  # 检查是否有任何数据集处理成功

  final_result <- list(beta_hat = beta_hat, beta_se = beta_se)
  return(final_result)
}


# %% FGWAS 简单版
fmr_trio_optimized <- function(data_test, snp_index = 1) {
  # 参数说明:
  # data_test: 输入数据框
  # snp_index: 指定使用哪个SNP位点进行分析 (默认为1，即第一个SNP)

  n <- nrow(data_test) # 获取样本量

  # --- 1. 数据预处理和列名匹配 ---
  # 处理exposure列名的大小写差异
  exposure_cols <- c()
  if ("father_expose" %in% names(data_test)) {
    exposure_cols <- c("father_expose", "mother_expose", "offspring_expose")
  } else if ("Father_expose" %in% names(data_test)) {
    exposure_cols <- c("Father_expose", "Mother_expose", "Offspring_expose")
  } else {
    stop("未找到exposure列，请检查列名")
  }

  # 构建SNP列名
  snp_cols <- c(
    paste0("father_snps.", snp_index),
    paste0("mother_snps.", snp_index),
    paste0("offspring_snps.", snp_index)
  )

  # 检查SNP列是否存在
  if (!all(snp_cols %in% names(data_test))) {
    stop(paste("未找到SNP列:", paste(snp_cols, collapse = ", ")))
  }

  # 提取数据
  exposures <- data_test[, exposure_cols, drop = FALSE]
  father_exp <- exposures[[1]]
  mother_exp <- exposures[[2]]
  offspring_exp <- exposures[[3]]

  # --- 2. 计算协方差矩阵 Sigma ---
  # 计算协方差/方差
  sa2 <- (cov(father_exp, offspring_exp) + cov(mother_exp, offspring_exp)) / 2
  ss2 <- cov(father_exp, mother_exp)
  se2_F <- var(father_exp)
  se2_M <- var(mother_exp)
  se2_O <- var(offspring_exp)

  Sigma <- matrix(c(
    se2_F, ss2,   sa2,
    ss2,   se2_M, sa2,
    sa2,   sa2,   se2_O
  ), nrow = 3, byrow = TRUE)

  # 计算 Sigma 的逆
  Omega_inv <- tryCatch(
    {
      chol2inv(chol(Sigma))
    },
    error = function(e) {
      warning("协方差矩阵 Sigma 求逆失败 (可能奇异): ", e$message)
      matrix(NaN, 3, 3)
    }
  )
  if (any(is.nan(Omega_inv))) {
    return(list(beta = matrix(NaN, 1, 2), informatrix = matrix(NaN, 2, 2), beta_se = NaN))
  }

  # --- 3. 准备基因型数据 ---
  # 提取指定SNP位点的数据
  snps <- data_test[, snp_cols, drop = FALSE]
  father_snps <- snps[[1]]
  mother_snps <- snps[[2]]
  offspring_snps <- snps[[3]]

  # 计算等位基因频率 p
  # 使用母亲和后代的SNP数据防止人群分层
  p <- mean(c(mother_snps, offspring_snps)) / 2

  # 计算 G1, G2, G3 (向量化)
  G1 <- ifelse(father_snps == 0, 2 * p,
    ifelse(father_snps == 1, 1 + 2 * p, 2 * (1 + p))
  )
  G2 <- ifelse(mother_snps == 0, 2 * p,
    ifelse(mother_snps == 1, 1 + 2 * p, 2 * (1 + p))
  )
  G3 <- father_snps + mother_snps

  # 构建矩阵
  X <- as.matrix(snps) # n x 3
  X_par <- cbind(G1, G2, G3) # n x 3
  Y <- as.matrix(exposures) # n x 3

  # --- 4. 向量化计算 fenzi 和 informatrix ---
  # 预计算常用乘积
  Y_OmInv <- Y %*% Omega_inv # n x 3
  X_OmInv <- X %*% Omega_inv # n x 3
  Xp_OmInv <- X_par %*% Omega_inv # n x 3

  # 计算 fenzi (1x2 矩阵)
  fenzi_1 <- sum(Y_OmInv * X)
  fenzi_2 <- sum(Y_OmInv * X_par)
  fenzi <- matrix(c(fenzi_1, fenzi_2), nrow = 1)

  # 计算 informatrix (2x2 矩阵)
  info_11 <- sum(X_OmInv * X)
  info_12 <- sum(X_OmInv * X_par)
  info_22 <- sum(Xp_OmInv * X_par)

  informatrix <- matrix(c(
    info_11, info_12,
    info_12, info_22
  ), nrow = 2, byrow = TRUE)

  # --- 5. 计算 Beta 和标准误 ---
  # 计算 informatrix 的逆
  informatrix_inv <- tryCatch(
    {
      chol2inv(chol(informatrix))
    },
    error = function(e) {
      warning("信息矩阵 informatrix 求逆失败 (可能奇异): ", e$message)
      matrix(NaN, 2, 2)
    }
  )

  if (any(is.nan(informatrix_inv))) {
    beta <- matrix(NaN, 1, 2)
    beta_se1 <- NaN
  } else {
    # 计算 beta (1x2 矩阵)
    beta <- fenzi %*% informatrix_inv

    # 计算第一个 beta 的标准误
    beta_se1 <- sqrt(pmax(informatrix_inv[1, 1], 0))
  }

  # --- 6. 返回结果 ---
  result <- list(
    beta = beta, # 1x2 矩阵
    informatrix = informatrix, # 2x2 信息矩阵
    beta_se = beta_se1, # 第一个 beta 的标准误 (标量)
    snp_used = snp_index, # 记录使用的SNP位点
    allele_freq = p # 记录等位基因频率
  )

  return(result)
}

# 批量分析多个SNP位点的辅助函数 - 合并输出版本
fmr_trio_multi_snp <- function(data_test, snp_indices = NULL) {
  # 如果未指定SNP位点，自动检测所有可用的SNP位点
  if (is.null(snp_indices)) {
    # 查找所有father_snps列
    father_snp_cols <- grep("^father_snps\\.", names(data_test), value = TRUE)
    snp_indices <- as.numeric(gsub("father_snps\\.", "", father_snp_cols))
    snp_indices <- sort(snp_indices[!is.na(snp_indices)])
  }

  n_snps <- length(snp_indices)

  # 初始化合并矩阵
  combined_beta <- matrix(NaN, nrow = n_snps, ncol = 2)
  combined_informatrix <- matrix(NaN, nrow = n_snps * 2, ncol = 2)
  combined_beta_se <- rep(NaN, n_snps)
  combined_allele_freq <- rep(NaN, n_snps)

  # 设置行名
  rownames(combined_beta) <- paste0("SNP_", snp_indices)
  colnames(combined_beta) <- c("Beta1", "Beta2")

  informatrix_rownames <- c()
  for (i in snp_indices) {
    informatrix_rownames <- c(
      informatrix_rownames,
      paste0("SNP_", i, "_Info1"),
      paste0("SNP_", i, "_Info2")
    )
  }
  rownames(combined_informatrix) <- informatrix_rownames
  colnames(combined_informatrix) <- c("Col1", "Col2")

  names(combined_beta_se) <- paste0("SNP_", snp_indices)
  names(combined_allele_freq) <- paste0("SNP_", snp_indices)

  # 记录成功分析的SNP
  successful_snps <- c()
  failed_snps <- c()

  # 为每个SNP位点运行分析
  for (idx in 1:n_snps) {
    i <- snp_indices[idx]
    cat("分析SNP位点", i, "...\n")

    tryCatch(
      {
        result <- fmr_trio_optimized(data_test, snp_index = i)

        # 检查结果是否有效
        if (!any(is.nan(result$beta))) {
          # 将结果填入合并矩阵
          combined_beta[idx, ] <- result$beta[1, ]

          # informatrix是2x2矩阵，需要按行填入
          start_row <- (idx - 1) * 2 + 1
          end_row <- idx * 2
          combined_informatrix[start_row:end_row, ] <- result$informatrix

          combined_beta_se[idx] <- result$beta_se
          combined_allele_freq[idx] <- result$allele_freq

          successful_snps <- c(successful_snps, i)
        } else {
          failed_snps <- c(failed_snps, i)
          cat("SNP位点", i, "分析结果包含NaN值\n")
        }
      },
      error = function(e) {
        cat("SNP位点", i, "分析失败:", e$message, "\n")
        failed_snps <<- c(failed_snps, i)
      }
    )
  }

  # 汇总信息
  cat("\n=== 分析汇总 ===\n")
  cat("总SNP数:", n_snps, "\n")
  cat("成功分析:", length(successful_snps), "个SNP\n")
  cat("失败分析:", length(failed_snps), "个SNP\n")
  if (length(failed_snps) > 0) {
    cat("失败的SNP:", paste(failed_snps, collapse = ", "), "\n")
  }

  temp <- combined_beta[, 1]
  combined_beta[, 1] <- combined_beta[, 2]
  combined_beta[, 2] <- temp
  # 返回合并结果
  result <- list(
    beta = combined_beta, # n_snps x 2 矩阵
    informatrix = combined_informatrix, # (n_snps*2) x 2 矩阵
    beta_se = combined_beta_se, # n_snps 向量
    allele_freq = combined_allele_freq, # n_snps 向量
    snp_indices = snp_indices, # 分析的SNP位点
    successful_snps = successful_snps, # 成功分析的SNP位点
    failed_snps = failed_snps, # 失败的SNP位点
    n_snps = n_snps # SNP总数
  )

  return(result)
}

# 提取有效结果的辅助函数
extract_valid_results <- function(combined_results) {
  # 找出非NaN的行
  valid_rows <- !is.nan(combined_results$beta[, 1])

  if (sum(valid_rows) == 0) {
    cat("警告: 没有有效的分析结果\n")
    return(NULL)
  }

  # 提取有效结果
  valid_beta <- combined_results$beta[valid_rows, , drop = FALSE]
  valid_beta_se <- combined_results$beta_se[valid_rows]
  valid_allele_freq <- combined_results$allele_freq[valid_rows]

  # 对于informatrix，需要提取对应的行
  valid_informatrix_rows <- c()
  for (i in which(valid_rows)) {
    valid_informatrix_rows <- c(
      valid_informatrix_rows,
      (i - 1) * 2 + 1,
      (i - 1) * 2 + 2
    )
  }
  valid_informatrix <- combined_results$informatrix[valid_informatrix_rows, , drop = FALSE]

  return(list(
    beta = valid_beta,
    informatrix = valid_informatrix,
    beta_se = valid_beta_se,
    allele_freq = valid_allele_freq,
    n_valid = sum(valid_rows)
  ))
}

# %% LMM 不用FGWAS
gls_trio_optimized <- function(data_test, snp_index = 1) {
  # 参数说明:
  # data_test: 输入数据框
  # snp_index: 指定使用哪个SNP位点进行分析 (默认为1，即第一个SNP)

  n <- nrow(data_test) # 获取样本量

  # --- 1. 数据预处理和列名匹配 ---
  # 处理exposure列名的大小写差异
  exposure_cols <- c()
  if ("father_expose" %in% names(data_test)) {
    exposure_cols <- c("father_expose", "mother_expose", "offspring_expose")
  } else if ("Father_expose" %in% names(data_test)) {
    exposure_cols <- c("Father_expose", "Mother_expose", "Offspring_expose")
  } else {
    stop("未找到exposure列，请检查列名")
  }

  # 构建SNP列名
  snp_cols <- c(
    paste0("father_snps.", snp_index),
    paste0("mother_snps.", snp_index),
    paste0("offspring_snps.", snp_index)
  )

  # 检查SNP列是否存在
  if (!all(snp_cols %in% names(data_test))) {
    stop(paste("未找到SNP列:", paste(snp_cols, collapse = ", ")))
  }

  # 提取数据
  exposures <- data_test[, exposure_cols, drop = FALSE]
  father_exp <- exposures[[1]]
  mother_exp <- exposures[[2]]
  offspring_exp <- exposures[[3]]

  # --- 2. 计算协方差矩阵 Sigma ---
  # 计算协方差/方差
  sa2 <- (cov(father_exp, offspring_exp) + cov(mother_exp, offspring_exp)) / 2
  ss2 <- cov(father_exp, mother_exp)
  se2_F <- var(father_exp)
  se2_M <- var(mother_exp)
  se2_O <- var(offspring_exp)

  Sigma <- matrix(c(
    se2_F, ss2,   sa2,
    ss2,   se2_M, sa2,
    sa2,   sa2,   se2_O
  ), nrow = 3, byrow = TRUE)

  # 计算 Sigma 的逆
  Omega_inv <- tryCatch(
    {
      chol2inv(chol(Sigma))
    },
    error = function(e) {
      warning("协方差矩阵 Sigma 求逆失败 (可能奇异): ", e$message)
      matrix(NaN, 3, 3)
    }
  )
  if (any(is.nan(Omega_inv))) {
    return(list(beta = matrix(NaN, 1, 1), informatrix = matrix(NaN, 1, 1), beta_se = NaN))
  }

  # --- 3. 准备基因型数据 ---
  # 提取指定SNP位点的数据
  snps <- data_test[, snp_cols, drop = FALSE]
  father_snps <- snps[[1]]
  mother_snps <- snps[[2]]
  offspring_snps <- snps[[3]]

  # 计算等位基因频率 p (用于记录，但不用于G1G2G3计算)
  # 使用母亲和后代的SNP数据防止人群分层
  p <- mean(c(mother_snps, offspring_snps)) / 2

  # 构建矩阵 - 直接使用原始SNP数据
  X <- as.matrix(snps) # n x 3 (father, mother, offspring SNPs)
  Y <- as.matrix(exposures) # n x 3 (father, mother, offspring exposures)

  # --- 4. 向量化计算 fenzi 和 informatrix ---
  # 预计算常用乘积
  Y_OmInv <- Y %*% Omega_inv # n x 3
  X_OmInv <- X %*% Omega_inv # n x 3

  # 计算 fenzi (1x1 矩阵 - 只有一个beta)
  fenzi <- sum(Y_OmInv * X)

  # 计算 informatrix (1x1 矩阵)
  informatrix <- sum(X_OmInv * X)

  # --- 5. 计算 Beta 和标准误 ---
  # 检查 informatrix 是否为零或接近零
  if (abs(informatrix) < .Machine$double.eps) {
    warning("信息矩阵接近零，无法计算beta")
    beta <- NaN
    beta_se <- NaN
  } else {
    # 计算 beta (标量)
    beta <- fenzi / informatrix

    # 计算标准误
    beta_se <- sqrt(1 / informatrix)
  }

  # --- 6. 返回结果 ---
  result <- list(
    beta = beta, # 标量
    informatrix = informatrix, # 标量
    beta_se = beta_se, # 标量
    snp_used = snp_index, # 记录使用的SNP位点
    allele_freq = p # 记录等位基因频率
  )

  return(result)
}

# 批量分析多个SNP位点的辅助函数 - 升级版
gls_trio_multi_snp <- function(data_test, snp_indices = NULL) {
  # 如果未指定SNP位点，自动检测所有可用的SNP位点
  if (is.null(snp_indices)) {
    # 查找所有father_snps列
    father_snp_cols <- grep("^father_snps\\.", names(data_test), value = TRUE)
    snp_indices <- as.numeric(gsub("father_snps\\.", "", father_snp_cols))
    snp_indices <- sort(snp_indices[!is.na(snp_indices)])
  }

  n_snps <- length(snp_indices)

  # 初始化结果向量 (现在beta是标量)
  combined_beta <- rep(NaN, n_snps)
  combined_informatrix <- rep(NaN, n_snps)
  combined_beta_se <- rep(NaN, n_snps)
  combined_allele_freq <- rep(NaN, n_snps)

  # 设置名称
  names(combined_beta) <- paste0("SNP_", snp_indices)
  names(combined_informatrix) <- paste0("SNP_", snp_indices)
  names(combined_beta_se) <- paste0("SNP_", snp_indices)
  names(combined_allele_freq) <- paste0("SNP_", snp_indices)

  # 记录成功分析的SNP
  successful_snps <- c()
  failed_snps <- c()

  # 为每个SNP位点运行分析
  for (idx in 1:n_snps) {
    i <- snp_indices[idx]
    cat("分析SNP位点", i, "...\n")

    tryCatch(
      {
        result <- gls_trio_optimized(data_test, snp_index = i)

        # 检查结果是否有效
        if (!is.nan(result$beta)) {
          # 将结果填入向量
          combined_beta[idx] <- result$beta
          combined_informatrix[idx] <- result$informatrix
          combined_beta_se[idx] <- result$beta_se
          combined_allele_freq[idx] <- result$allele_freq

          successful_snps <- c(successful_snps, i)
        } else {
          failed_snps <- c(failed_snps, i)
          cat("SNP位点", i, "分析结果包含NaN值\n")
        }
      },
      error = function(e) {
        cat("SNP位点", i, "分析失败:", e$message, "\n")
        failed_snps <<- c(failed_snps, i)
      }
    )
  }

  # 汇总信息
  cat("\n=== 分析汇总 ===\n")
  cat("总SNP数:", n_snps, "\n")
  cat("成功分析:", length(successful_snps), "个SNP\n")
  cat("失败分析:", length(failed_snps), "个SNP\n")
  if (length(failed_snps) > 0) {
    cat("失败的SNP:", paste(failed_snps, collapse = ", "), "\n")
  }


  # 返回合并结果
  result <- list(
    beta = combined_beta, # n_snps 向量
    informatrix = combined_informatrix, # n_snps 向量
    beta_se = combined_beta_se, # n_snps 向量
    allele_freq = combined_allele_freq, # n_snps 向量
    snp_indices = snp_indices, # 分析的SNP位点
    successful_snps = successful_snps, # 成功分析的SNP位点
    failed_snps = failed_snps, # 失败的SNP位点
    n_snps = n_snps # SNP总数
  )

  return(result)
}

# 提取有效结果的辅助函数 - 升级版
gls_extract_valid_results <- function(combined_results) {
  # 找出非NaN的结果
  valid_indices <- !is.nan(combined_results$beta)

  if (sum(valid_indices) == 0) {
    cat("警告: 没有有效的分析结果\n")
    return(NULL)
  }

  # 提取有效结果
  valid_beta <- combined_results$beta[valid_indices]
  valid_informatrix <- combined_results$informatrix[valid_indices]
  valid_beta_se <- combined_results$beta_se[valid_indices]
  valid_allele_freq <- combined_results$allele_freq[valid_indices]

  return(list(
    beta = valid_beta,
    informatrix = valid_informatrix,
    beta_se = valid_beta_se,
    allele_freq = valid_allele_freq,
    n_valid = sum(valid_indices),
    valid_snp_names = names(combined_results$beta)[valid_indices]
  ))
}

# 结果汇总和可视化函数
summarize_results <- function(combined_results) {
  valid_results <- extract_valid_results(combined_results)

  if (is.null(valid_results)) {
    return(NULL)
  }

  cat("\n=== 有效结果汇总 ===\n")
  cat("有效SNP数量:", valid_results$n_valid, "\n")
  cat("Beta统计:\n")
  cat("  均值:", round(mean(valid_results$beta), 6), "\n")
  cat("  中位数:", round(median(valid_results$beta), 6), "\n")
  cat("  标准差:", round(sd(valid_results$beta), 6), "\n")
  cat("  范围:", round(range(valid_results$beta), 6), "\n")

  cat("\n标准误统计:\n")
  cat("  均值:", round(mean(valid_results$beta_se), 6), "\n")
  cat("  中位数:", round(median(valid_results$beta_se), 6), "\n")
  cat("  范围:", round(range(valid_results$beta_se), 6), "\n")

  # 计算显著性
  z_scores <- valid_results$beta / valid_results$beta_se
  p_values <- 2 * pnorm(-abs(z_scores))
  significant <- p_values < 0.05

  cat("\n显著性结果 (p < 0.05):\n")
  cat("  显著SNP数量:", sum(significant), "/", valid_results$n_valid, "\n")

  if (sum(significant) > 0) {
    cat("  显著的SNP:\n")
    sig_results <- data.frame(
      SNP = valid_results$valid_snp_names[significant],
      Beta = round(valid_results$beta[significant], 6),
      SE = round(valid_results$beta_se[significant], 6),
      Z_score = round(z_scores[significant], 3),
      P_value = format(p_values[significant], scientific = TRUE, digits = 3)
    )
    print(sig_results)
  }

  return(list(
    valid_results = valid_results,
    z_scores = z_scores,
    p_values = p_values,
    significant_indices = which(significant),
    summary_stats = list(
      n_valid = valid_results$n_valid,
      n_significant = sum(significant),
      beta_mean = mean(valid_results$beta),
      beta_median = median(valid_results$beta),
      beta_sd = sd(valid_results$beta)
    )
  ))
}
# %% mr普通版
fgwas_to_mr <- function(results_of_fgwas) {
  results_of_fgwas_beta <- matrix(results_of_fgwas$beta_hat, ncol = 2)[, 1]
  beta_se <- rep(0, nrow(results_of_fgwas$beta_hat))
  beta_se_p <- rep(0, nrow(results_of_fgwas$beta_hat))
  for (i in 1:nrow(results_of_fgwas$beta_hat)) {
    current_sigma_i <- (2 * i - 1):(2 * i)
    current_sigma <- results_of_fgwas$Sigma_inv[current_sigma_i, ]
    current_sigma_inv <- solve(current_sigma)
    beta_se[i] <- sqrt(current_sigma_inv[1, 1])
    beta_se_p[i] <- sqrt(current_sigma_inv[2, 2])
  }
  return(list(results_of_fgwas_beta = results_of_fgwas_beta, beta_se = beta_se, beta_se_p = beta_se_p))
}
# a <- generate_multiple_datasets_v3()
# results_of_fgwas <- cML_hat_expose_optimized(a[[1]])
# fgwas_to_MR(results_of_fgwas)


# %% 新增函数：封装 FGWAS 分析过程 一类错误膨胀

perform_fgwas_analysis_lmm <- function(data_set, n_snps) {
  # 确保加载了必要的包
  suppressPackageStartupMessages({
    library(dplyr)
    library(nlme)
  })


  # 动态确定 n_snps
  # 假设 snps 列的命名模式是 "snps.1", "snps.2", ...


  if (n_snps == 0) {
    stop("未在 data_set[[1]] (data_independent_exp) 中找到 snps 列，无法确定 n_snps。")
  }

  beta_hat_exp <- matrix(0, nrow = n_snps, ncol = 2)
  beta_hat_out <- matrix(0, nrow = n_snps, ncol = 2)

  beta_sigma_exp <- matrix(0, nrow = 2 * n_snps, ncol = 2)
  beta_sigma_out <- matrix(0, nrow = 2 * n_snps, ncol = 2)
  p_values_exp <- matrix(0, nrow = n_snps, ncol = 2) # 新增：用于存储P值
  f_statistic_exp <- matrix(0, nrow = n_snps, ncol = 2) # 新增：用于存储F统计量

  for (i in 1:n_snps) {
    snps_matrix_indicator <- (2 * i - 1):(2 * i)
    snps_indicator <- paste0("snps.", i)

    # 处理暴露数据
    data_independent_exp_i <- data_set[[1]] %>%
      dplyr::select(
        snps = ends_with(snps_indicator),
        expose = ends_with("expose")
      ) %>%
      dplyr::mutate(
        family_id = 1:n(),
        family_role = 4
      )

    data_trio_exp_i <- data_set[[2]] %>%
      dplyr::select(ends_with(snps_indicator), ends_with("expose"))
    names(data_trio_exp_i)[1] <- "father_snps"
    names(data_trio_exp_i)[2] <- "mother_snps"
    names(data_trio_exp_i)[3] <- "offspring_snps"
    # 求一下次等位基因频率
    f <- (mean(data_independent_exp_i$snps) / 2 +
      mean(c(
        data_trio_exp_i$offspring_snps,
        data_trio_exp_i$mother_snps
      )) / 2) / 2

    # 进行填补
    data_independent_exp_i <- data_independent_exp_i %>%
      dplyr::mutate(g_snps = snps + 2 * f)

    # 提取3个人并分开
    data_trio_exp_i_father <- data_trio_exp_i %>%
      dplyr::select("father_snps", "father_expose")
    data_trio_exp_i_father$family_id <-
      (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) + nrow(data_trio_exp_i_father))
    data_trio_exp_i_mother <- data_trio_exp_i %>%
      dplyr::select("mother_snps", "mother_expose")
    data_trio_exp_i_mother$family_id <-
      (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) + nrow(data_trio_exp_i_father))
    data_trio_exp_i_offspring <- data_trio_exp_i %>%
      dplyr::select("offspring_snps", "offspring_expose")
    data_trio_exp_i_offspring$family_id <-
      (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) + nrow(data_trio_exp_i_father))
    # 给他们都变成长数据然后拼起来
    names(data_trio_exp_i_father) <- c("snps", "expose", "family_id")
    data_trio_exp_i_father$family_role <- 1
    names(data_trio_exp_i_mother) <- c("snps", "expose", "family_id")
    data_trio_exp_i_mother$family_role <- 2
    names(data_trio_exp_i_offspring) <- c("snps", "expose", "family_id")
    data_trio_exp_i_offspring$family_role <- 3
    # 给这个几个长数据生成，协变量
    data_trio_exp_i_offspring$g_snps <- data_trio_exp_i_mother$snps +
      data_trio_exp_i_father$snps
    data_trio_exp_i_father$g_snps <- data_trio_exp_i_father$snps + 2 * f
    data_trio_exp_i_mother$g_snps <- data_trio_exp_i_mother$snps + 2 * f
    # 合成起来
    data_trio_exp_i_long <- bind_rows(
      data_trio_exp_i_father,
      data_trio_exp_i_mother, data_trio_exp_i_offspring
    )
    data_trio_exp_i_long <- arrange(data_trio_exp_i_long, family_id)
    # 重命名然后合并两种数据
    data_i_long <- bind_rows(
      data_independent_exp_i %>% dplyr::select(
        "family_id",
        "family_role", "snps", "g_snps", "expose"
      ),
      data_trio_exp_i_long %>% dplyr::select(
        "family_id",
        "family_role", "snps", "g_snps", "expose"
      )
    )
    # 得到了这个长数据
    data_i_long <- arrange(data_i_long, family_id)
    data_i_long <- na.omit(data_i_long)

    lmm_results_exp <- perform_lmm_analysis_v2(data_i_long,
      response_var_name = "expose"
    )

    if (!is.null(lmm_results_exp$model)) {
      model_summary <- summary(lmm_results_exp$model)
      coeffs <- model_summary$coefficients

      beta_hat_exp[i, ] <- coeffs[c("snps", "g_snps"), "Estimate"]
      vcov_mat <- lmm_results_exp$vcov_matrix
      beta_sigma_exp[snps_matrix_indicator, ] <-
        solve(vcov_mat[c("snps", "g_snps"), c("snps", "g_snps")])
    }

    data_independent_out_i <- data_set[[3]] %>%
      dplyr::select(
        snps = ends_with(snps_indicator),
        outcome = ends_with("_outcome")
      ) %>%
      dplyr::mutate(
        family_id = 1:n(),
        family_role = 4
      )

    data_trio_out_i <- data_set[[4]] %>%
      dplyr::select(ends_with(snps_indicator), ends_with("_outcome"))
    names(data_trio_out_i)[1] <- "father_snps"
    names(data_trio_out_i)[2] <- "mother_snps"
    names(data_trio_out_i)[3] <- "offspring_snps"
    # 求一下次等位基因频率
    f_out <- (mean(data_independent_out_i$snps) / 2 +
      mean(c(
        data_trio_out_i$offspring_snps,
        data_trio_out_i$mother_snps
      )) / 2) / 2

    # 进行填补
    data_independent_out_i <- data_independent_out_i %>%
      dplyr::mutate(g_snps = snps + 2 * f_out)

    # 提取3个人并分开
    data_trio_out_i_father <- data_trio_out_i %>%
      dplyr::select("father_snps", "father_outcome")
    data_trio_out_i_father$family_id <-
      (nrow(data_independent_out_i) + 1):(nrow(data_independent_out_i) +
        nrow(data_trio_out_i_father))
    data_trio_out_i_mother <- data_trio_out_i %>%
      dplyr::select("mother_snps", "mother_outcome")
    data_trio_out_i_mother$family_id <-
      (nrow(data_independent_out_i) + 1):(nrow(data_independent_out_i) +
        nrow(data_trio_out_i_father))
    data_trio_out_i_offspring <- data_trio_out_i %>%
      dplyr::select("offspring_snps", "offspring_outcome")
    data_trio_out_i_offspring$family_id <-
      (nrow(data_independent_out_i) + 1):(nrow(data_independent_out_i) +
        nrow(data_trio_out_i_father))
    # 给他们都变成长数据然后拼起来
    names(data_trio_out_i_father) <- c("snps", "outcome", "family_id")
    data_trio_out_i_father$family_role <- 1
    names(data_trio_out_i_mother) <- c("snps", "outcome", "family_id")
    data_trio_out_i_mother$family_role <- 2
    names(data_trio_out_i_offspring) <- c("snps", "outcome", "family_id")
    data_trio_out_i_offspring$family_role <- 3
    # 给这个几个长数据生成，协变量
    data_trio_out_i_offspring$g_snps <- data_trio_out_i_mother$snps +
      data_trio_out_i_father$snps
    data_trio_out_i_father$g_snps <- data_trio_out_i_father$snps + 2 * f_out
    data_trio_out_i_mother$g_snps <- data_trio_out_i_mother$snps + 2 * f_out
    # 合成起来
    data_trio_out_i_long <- bind_rows(
      data_trio_out_i_father, data_trio_out_i_mother,
      data_trio_out_i_offspring
    )
    data_trio_out_i_long <- arrange(data_trio_out_i_long, desc(family_id))
    # 重命名然后合并两种数据
    data_i_long_out <- bind_rows(
      data_independent_out_i %>% dplyr::select(
        "family_id",
        "family_role", "snps", "g_snps", "outcome"
      ),
      data_trio_out_i_long %>% dplyr::select(
        "family_id",
        "family_role", "snps", "g_snps", "outcome"
      )
    )
    # 得到了这个长数据
    data_i_long_out <- arrange(data_i_long_out, family_id)
    data_i_long_out <- na.omit(data_i_long_out)
    lmm_results_out <- perform_lmm_analysis_v2(data_i_long_out,
      include_intercept = TRUE,
      response_var_name = "outcome"
    )
    if (!is.null(lmm_results_exp$model)) {
      model_summary <- summary(lmm_results_out$model)
      coeffs <- model_summary$coefficients

      beta_hat_out[i, ] <- coeffs[c("snps", "g_snps"), "Estimate"]

      vcov_mat <- lmm_results_out$vcov_matrix
      beta_sigma_out[snps_matrix_indicator, ] <-
        solve(vcov_mat[c("snps", "g_snps"), c("snps", "g_snps")])
    }
  }

  return(list(
    beta_hat_exp = beta_hat_exp,
    beta_hat_out = beta_hat_out,
    beta_sigma_exp = beta_sigma_exp,
    beta_sigma_out = beta_sigma_out,
    p_values_exp = p_values_exp, # 新增的P值矩阵
    f_statistic_exp = f_statistic_exp # 新增的F统计量矩阵
  ))
}

perform_lmm_analysis_v2 <- function(
    data, response_var_name = "expose",
    include_intercept = TRUE) {
  # 创建精确的“连接”分组因子
  data <- data %>%
    mutate(
      family_id_char = as.character(family_id),
      father_child_link_id = case_when(
        family_role %in% c(1, 3) ~ paste0("fc_", family_id_char),
        TRUE ~ NA_character_
      ),
      mother_child_link_id = case_when(
        family_role %in% c(2, 3) ~ paste0("mc_", family_id_char),
        TRUE ~ NA_character_
      )
    )


  # 构建模型公式
  fixed_effects <- "snps + g_snps"
  if (!include_intercept) {
    fixed_effects <- paste(fixed_effects, "- 1")
  }

  random_effects <- "(1 | family_id) + (1 | father_child_link_id) +
  (1 | mother_child_link_id)"
  formula_str <- paste(response_var_name, "~", fixed_effects, "+", random_effects)
  model_formula <- as.formula(formula_str)

  # 初始化返回值
  lmm_model <- NULL
  vcov_matrix <- matrix(NA,
    nrow = 3, ncol = 3,
    dimnames = list(
      c("(Intercept)", "snps", "g_snps"),
      c("(Intercept)", "snps", "g_snps")
    )
  )

  tryCatch(
    {
      lmm_fit <- lme4::lmer(
        formula = model_formula,
        data = data,
        control = lmerControl(
          check.nobs.vs.nlev = "ignore",
          check.nobs.vs.nRE = "ignore",
          optimizer = "bobyqa"
        )
      )
      lmm_model <- lmm_fit
      if (!is.null(lmm_model)) {
        vcov_matrix <- as.matrix(vcov(lmm_model))
      }
    },
    error = function(e) {
      cat("LMM (v2) model failed for SNP", i, "with error:", conditionMessage(e), "\n")
    }
  )

  return(list(model = lmm_model, vcov_matrix = vcov_matrix))
}
# %% 鲁棒版本的FGWAS
perform_fgwas_analysis_lmm_updated <- function(data_set, n_snps) {
  # 确保加载了必要的包
  suppressPackageStartupMessages({
    library(dplyr)
    library(lme4)
  })

  if (n_snps == 0) {
    stop("未在 data_set[[1]] (data_independent_exp) 中找到 snps 列，无法确定 n_snps。")
  }

  # 使用列表来动态存储成功的结果
  results_list <- list()

  for (i in 1:n_snps) {
    snps_matrix_indicator <- (2 * i - 1):(2 * i)
    snps_indicator <- paste0("snps.", i)

    cat(paste0("正在处理 SNP #", i, "...\n")) # 增加进度提示

    # 1. 处理暴露（Exposure）数据
    #----------------------------------
    data_independent_exp_i <- data_set[[1]] %>%
      dplyr::select(snps = ends_with(snps_indicator), expose = ends_with("expose")) %>%
      dplyr::mutate(family_id = 1:n(), family_role = 4)

    data_trio_exp_i <- data_set[[2]] %>%
      dplyr::select(ends_with(snps_indicator), ends_with("expose"))
    names(data_trio_exp_i) <- c(
      "father_snps", "mother_snps", "offspring_snps",
      "father_expose", "mother_expose", "offspring_expose"
    )

    f <- (mean(data_independent_exp_i$snps, na.rm = TRUE) / 2 +
      mean(c(data_trio_exp_i$offspring_snps, data_trio_exp_i$mother_snps), na.rm = TRUE) / 2) / 2

    data_independent_exp_i <- data_independent_exp_i %>% dplyr::mutate(g_snps = snps + 2 * f)

    data_trio_exp_i_father <- data_trio_exp_i %>%
      dplyr::select("father_snps", "father_expose") %>%
      mutate(family_id = (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) + n()))
    data_trio_exp_i_mother <- data_trio_exp_i %>%
      dplyr::select("mother_snps", "mother_expose") %>%
      mutate(family_id = (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) + n()))
    data_trio_exp_i_offspring <- data_trio_exp_i %>%
      dplyr::select("offspring_snps", "offspring_expose") %>%
      mutate(family_id = (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) + n()))

    names(data_trio_exp_i_father) <- c("snps", "expose", "family_id")
    data_trio_exp_i_father$family_role <- 1
    names(data_trio_exp_i_mother) <- c("snps", "expose", "family_id")
    data_trio_exp_i_mother$family_role <- 2
    names(data_trio_exp_i_offspring) <- c("snps", "expose", "family_id")
    data_trio_exp_i_offspring$family_role <- 3

    data_trio_exp_i_offspring$g_snps <- data_trio_exp_i_mother$snps + data_trio_exp_i_father$snps
    data_trio_exp_i_father$g_snps <- data_trio_exp_i_father$snps + 2 * f
    data_trio_exp_i_mother$g_snps <- data_trio_exp_i_mother$snps + 2 * f

    data_trio_exp_i_long <- bind_rows(data_trio_exp_i_father, data_trio_exp_i_mother, data_trio_exp_i_offspring)

    data_i_long <- bind_rows(
      data_independent_exp_i %>% dplyr::select("family_id", "family_role", "snps", "g_snps", "expose"),
      data_trio_exp_i_long %>% dplyr::select("family_id", "family_role", "snps", "g_snps", "expose")
    ) %>%
      arrange(family_id) %>%
      na.omit()

    lmm_results_exp <- perform_lmm_analysis_v2_updated(data_i_long, response_var_name = "expose", snp_index = i)


    # 2. 处理结局（Outcome）数据
    #----------------------------------
    data_independent_out_i <- data_set[[3]] %>%
      dplyr::select(snps = ends_with(snps_indicator), outcome = ends_with("_outcome")) %>%
      dplyr::mutate(family_id = 1:n(), family_role = 4)

    data_trio_out_i <- data_set[[4]] %>%
      dplyr::select(ends_with(snps_indicator), ends_with("_outcome"))
    names(data_trio_out_i) <- c(
      "father_snps", "mother_snps", "offspring_snps",
      "father_outcome", "mother_outcome", "offspring_outcome"
    )

    f_out <- (mean(data_independent_out_i$snps, na.rm = TRUE) / 2 +
      mean(c(data_trio_out_i$offspring_snps, data_trio_out_i$mother_snps), na.rm = TRUE) / 2) / 2

    data_independent_out_i <- data_independent_out_i %>% dplyr::mutate(g_snps = snps + 2 * f_out)

    data_trio_out_i_father <- data_trio_out_i %>%
      dplyr::select("father_snps", "father_outcome") %>%
      mutate(family_id = (nrow(data_independent_out_i) + 1):(nrow(data_independent_out_i) + n()))
    data_trio_out_i_mother <- data_trio_out_i %>%
      dplyr::select("mother_snps", "mother_outcome") %>%
      mutate(family_id = (nrow(data_independent_out_i) + 1):(nrow(data_independent_out_i) + n()))
    data_trio_out_i_offspring <- data_trio_out_i %>%
      dplyr::select("offspring_snps", "offspring_outcome") %>%
      mutate(family_id = (nrow(data_independent_out_i) + 1):(nrow(data_independent_out_i) + n()))

    names(data_trio_out_i_father) <- c("snps", "outcome", "family_id")
    data_trio_out_i_father$family_role <- 1
    names(data_trio_out_i_mother) <- c("snps", "outcome", "family_id")
    data_trio_out_i_mother$family_role <- 2
    names(data_trio_out_i_offspring) <- c("snps", "outcome", "family_id")
    data_trio_out_i_offspring$family_role <- 3

    data_trio_out_i_offspring$g_snps <- data_trio_out_i_mother$snps + data_trio_out_i_father$snps
    data_trio_out_i_father$g_snps <- data_trio_out_i_father$snps + 2 * f_out
    data_trio_out_i_mother$g_snps <- data_trio_out_i_mother$snps + 2 * f_out

    data_trio_out_i_long <- bind_rows(data_trio_out_i_father, data_trio_out_i_mother, data_trio_out_i_offspring)

    data_i_long_out <- bind_rows(
      data_independent_out_i %>% dplyr::select("family_id", "family_role", "snps", "g_snps", "outcome"),
      data_trio_out_i_long %>% dplyr::select("family_id", "family_role", "snps", "g_snps", "outcome")
    ) %>%
      arrange(family_id) %>%
      na.omit()

    lmm_results_out <- perform_lmm_analysis_v2_updated(data_i_long_out, response_var_name = "outcome", include_intercept = TRUE, snp_index = i)

    # 3. 检查模型是否都成功拟合，如果成功则保存结果
    #----------------------------------------------------
    if (!is.null(lmm_results_exp$model) && !is.null(lmm_results_out$model)) {
      # 提取暴露模型结果
      summary_exp <- summary(lmm_results_exp$model)
      coeffs_exp <- summary_exp$coefficients
      beta_hat_exp_i <- coeffs_exp[c("snps", "g_snps"), "Estimate"]
      vcov_mat_exp <- lmm_results_exp$vcov_matrix
      beta_sigma_exp_i <- solve(vcov_mat_exp[c("snps", "g_snps"), c("snps", "g_snps")])

      # 提取结局模型结果
      summary_out <- summary(lmm_results_out$model)
      coeffs_out <- summary_out$coefficients
      beta_hat_out_i <- coeffs_out[c("snps", "g_snps"), "Estimate"]
      vcov_mat_out <- lmm_results_out$vcov_matrix
      beta_sigma_out_i <- solve(vcov_mat_out[c("snps", "g_snps"), c("snps", "g_snps")])

      # 将当前成功位点的所有结果存入列表
      current_result <- list(
        snp_index = i,
        beta_hat_exp = beta_hat_exp_i,
        beta_hat_out = beta_hat_out_i,
        beta_sigma_exp = beta_sigma_exp_i,
        beta_sigma_out = beta_sigma_out_i
      )

      results_list[[length(results_list) + 1]] <- current_result
      cat(paste0("SNP #", i, " 拟合成功并已保存。\n"))
    } else {
      cat(paste0("SNP #", i, " 拟合失败，已跳过。\n"))
    }
  }

  # 4. 循环结束后，将列表中的结果重组成最终格式
  #----------------------------------------------------
  if (length(results_list) == 0) {
    warning("所有SNP位点的模型拟合均失败，返回空列表。")
    return(list())
  }

  # 提取成功位点的索引
  successful_snps_index <- sapply(results_list, function(x) x$snp_index)

  # 使用do.call和rbind将列表中的向量/矩阵堆叠起来
  beta_hat_exp <- do.call(rbind, lapply(results_list, function(x) x$beta_hat_exp))
  beta_hat_out <- do.call(rbind, lapply(results_list, function(x) x$beta_hat_out))
  beta_sigma_exp <- do.call(rbind, lapply(results_list, function(x) x$beta_sigma_exp))
  beta_sigma_out <- do.call(rbind, lapply(results_list, function(x) x$beta_sigma_out))

  # 为结果矩阵的行命名，以便追溯
  rownames(beta_hat_exp) <- successful_snps_index
  rownames(beta_hat_out) <- successful_snps_index

  # 对于sigma矩阵，每2行对应一个SNP
  rownames(beta_sigma_exp) <- rep(successful_snps_index, each = 2)
  rownames(beta_sigma_out) <- rep(successful_snps_index, each = 2)

  return(list(
    successful_snps_index = successful_snps_index,
    beta_hat_exp = beta_hat_exp,
    beta_hat_out = beta_hat_out,
    beta_sigma_exp = beta_sigma_exp,
    beta_sigma_out = beta_sigma_out
  ))
}



perform_lmm_analysis_v2_updated <- function(
    data,
    response_var_name = "expose",
    include_intercept = TRUE,
    snp_index = "N/A") {
  # (数据准备部分与之前相同，保持不变)
  data <- data %>%
    mutate(
      family_id_char = as.character(family_id),
      father_child_link_id = case_when(
        family_role %in% c(1, 3) ~ paste0("fc_", family_id_char),
        TRUE ~ NA_character_
      ),
      mother_child_link_id = case_when(
        family_role %in% c(2, 3) ~ paste0("mc_", family_id_char),
        TRUE ~ NA_character_
      )
    )

  # (模型公式构建部分与之前相同，保持不变)
  fixed_effects <- "snps + g_snps"
  if (!include_intercept) {
    fixed_effects <- paste0(response_var_name, " ~ ", fixed_effects, " - 1")
    model_formula <- as.formula(fixed_effects)
    vcov_dimnames <- list(c("snps", "g_snps"), c("snps", "g_snps"))
    vcov_size <- 2
  } else {
    random_effects <- "(1 | family_id) + (1 | father_child_link_id) + (1 | mother_child_link_id)"
    formula_str <- paste(response_var_name, "~", fixed_effects, "+", random_effects)
    model_formula <- as.formula(formula_str)
    vcov_dimnames <- list(c("(Intercept)", "snps", "g_snps"), c("(Intercept)", "snps", "g_snps"))
    vcov_size <- 3
  }

  # 初始化返回值
  lmm_model <- NULL
  vcov_matrix <- matrix(NA, nrow = vcov_size, ncol = vcov_size, dimnames = vcov_dimnames)

  tryCatch(
    {
      # 1. 运行lmer模型
      lmm_fit <- lme4::lmer(
        formula = model_formula,
        data = data,
        control = lmerControl(
          check.nobs.vs.nlev = "ignore",
          check.nobs.vs.nRE = "ignore",
          optimizer = "bobyqa"
        )
      )

      # 2. **【新增】检查模型的健康状况（收敛性和奇异性）**
      # 从模型对象中提取收敛信息
      conv_msg <- lmm_fit@optinfo$conv$lme4$messages
      is_singular <- isSingular(lmm_fit)

      # 如果存在任何收敛警告信息或模型是奇异的，则判定为失败
      if (!is.null(conv_msg) || is_singular) {
        cat(
          "LMM (v2) model for SNP index", snp_index, "for response '", response_var_name,
          "' fit with issues and will be DISCARDED.\n"
        )
        if (!is.null(conv_msg)) {
          cat("  - Convergence Message:", conv_msg, "\n")
        }
        if (is_singular) {
          cat("  - Singularity Issue: Model is singular.\n")
        }
        # 关键：不将 lmm_fit 赋值给 lmm_model，使得 lmm_model 保持为 NULL
      } else {
        # 3. 只有完全“健康”的模型才被接受
        lmm_model <- lmm_fit
        vcov_matrix <- as.matrix(vcov(lmm_model))
      }
    },
    error = function(e) {
      # 这个error处理器保持不变，用于捕捉更严重的、导致程序中断的错误
      cat(
        "LMM (v2) model FAILED (error) for SNP index", snp_index,
        "for response '", response_var_name, "' with error: ", conditionMessage(e), "\n"
      )
    }
  )

  return(list(model = lmm_model, vcov_matrix = vcov_matrix))
}

# %% 没有用FGWAS而是普通的填补
perform_normal_analysis_lmm <- function(data_set, n_snps) {
  # 确保加载了必要的包
  suppressPackageStartupMessages({
    library(dplyr)
    library(nlme)
  })

  if (n_snps == 0) {
    stop("未在 data_set[[1]] (data_independent_exp) 中找到 snps 列，无法确定 n_snps。")
  }

  # --- 修改：调整结果矩阵/向量的维度，因为现在只有一个固定效应 "snps" ---
  beta_hat_exp <- matrix(0, nrow = n_snps, ncol = 1)
  beta_hat_out <- matrix(0, nrow = n_snps, ncol = 1)
  # --- 修改：重命名并调整为向量，用于存储 snps 系数的方差 ---
  beta_var_exp <- numeric(n_snps)
  beta_var_out <- numeric(n_snps)

  colnames(beta_hat_exp) <- "snps"
  colnames(beta_hat_out) <- "snps"


  for (i in 1:n_snps) {
    snps_indicator <- paste0("snps.", i)

    # --- 处理暴露数据 ---
    data_independent_exp_i <- data_set[[1]] %>%
      dplyr::select(
        snps = ends_with(snps_indicator),
        expose = ends_with("expose")
      ) %>%
      dplyr::mutate(
        family_id = 1:n(),
        family_role = 4
      )

    data_trio_exp_i <- data_set[[2]] %>%
      dplyr::select(ends_with(snps_indicator), ends_with("expose"))
    # 为了代码清晰，重命名列
    names(data_trio_exp_i) <- c(
      "father_snps", "mother_snps", "offspring_snps",
      "father_expose", "mother_expose", "offspring_expose"
    )

    # 提取3个人并分开
    base_family_id <- nrow(data_independent_exp_i)

    data_trio_exp_i_father <- data_trio_exp_i %>%
      dplyr::select("father_snps", "father_expose") %>%
      dplyr::mutate(family_id = (base_family_id + 1):(base_family_id + n()), family_role = 1)

    data_trio_exp_i_mother <- data_trio_exp_i %>%
      dplyr::select("mother_snps", "mother_expose") %>%
      dplyr::mutate(family_id = (base_family_id + 1):(base_family_id + n()), family_role = 2)

    data_trio_exp_i_offspring <- data_trio_exp_i %>%
      dplyr::select("offspring_snps", "offspring_expose") %>%
      dplyr::mutate(family_id = (base_family_id + 1):(base_family_id + n()), family_role = 3)

    # 重命名列以匹配长数据格式
    names(data_trio_exp_i_father)[1:2] <- c("snps", "expose")
    names(data_trio_exp_i_mother)[1:2] <- c("snps", "expose")
    names(data_trio_exp_i_offspring)[1:2] <- c("snps", "expose")

    # --- 移除：所有关于 g_snps 的计算和处理 ---

    # 合并Trio数据
    data_trio_exp_i_long <- bind_rows(
      data_trio_exp_i_father,
      data_trio_exp_i_mother,
      data_trio_exp_i_offspring
    )

    # 合并独立样本和Trio数据
    data_i_long <- bind_rows(
      data_independent_exp_i,
      data_trio_exp_i_long
    ) %>%
      arrange(family_id, family_role) %>%
      na.omit()

    # 对暴露数据运行LMM
    lmm_results_exp <- perform_lmm_analysis_normal(data_i_long, response_var_name = "expose")

    if (!is.null(lmm_results_exp$model)) {
      model_summary <- summary(lmm_results_exp$model)
      coeffs <- coef(model_summary)

      # --- 修改：只提取 "snps" 的系数和方差 ---
      beta_hat_exp[i, 1] <- coeffs["snps", "Estimate"]
      vcov_mat <- lmm_results_exp$vcov_matrix
      beta_var_exp[i] <- vcov_mat["snps", "snps"] # 提取方差
    }

    # --- 处理结局数据 (逻辑与暴露数据相同) ---
    data_independent_out_i <- data_set[[3]] %>%
      dplyr::select(
        snps = ends_with(snps_indicator),
        outcome = ends_with("_outcome")
      ) %>%
      dplyr::mutate(
        family_id = 1:n(),
        family_role = 4
      )

    data_trio_out_i <- data_set[[4]] %>%
      dplyr::select(ends_with(snps_indicator), ends_with("_outcome"))
    names(data_trio_out_i) <- c(
      "father_snps", "mother_snps", "offspring_snps",
      "father_outcome", "mother_outcome", "offspring_outcome"
    )

    base_family_id_out <- nrow(data_independent_out_i)

    data_trio_out_i_father <- data_trio_out_i %>%
      dplyr::select("father_snps", "father_outcome") %>%
      dplyr::mutate(family_id = (base_family_id_out + 1):(base_family_id_out + n()), family_role = 1)

    data_trio_out_i_mother <- data_trio_out_i %>%
      dplyr::select("mother_snps", "mother_outcome") %>%
      dplyr::mutate(family_id = (base_family_id_out + 1):(base_family_id_out + n()), family_role = 2)

    data_trio_out_i_offspring <- data_trio_out_i %>%
      dplyr::select("offspring_snps", "offspring_outcome") %>%
      dplyr::mutate(family_id = (base_family_id_out + 1):(base_family_id_out + n()), family_role = 3)

    names(data_trio_out_i_father)[1:2] <- c("snps", "outcome")
    names(data_trio_out_i_mother)[1:2] <- c("snps", "outcome")
    names(data_trio_out_i_offspring)[1:2] <- c("snps", "outcome")

    data_trio_out_i_long <- bind_rows(
      data_trio_out_i_father, data_trio_out_i_mother,
      data_trio_out_i_offspring
    )

    data_i_long_out <- bind_rows(
      data_independent_out_i,
      data_trio_out_i_long
    ) %>%
      arrange(family_id, family_role) %>%
      na.omit()

    # 对结局数据运行LMM
    lmm_results_out <- perform_lmm_analysis_normal(data_i_long_out,
      response_var_name = "outcome",
      include_intercept = TRUE
    )

    if (!is.null(lmm_results_out$model)) {
      model_summary <- summary(lmm_results_out$model)
      coeffs <- coef(model_summary)

      # --- 修改：只提取 "snps" 的系数和方差 ---
      beta_hat_out[i, 1] <- coeffs["snps", "Estimate"]
      vcov_mat <- lmm_results_out$vcov_matrix
      beta_var_out[i] <- vcov_mat["snps", "snps"] # 提取方差
    }
  }

  # --- 修改：返回更新后的结果列表 ---
  return(list(
    beta_hat_exp = beta_hat_exp,
    beta_hat_out = beta_hat_out,
    beta_var_exp = beta_var_exp,
    beta_var_out = beta_var_out
  ))
}

perform_lmm_analysis_normal <- function(
    data,
    response_var_name = "expose",
    include_intercept = TRUE) {
  # 创建用于定义随机效应的“连接”ID
  data <- data %>%
    mutate(
      family_id_char = as.character(family_id),
      father_child_link_id = case_when(
        family_role %in% c(1, 3) ~ paste0("fc_", family_id_char),
        TRUE ~ NA_character_
      ),
      mother_child_link_id = case_when(
        family_role %in% c(2, 3) ~ paste0("mc_", family_id_char),
        TRUE ~ NA_character_
      )
    )

  # 构建模型公式
  # 固定效应只包含 snps
  fixed_effects <- "snps"
  if (!include_intercept) {
    # 如果不要截距，则在公式中加入 -1
    fixed_effects <- "snps - 1"
  }

  # 定义随机效应部分
  random_effects <- "(1 | family_id) + (1 | father_child_link_id) + (1 | mother_child_link_id)"

  # 组合成完整公式
  formula_str <- paste(response_var_name, "~", fixed_effects, "+", random_effects)
  model_formula <- as.formula(formula_str)

  # 初始化返回值
  lmm_model <- NULL
  vcov_matrix <- NULL

  tryCatch(
    {
      lmm_fit <- lme4::lmer(
        formula = model_formula,
        data = data,
        control = lmerControl(
          check.nobs.vs.nlev = "ignore",
          check.nobs.vs.nRE = "ignore",
          optimizer = "bobyqa"
        )
      )
      lmm_model <- lmm_fit
      if (!is.null(lmm_model)) {
        # 直接从拟合好的模型中获取方差-协方差矩阵
        vcov_matrix <- as.matrix(vcov(lmm_model))
      }
    },
    error = function(e) {
      # 使用 `i` 可能会引发错误，因为 `i` 在此函数作用域内不存在
      # 建议传递 SNP 的标识符或者在外部处理错误信息
      cat("LMM model failed with error:", conditionMessage(e), "\n")
    }
  )

  return(list(model = lmm_model, vcov_matrix = vcov_matrix))
}
# %% 鲁棒性版本的普通lmm
perform_normal_analysis_lmm_updated <- function(data_set, n_snps) {
  # 确保加载了必要的包
  suppressPackageStartupMessages({
    library(dplyr)
    library(lme4) # lmer函数需要lme4包
  })

  if (n_snps == 0) {
    stop("未在 data_set[[1]] 中找到 snps 列，无法确定 n_snps。")
  }

  # 1. 初始化一个空列表来存储成功的结果
  results_list <- list()

  for (i in 1:n_snps) {
    snps_indicator <- paste0("snps.", i)
    cat(paste0("正在处理常规LMM SNP #", i, "...\n"))

    # --- 2. 处理暴露数据 ---
    data_independent_exp_i <- data_set[[1]] %>%
      dplyr::select(snps = ends_with(snps_indicator), expose = ends_with("expose")) %>%
      dplyr::mutate(family_id = 1:n(), family_role = 4)

    data_trio_exp_i <- data_set[[2]] %>%
      dplyr::select(ends_with(snps_indicator), ends_with("expose"))
    names(data_trio_exp_i) <- c(
      "father_snps", "mother_snps", "offspring_snps",
      "father_expose", "mother_expose", "offspring_expose"
    )

    base_family_id <- nrow(data_independent_exp_i)
    data_trio_exp_i_father <- data_trio_exp_i %>%
      dplyr::select("father_snps", "father_expose") %>%
      dplyr::mutate(family_id = (base_family_id + 1):(base_family_id + n()), family_role = 1, snps = father_snps, expose = father_expose)
    data_trio_exp_i_mother <- data_trio_exp_i %>%
      dplyr::select("mother_snps", "mother_expose") %>%
      dplyr::mutate(family_id = (base_family_id + 1):(base_family_id + n()), family_role = 2, snps = mother_snps, expose = mother_expose)
    data_trio_exp_i_offspring <- data_trio_exp_i %>%
      dplyr::select("offspring_snps", "offspring_expose") %>%
      dplyr::mutate(family_id = (base_family_id + 1):(base_family_id + n()), family_role = 3, snps = offspring_snps, expose = offspring_expose)

    data_trio_exp_i_long <- bind_rows(data_trio_exp_i_father, data_trio_exp_i_mother, data_trio_exp_i_offspring) %>%
      dplyr::select(family_id, family_role, snps, expose)

    data_i_long <- bind_rows(data_independent_exp_i, data_trio_exp_i_long) %>%
      arrange(family_id, family_role) %>%
      na.omit()

    # 调用升级版的LMM辅助函数
    lmm_results_exp <- perform_lmm_analysis_normal_updated(
      data_i_long,
      response_var_name = "expose",
      snp_index = i
    )

    # --- 3. 处理结局数据 (逻辑与暴露数据相同) ---
    data_independent_out_i <- data_set[[3]] %>%
      dplyr::select(snps = ends_with(snps_indicator), outcome = ends_with("_outcome")) %>%
      dplyr::mutate(family_id = 1:n(), family_role = 4)

    data_trio_out_i <- data_set[[4]] %>%
      dplyr::select(ends_with(snps_indicator), ends_with("_outcome"))
    names(data_trio_out_i) <- c(
      "father_snps", "mother_snps", "offspring_snps",
      "father_outcome", "mother_outcome", "offspring_outcome"
    )

    base_family_id_out <- nrow(data_independent_out_i)
    data_trio_out_i_father <- data_trio_out_i %>%
      dplyr::select("father_snps", "father_outcome") %>%
      dplyr::mutate(family_id = (base_family_id_out + 1):(base_family_id_out + n()), family_role = 1, snps = father_snps, outcome = father_outcome)
    data_trio_out_i_mother <- data_trio_out_i %>%
      dplyr::select("mother_snps", "mother_outcome") %>%
      dplyr::mutate(family_id = (base_family_id_out + 1):(base_family_id_out + n()), family_role = 2, snps = mother_snps, outcome = mother_outcome)
    data_trio_out_i_offspring <- data_trio_out_i %>%
      dplyr::select("offspring_snps", "offspring_outcome") %>%
      dplyr::mutate(family_id = (base_family_id_out + 1):(base_family_id_out + n()), family_role = 3, snps = offspring_snps, outcome = offspring_outcome)

    data_trio_out_i_long <- bind_rows(data_trio_out_i_father, data_trio_out_i_mother, data_trio_out_i_offspring) %>%
      dplyr::select(family_id, family_role, snps, outcome)

    data_i_long_out <- bind_rows(data_independent_out_i, data_trio_out_i_long) %>%
      arrange(family_id, family_role) %>%
      na.omit()

    lmm_results_out <- perform_lmm_analysis_normal_updated(
      data_i_long_out,
      response_var_name = "outcome",
      include_intercept = TRUE,
      snp_index = i
    )

    # 4. 检查两个模型是否都成功拟合，如果成功则保存结果
    if (!is.null(lmm_results_exp$model) && !is.null(lmm_results_out$model)) {
      # 提取暴露模型结果
      coeffs_exp <- coef(summary(lmm_results_exp$model))
      beta_hat_exp_i <- coeffs_exp["snps", "Estimate"]
      beta_var_exp_i <- vcov(lmm_results_exp$model)["snps", "snps"]

      # 提取结局模型结果
      coeffs_out <- coef(summary(lmm_results_out$model))
      beta_hat_out_i <- coeffs_out["snps", "Estimate"]
      beta_var_out_i <- vcov(lmm_results_out$model)["snps", "snps"]

      # 将当前成功位点的所有结果存入列表
      current_result <- list(
        snp_index = i,
        beta_hat_exp = beta_hat_exp_i,
        beta_var_exp = beta_var_exp_i,
        beta_hat_out = beta_hat_out_i,
        beta_var_out = beta_var_out_i
      )

      results_list[[length(results_list) + 1]] <- current_result
      cat(paste0("SNP #", i, " 常规LMM拟合成功并已保存。\n"))
    } else {
      cat(paste0("SNP #", i, " 常规LMM拟合失败或模型不可靠，已跳过。\n"))
    }
  }

  # 5. 循环结束后，将列表中的结果重组成最终格式
  if (length(results_list) == 0) {
    warning("所有SNP位点的常规LMM模型拟合均失败或不可靠，返回空列表。")
    return(list())
  }

  successful_snps_index <- sapply(results_list, function(x) x$snp_index)

  beta_hat_exp <- matrix(sapply(results_list, function(x) x$beta_hat_exp), ncol = 1)
  beta_hat_out <- matrix(sapply(results_list, function(x) x$beta_hat_out), ncol = 1)
  beta_var_exp <- sapply(results_list, function(x) x$beta_var_exp)
  beta_var_out <- sapply(results_list, function(x) x$beta_var_out)

  # 为结果命名，方便追溯
  rownames(beta_hat_exp) <- successful_snps_index
  colnames(beta_hat_exp) <- "snps"
  rownames(beta_hat_out) <- successful_snps_index
  colnames(beta_hat_out) <- "snps"
  names(beta_var_exp) <- successful_snps_index
  names(beta_var_out) <- successful_snps_index

  return(list(
    successful_snps_index = successful_snps_index,
    beta_hat_exp = beta_hat_exp,
    beta_hat_out = beta_hat_out,
    beta_var_exp = beta_var_exp,
    beta_var_out = beta_var_out
  ))
}

perform_lmm_analysis_normal_updated <- function(
    data,
    response_var_name = "expose",
    include_intercept = TRUE,
    snp_index = "N/A" # 新增参数用于日志
    ) {
  data <- data %>%
    mutate(
      family_id_char = as.character(family_id),
      father_child_link_id = case_when(
        family_role %in% c(1, 3) ~ paste0("fc_", family_id_char),
        TRUE ~ NA_character_
      ),
      mother_child_link_id = case_when(
        family_role %in% c(2, 3) ~ paste0("mc_", family_id_char),
        TRUE ~ NA_character_
      )
    )

  fixed_effects <- "snps"
  if (!include_intercept) {
    fixed_effects <- "snps - 1"
  }
  random_effects <- "(1 | family_id) + (1 | father_child_link_id) + (1 | mother_child_link_id)"
  formula_str <- paste(response_var_name, "~", fixed_effects, "+", random_effects)
  model_formula <- as.formula(formula_str)

  lmm_model <- NULL
  vcov_matrix <- NULL

  tryCatch(
    {
      lmm_fit <- lme4::lmer(
        formula = model_formula,
        data = data,
        control = lmerControl(
          check.nobs.vs.nlev = "ignore",
          check.nobs.vs.nRE = "ignore",
          optimizer = "bobyqa"
        )
      )

      # **健康检查**
      conv_msg <- lmm_fit@optinfo$conv$lme4$messages
      is_singular <- isSingular(lmm_fit)

      if (!is.null(conv_msg) || is_singular) {
        cat(
          "常规LMM (normal) model for SNP index", snp_index, "for response '", response_var_name,
          "' fit with issues and will be DISCARDED.\n"
        )
        if (!is.null(conv_msg)) cat("  - Convergence Message:", conv_msg, "\n")
        if (is_singular) cat("  - Singularity Issue: Model is singular.\n")
      } else {
        # 只有“健康”的模型才被接受
        lmm_model <- lmm_fit
        vcov_matrix <- as.matrix(vcov(lmm_model))
      }
    },
    error = function(e) {
      cat(
        "常规LMM (normal) model FAILED (error) for SNP index", snp_index,
        "for response '", response_var_name, "' with error: ", conditionMessage(e), "\n"
      )
    }
  )

  return(list(model = lmm_model, vcov_matrix = vcov_matrix))
}

# %% 测试一下这个gls是不是对的
