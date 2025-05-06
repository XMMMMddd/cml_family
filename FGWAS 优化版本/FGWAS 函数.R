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
  p <- mean(c(father_snps, mother_snps, offspring_snps)) / 2
  
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
da