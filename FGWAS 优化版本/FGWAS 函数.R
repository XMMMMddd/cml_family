# %% 需要用到的包
library(lme4)
library(dplyr)
library(MASS)
library(nlme) # <-- 已添加
library(parallel)
library(pbapply)
# %% FGWAS 多个数据

# GLS 的 FGWAS 函数（对于三联体家庭用均值填补）
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

    # 检查子集是否为空
    if (nrow(subset_data) == 0) {
      warning(paste0("数据集 '", current_id, "' 为空，跳过处理。"))
      next # 跳到下一个数据集
    }

    # 调用处理函数 (例如 FMR_trio) **一次**
    # 使用 tryCatch 捕获该次调用的错误
    processed_result <- tryCatch(
      {
        processing_func(subset_data)
      },
      error = function(e) {
        warning(paste0("处理数据集 '", current_id, "' 时出错: ", e$message))
        return(NULL) # 返回 NULL 表示失败
      }
    )

    # --- 检查并存储结果 ---
    # 检查返回结果是否符合预期结构
    expected_structure <-
      !is.null(processed_result) &&
        is.list(processed_result) &&
        length(processed_result) >= 2 &&
        is.numeric(processed_result[[1]]) &&
        length(processed_result[[1]]) == 2 &&
        is.matrix(processed_result[[2]]) &&
        all(dim(processed_result[[2]]) == c(2, 2)) &&
        !anyNA(processed_result[[1]]) && # 检查 NA
        !anyNA(processed_result[[2]])

    if (expected_structure) {
      # 存储 beta 估计值
      beta_hat_matrix[i, ] <- processed_result[[1]]
      # 存储 2x2 的 Sigma_inv 块
      row_indices <- (2 * i - 1):(2 * i)
      Sigma_inv_matrix[row_indices, ] <- processed_result[[2]]
    } else {
      warning(paste0(
        "'processing_func' 未对数据集 '", current_id,
        "' 返回预期的结构或包含NA值。该数据集的结果将被设为 NA。"
      ))
      # beta_hat_matrix[i,] 和 Sigma_inv_matrix[row_indices,] 保持 NA
    }
  } # 结束 for 循环

  cat("数据集处理完成。\n")

  # --- 返回结果 ---
  # 检查是否有任何数据集处理成功
  successful_runs <- sum(!is.na(beta_hat_matrix[, 1])) # 检查第一列非NA的数量
  if (successful_runs == 0 && n_datasets > 0) {
    warning("所有数据集的处理均失败或返回无效结果。")
  } else if (successful_runs < n_datasets) {
    warning(paste("有", n_datasets - successful_runs, "个数据集处理失败或返回无效结果。"))
  }

  final_result <- list(beta_hat = beta_hat_matrix, Sigma_inv = Sigma_inv_matrix)
  return(final_result)
}


# 把 fgwas 的结果弄成好处理的
fgwas_to_mr <- function(results_of_fgwas) {
  results_of_fgwas_beta <- matrix(results_of_fgwas$beta_hat, ncol = 2)[, 1]
  beta_se <- rep(0, nrow(results_of_fgwas$beta_hat))
  for (i in 1:nrow(results_of_fgwas$beta_hat)) {
    current_sigma_i <- (2 * i - 1):(2 * i)
    current_sigma <- results_of_fgwas$sigma_inv[current_sigma_i, ]
    beta_se[i] <- sqrt(solve(current_sigma)[1, 1])
  }
  return(list(results_of_fgwas_beta = results_of_fgwas_beta, beta_se = beta_se))
}
# a <- generate_multiple_datasets_v3()
# results_of_fgwas <- cML_hat_expose_optimized(a[[1]])
# fgwas_to_MR(results_of_fgwas)


# %% 新增函数：封装 FGWAS 分析过程 一类错误膨胀

perform_fgwas_analysis <- function(data_set, n_snps) {
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

  for (i in 1:n_snps) {
    snps_matrix_indicator <- (2 * i - 1):(2 * i)
    snps_indicator <- paste0("snps.", i)
    data_independent_exp_i <- data_set[[1]] %>%
      dplyr::select(snps = ends_with(snps_indicator), expose = ends_with("expose")) %>%
      dplyr::mutate(
        family_id = 1:n(),
        family_role = "independent"
      )

    data_trio_exp_i <- data_set[[2]] %>%
      dplyr::select(ends_with(snps_indicator), ends_with("expose"))
    names(data_trio_exp_i)[1] <- "father_snps"
    names(data_trio_exp_i)[2] <- "mother_snps"
    names(data_trio_exp_i)[3] <- "offspring_snps"
    # 求一下次等位基因频率
    f <- (mean(data_independent_exp_i$snps) / 2 +
      mean(c(data_trio_exp_i$offspring_snps, data_trio_exp_i$mother_snps)) / 2) / 2

    # 进行填补
    data_independent_exp_i <- data_independent_exp_i %>%
      dplyr::mutate(g_snps = snps + 2 * f)

    # 提取3个人并分开
    data_trio_exp_i_father <- data_trio_exp_i %>% dplyr::select("father_snps", "father_expose")
    data_trio_exp_i_father$family_id <-
      (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) + nrow(data_trio_exp_i_father))
    data_trio_exp_i_mother <- data_trio_exp_i %>% dplyr::select("mother_snps", "mother_expose")
    data_trio_exp_i_mother$family_id <-
      (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) + nrow(data_trio_exp_i_father))
    data_trio_exp_i_offspring <- data_trio_exp_i %>% dplyr::select("offspring_snps", "offspring_expose")
    data_trio_exp_i_offspring$family_id <-
      (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) + nrow(data_trio_exp_i_father))
    # 给他们都变成长数据然后拼起来
    names(data_trio_exp_i_father) <- c("snps", "expose", "family_id")
    data_trio_exp_i_father$family_role <- "father"
    names(data_trio_exp_i_mother) <- c("snps", "expose", "family_id")
    data_trio_exp_i_mother$family_role <- "mother"
    names(data_trio_exp_i_offspring) <- c("snps", "expose", "family_id")
    data_trio_exp_i_offspring$family_role <- "offspring"
    # 给这个几个长数据生成，协变量
    data_trio_exp_i_offspring$g_snps <- data_trio_exp_i_mother$snps +
      data_trio_exp_i_father$snps
    data_trio_exp_i_father$g_snps <- data_trio_exp_i_father$snps + 2 * f
    data_trio_exp_i_mother$g_snps <- data_trio_exp_i_mother$snps + 2 * f
    # 合成起来
    data_trio_exp_i_long <- bind_rows(data_trio_exp_i_father, data_trio_exp_i_mother, data_trio_exp_i_offspring)
    data_trio_exp_i_long <- arrange(data_trio_exp_i_long, family_id)
    # 重命名然后合并两种数据
    data_i_long <- bind_rows(
      data_independent_exp_i %>% dplyr::select("family_id", "family_role", "snps", "g_snps", "expose"),
      data_trio_exp_i_long %>% dplyr::select("family_id", "family_role", "snps", "g_snps", "expose")
    )
    # 得到了这个长数据
    data_i_long <- arrange(data_i_long, family_id)
    data_lmm <- perform_lmm_analysis(data_i_long, include_intercept = TRUE)

    beta_hat_exp[i, ] <- fixef(data_lmm$model)[c(2, 3)]
    beta_sigma_exp[snps_matrix_indicator, ] <- data_lmm$information_matrix[2:3, 2:3]

    # 针对 data_set 中 out 的部分进行相同的处理
    data_independent_out_i <- data_set[[3]] %>%
      dplyr::select(snps = ends_with(snps_indicator), outcome = ends_with("_outcome")) %>%
      dplyr::mutate(
        family_id = 1:n(),
        family_role = "independent"
      )

    data_trio_out_i <- data_set[[4]] %>%
      dplyr::select(ends_with(snps_indicator), ends_with("_outcome"))
    names(data_trio_out_i)[1] <- "father_snps"
    names(data_trio_out_i)[2] <- "mother_snps"
    names(data_trio_out_i)[3] <- "offspring_snps"
    # 求一下次等位基因频率
    f_out <- (mean(data_independent_out_i$snps) / 2 +
      mean(c(data_trio_out_i$offspring_snps, data_trio_out_i$mother_snps)) / 2) / 2

    # 进行填补
    data_independent_out_i <- data_independent_out_i %>%
      dplyr::mutate(g_snps = snps + 2 * f_out)

    # 提取3个人并分开
    data_trio_out_i_father <- data_trio_out_i %>% dplyr::select("father_snps", "father_outcome")
    data_trio_out_i_father$family_id <-
      (nrow(data_independent_out_i) + 1):(nrow(data_independent_out_i) +
        nrow(data_trio_out_i_father))
    data_trio_out_i_mother <- data_trio_out_i %>%
      dplyr::select("mother_snps", "mother_outcome")
    data_trio_out_i_mother$family_id <-
      (nrow(data_independent_out_i) + 1):(nrow(data_independent_out_i) + nrow(data_trio_out_i_father))
    data_trio_out_i_offspring <- data_trio_out_i %>%
      dplyr::select("offspring_snps", "offspring_outcome")
    data_trio_out_i_offspring$family_id <-
      (nrow(data_independent_out_i) + 1):(nrow(data_independent_out_i) +
        nrow(data_trio_out_i_father))
    # 给他们都变成长数据然后拼起来
    names(data_trio_out_i_father) <- c("snps", "outcome", "family_id")
    data_trio_out_i_father$family_role <- "father"
    names(data_trio_out_i_mother) <- c("snps", "outcome", "family_id")
    data_trio_out_i_mother$family_role <- "mother"
    names(data_trio_out_i_offspring) <- c("snps", "outcome", "family_id")
    data_trio_out_i_offspring$family_role <- "offspring"
    # 给这个几个长数据生成，协变量
    data_trio_out_i_offspring$g_snps <- data_trio_out_i_mother$snps +
      data_trio_out_i_father$snps
    data_trio_out_i_father$g_snps <- data_trio_out_i_father$snps + 2 * f_out
    data_trio_out_i_mother$g_snps <- data_trio_out_i_mother$snps + 2 * f_out
    # 合成起来
    data_trio_out_i_long <- bind_rows(data_trio_out_i_father, data_trio_out_i_mother, data_trio_out_i_offspring)
    data_trio_out_i_long <- arrange(data_trio_out_i_long, desc(family_id))
    # 重命名然后合并两种数据
    data_i_long_out <- bind_rows(
      data_independent_out_i %>% dplyr::select("family_id", "family_role", "snps", "g_snps", "outcome"),
      data_trio_out_i_long %>% dplyr::select("family_id", "family_role", "snps", "g_snps", "outcome")
    )
    # 得到了这个长数据
    data_i_long_out <- arrange(data_i_long_out, family_id)
    data_lmm_out <- perform_lmm_analysis(data_i_long_out, include_intercept = TRUE, response_var_name = "outcome")

    beta_hat_out[i, ] <- fixef(data_lmm_out$model)[c(2, 3)]
    beta_sigma_out[snps_matrix_indicator, ] <-
      data_lmm_out$information_matrix[2:3, 2:3]
  }

  return(list(
    beta_hat_exp = beta_hat_exp,
    beta_hat_out = beta_hat_out,
    beta_sigma_exp = beta_sigma_exp,
    beta_sigma_out = beta_sigma_out
  ))
}
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

# %% 测试一下协方差结构
phase_two_data_full <- generate_mr_trio_data_ultra_updata(
  n_independent = 2000,
  p_trio = 0.9,
  beta_fs_oe_exp = 0.9,
  beta_ms_oe_exp = 0.9, beta_os_oe_exp = 0.9, beta_exp_to_out = 0
)
test_gls <- perform_fgwas_analysis_lmm(phase_two_data_full, n_snps = 3)
test_gls$beta_hat_out
test_gls$beta_hat_exp


phase_two_data_analysis_out <- list(
  beta_hat =
    test_gls$beta_hat_out,
  sigma_inv =
    test_gls$beta_sigma_out
)

phase_two_data_analysis_out <- fgwas_to_mr(
  phase_two_data_analysis_out
)
phase_two_data_analysis_out$p_value <-
  2 * pnorm(abs(
    phase_two_data_analysis_out$results_of_fgwas_beta /
      phase_two_data_analysis_out$beta_se
  ), lower.tail = FALSE)
phase_two_data_analysis_out$p_value

# %% 测试一下这个gls是不是对的


n_simulations <- 1000 # 您想运行的总次数
n_p_values <- 7 # 您要统计的p值的个数
set.seed(123) # 设置随机数种子以保证结果可复现

run_one_simulation <- function() {
  phase_two_data_full <- generate_mr_trio_data_ultra_updata(
    n_independent = 2000, p_trio = 0.9, beta_fs_oe_exp = 0.9,
    beta_ms_oe_exp = 0.9, beta_os_oe_exp = 0.9, beta_exp_to_out = 0
  )
  test_gls <- perform_fgwas_analysis_lmm(phase_two_data_full, n_snps = 3)
  phase_two_data_analysis_out <- list(
    beta_hat = test_gls$beta_hat_out, sigma_inv = test_gls$beta_sigma_out
  )
  phase_two_data_analysis_out <- fgwas_to_mr(phase_two_data_analysis_out)
  phase_two_data_analysis_out$p_value <-
    2 * pnorm(abs(
      phase_two_data_analysis_out$results_of_fgwas_beta /
        phase_two_data_analysis_out$beta_se
    ), lower.tail = FALSE)

  f_out <- mean(c(
    phase_two_data_full$data_trio_out$offspring_snps.1,
    phase_two_data_full$data_trio_out$mother_snps.1
  )) / 2
  g_snps_1 <- phase_two_data_full$data_trio_out$mother_snps.1 +
    phase_two_data_full$data_trio_out$father_snps.1
  g_snps_2 <- phase_two_data_full$data_trio_out$mother_snps.1 + 2 * f_out
  g_snps_3 <- phase_two_data_full$data_trio_out$father_snps.1 + 2 * f_out

  reference_p_value_1 <- summary(lm(
    phase_two_data_full$data_trio_out$offspring_outcome ~
      phase_two_data_full$data_trio_out$offspring_snps.1
  ))$coefficients[8]

  reference_p_value_2 <- summary(lm(
    phase_two_data_full$data_trio_out$mother_outcome ~
      phase_two_data_full$data_trio_out$mother_snps.1
  ))$coefficients[8]
  reference_p_value_3 <- summary(lm(
    phase_two_data_full$data_trio_out$father_outcome ~
      phase_two_data_full$data_trio_out$father_snps.1
  ))$coefficients[8]

  reference_p_value_5 <- summary(lm(
    phase_two_data_full$data_independent_out$offspring_outcome ~
      phase_two_data_full$data_independent_out$offspring_snps.1
  ))$coefficients[8]

  phase_two_data_analysis_out$p_value <- c(
    phase_two_data_analysis_out$p_value,
    reference_p_value_1, reference_p_value_2,
    reference_p_value_3, reference_p_value_5
  )
  return(phase_two_data_analysis_out$p_value < 0.05)
}

n_cores <- detectCores() - 1
if (n_cores < 1) n_cores <- 1
cat("使用", n_cores, "个核心进行并行计算...\n")

# 创建集群
cl <- makeCluster(n_cores)

# 定义需要在每个核心上加载的包和导出的函数
packages_to_load <- c("dplyr", "MASS", "nlme", "lme4") # <-- 已根据您的要求更新
functions_to_export <- c(
  "generate_mr_trio_data_ultra_updata",
  "perform_fgwas_analysis_lmm",
  "fgwas_to_mr",
  "run_one_simulation",
  "perform_lmm_analysis_v2"
)

tryCatch({
  # 第一步：将 'packages_to_load' 向量本身导出到每个工作核心
  clusterExport(cl, "packages_to_load")

  # 第二步：现在，每个核心都有了 packages_to_load 向量，可以执行加载命令了
  clusterEvalQ(cl, lapply(packages_to_load, library, character.only = TRUE))
  cat("依赖包 (dplyr, MASS, nlme) 已在所有核心上加载。\n")

  # 第三步：导出您自己的函数到每个工作进程
  clusterExport(cl, varlist = functions_to_export)
  cat("自定义函数已成功导出到所有核心。\n")

  # 第四步：使用 pblapply 进行并行计算
  cat("开始并行模拟...\n")
  results_list <- pblapply(1:n_simulations, function(i) {
    run_one_simulation()
  }, cl = cl)

  # 汇总结果
  significant_counts <- Reduce("+", results_list)
  names(significant_counts) <- paste0("p_value_", 1:n_p_values)

  # 打印最终结果
  cat("\n------ 模拟全部完成 ------\n")
  cat("总模拟次数:", n_simulations, "\n")
  cat("P值小于0.05的次数统计：\n")
  print(significant_counts)
}, error = function(e) {
  # 如果任何步骤出错，打印错误信息
  cat("\n并行计算过程中发生错误:\n")
  print(e)
}, finally = {
  # 【重要】: 无论成功与否，最后都要停止集群，释放资源
  stopCluster(cl)
  cat("\n并行集群已关闭。\n")
})

significant_counts / n_simulations

# %% 重构一下gls的代码

data_set <- generate_mr_trio_data_ultra_updata(
  n_independent = 2000, p_trio = 0.9, beta_fs_oe_exp = 0.9,
  beta_ms_oe_exp = 0.9, beta_os_oe_exp = 0.9, beta_exp_to_out = 0
)
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
  data_independent_exp_i <- data_set[[1]] %>%
    dplyr::select(snps = ends_with(snps_indicator), expose = ends_with("expose")) %>%
    dplyr::mutate(
      family_id = 1:n(),
      family_role = "independent"
    )

  data_trio_exp_i <- data_set[[2]] %>%
    dplyr::select(ends_with(snps_indicator), ends_with("expose"))
  names(data_trio_exp_i)[1] <- "father_snps"
  names(data_trio_exp_i)[2] <- "mother_snps"
  names(data_trio_exp_i)[3] <- "offspring_snps"
  str(data_trio_exp_i)
  # 求一下次等位基因频率
  f <- (mean(data_independent_exp_i$snps) / 2 +
    mean(c(data_trio_exp_i$offspring_snps, data_trio_exp_i$mother_snps)) / 2) / 2

  # 进行填补
  data_independent_exp_i <- data_independent_exp_i %>%
    dplyr::mutate(g_snps = snps + 2 * f)

  # 提取3个人并分开
  data_trio_exp_i_father <- data_trio_exp_i %>% dplyr::select("father_snps", "father_expose")
  data_trio_exp_i_father$family_id <-
    (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) +
      nrow(data_trio_exp_i_father))

  data_trio_exp_i_mother <- data_trio_exp_i %>% dplyr::select("mother_snps", "mother_expose")

  data_trio_exp_i_mother$family_id <-
    (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) +
      nrow(data_trio_exp_i_father))

  data_trio_exp_i_offspring <- data_trio_exp_i %>% dplyr::select(
    "offspring_snps",
    "offspring_expose"
  )

  data_trio_exp_i_offspring$family_id <-
    (nrow(data_independent_exp_i) + 1):(nrow(data_independent_exp_i) +
      nrow(data_trio_exp_i_father))


  # 给他们都变成长数据然后拼起来
  names(data_trio_exp_i_father) <- c("snps", "expose", "family_id")
  data_trio_exp_i_father$family_role <- "father"
  names(data_trio_exp_i_mother) <- c("snps", "expose", "family_id")
  data_trio_exp_i_mother$family_role <- "mother"
  names(data_trio_exp_i_offspring) <- c("snps", "expose", "family_id")
  data_trio_exp_i_offspring$family_role <- "offspring"
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
    data_independent_exp_i %>% dplyr::select("family_id", "family_role", "snps", "g_snps", "expose"),
    data_trio_exp_i_long %>% dplyr::select("family_id", "family_role", "snps", "g_snps", "expose")
  )
  # 得到了这个长数据
  data_i_long <- arrange(data_i_long, family_id)
  data_gls <- perform_gls_analysis(data_i_long, include_intercept = TRUE)

  beta_hat_exp[i, ] <- coef(data_gls$model)[c(2, 3)]
  beta_sigma_exp[snps_matrix_indicator, ] <-
    data_gls$information_matrix[2:3, 2:3]
  if (!is.null(data_gls$model)) {
    # 从模型摘要中提取P值
    model_summary <- summary(data_gls$model)
    p_values_exp[i, ] <- model_summary$tTable[c("snps", "g_snps"), "p-value"]

    # 从ANOVA表中提取F统计量
    model_anova <- anova(data_gls$model)
    f_statistic_exp[i, ] <- model_anova[c("snps", "g_snps"), "F-value"]
  }

  # 对汇总数据的进行
  data_independent_out_i <- data_set[[3]] %>%
    dplyr::select(snps = ends_with(snps_indicator), outcome = ends_with("_outcome")) %>%
    dplyr::mutate(
      family_id = 1:n(),
      family_role = "independent"
    )
  str(data_independent_out_i)
  data_trio_out_i <- data_set[[4]] %>%
    dplyr::select(ends_with(snps_indicator), ends_with("_outcome"))
  names(data_trio_out_i)[1] <- "father_snps"
  names(data_trio_out_i)[2] <- "mother_snps"
  names(data_trio_out_i)[3] <- "offspring_snps"
  str(data_trio_out_i)
  # 求一下次等位基因频率
  f_out <- (mean(data_independent_out_i$snps) / 2 +
    mean(c(data_trio_out_i$offspring_snps, data_trio_out_i$mother_snps)) / 2) / 2

  # 进行填补
  data_independent_out_i <- data_independent_out_i %>%
    dplyr::mutate(g_snps = snps + 2 * f_out)
  str(data_independent_out_i)
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
  str(data_trio_out_i_offspring)
  # 给他们都变成长数据然后拼起来
  names(data_trio_out_i_father) <- c("snps", "outcome", "family_id")
  data_trio_out_i_father$family_role <- "father"
  names(data_trio_out_i_mother) <- c("snps", "outcome", "family_id")
  data_trio_out_i_mother$family_role <- "mother"
  names(data_trio_out_i_offspring) <- c("snps", "outcome", "family_id")
  data_trio_out_i_offspring$family_role <- "offspring"
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
  data_i_long_out$family_id
  data_gls_out <- perform_gls_analysis(data_i_long_out,
    include_intercept = TRUE,
    response_var_name = "outcome"
  )

  beta_hat_out[i, ] <- coef(data_gls_out$model)[c(2, 3)]
  beta_sigma_out[snps_matrix_indicator, ] <-
    data_gls_out$information_matrix[2:3, 2:3]
}
