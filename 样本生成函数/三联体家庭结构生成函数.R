# %% 载入要用到的包
library(dplyr)
library(MASS)

# %% 多数据集生成函数

generate_snp_hwe <- function(n_samples, maf) {
  # 根据次要等位基因频率(maf)计算基因型频率
  # AA (记为2), Aa (记为1), aa (记为0)
  # 假设'a'是次要等位基因, 'A'是主要等位基因
  # p = maf (a的频率), q = 1-maf (A的频率)
  # 基因型频率: aa (p^2), Aa (2pq), AA (q^2)
  genotype_freqs <- c(maf^2, 2 * maf * (1 - maf), (1 - maf)^2)
  genotypes <- sample(c(0, 1, 2), size = n_samples, replace = TRUE, prob = genotype_freqs)
  return(genotypes)
}
generate_mr_trio_data_sample <- function(
    n = 1000,
    p_f = 0.3, p_m = 0.3, # p_m 当前未在SNP生成中使用
    # 暴露效应
    beta_FStoOE_exp = 0.3, beta_MStoOE_exp = 0.3,
    beta_OStoOE_exp = 0.3,
    # 结局效应 (直接多效性 / 遗传叠加效应)
    beta_FStoOE_out = 0, beta_MStoOE_out = 0,
    beta_OStoOE_out = 0,
    # 因果效应
    beta_exp_to_out = 0.4,
    # 混杂效应
    confounding_exp = 0.2, confounding_out = 0.2,
    # 其他参数
    correlation = 0.2, seed = NULL,
    # 选型婚配强度
    assortative_mating_strength = 1000,
    sample.outcome.betas.from.range = FALSE) {
  # --- 0. 设置随机种子 (确保结果可重复性) ---
  set.seed(seed)



  #' @return 一个长度为N的数值向量，包含生成的基因型 (0, 1, 或 2)。
  internal_generate_hwe_snps <- function(N, p) {
    if (!is.numeric(p) || length(p) != 1 || p < 0 || p > 1) {
      stop("内部错误: internal_generate_hwe_snps 函数的 'p' 参数必须是0到1之间的单个数值。")
    }
    # 基因型频率: AA (记为2), Aa (记为1), aa (记为0)
    # p 是等位基因 A 的频率, (1-p) 是等位基因 a 的频率
    genotype_freqs <- c(p^2, 2 * p * (1 - p), (1 - p)^2) # 对应基因型 AA, Aa, aa
    # 使用 R 的 sample 函数生成基因型
    # 注意：提供的基因型值顺序与频率顺序对应 (2对应p^2, 1对应2pq, 0对应(1-p)^2)
    genotypes <- sample(c(2, 1, 0), size = N, replace = TRUE, prob = genotype_freqs)
    return(genotypes)
  }


  internal_get_transmitted_allele <- function(parent_snps) {
    transmitted_alleles <- numeric(length(parent_snps))
    for (i in seq_along(parent_snps)) {
      genotype <- parent_snps[i]
      if (genotype == 0) { # 亲本基因型 aa, 必定传递 a (编码为0)
        transmitted_alleles[i] <- 0
      } else if (genotype == 2) { # 亲本基因型 AA, 必定传递 A (编码为1)
        transmitted_alleles[i] <- 1
      } else if (genotype == 1) { # 亲本基因型 Aa, 传递 A 或 a 的概率各为50%
        transmitted_alleles[i] <- rbinom(1, 1, 0.5) # 随机传递0或1
      } else {
        warning(paste("在 internal_get_transmitted_allele 中发现无效的亲本基因型:", genotype, "位于索引", i, "- 将传递NA值"))
        transmitted_alleles[i] <- NA
      }
    }
    return(transmitted_alleles)
  }



  # 父方祖父母 (爷爷、奶奶)
  Grandfather_Father_SNPs <- generate_snp_hwe(n_samples = n, maf = p_f)
  Grandmother_Father_SNPs <- generate_snp_hwe(n_samples = n, maf = p_f)

  # 母方祖父母 (外公、外婆)
  Grandfather_Mother_SNPs <- generate_snp_hwe(n_samples = n, maf = p_f)
  Grandmother_Mother_SNPs <- generate_snp_hwe(n_samples = n, maf = p_f)


  # --- 2. 定义暴露和结局的生成函数 ---


  generate_expose_function <- function(beta_StoE_exp,
                                       beta_FStoOE_exp, beta_MStoOE_exp,
                                       SNPs_all,
                                       SNPs_Father,
                                       SNPs_Mother,
                                       confounder, correlation_factor) {
    # 计算暴露的确定性部分 (遗传效应 + 混杂 + 共享环境)
    results <- beta_StoE_exp * SNPs_all +
      beta_FStoOE_exp * SNPs_Father +
      beta_MStoOE_exp * SNPs_Mother +
      confounder + correlation_factor

    # 计算确定性部分的方差
    var_results <- var(results, na.rm = TRUE)

    # 计算随机误差的标准差，以使最终暴露的方差为1
    # 如果 var_results >= 1, 则误差标准差为0
    sd_err <- sqrt(max(0, 1 - var_results))

    # 添加随机误差，生成最终的暴露值
    expose_all <- results + rnorm(length(SNPs_all), mean = 0, sd = sd_err) # 注意：rnorm的N应与SNPs_all长度一致
    return(expose_all)
  }
  generate_outcome_function <- function(beta_exp_to_out_param,
                                        beta_OStoOE_out_param,
                                        beta_FStoOE_out_param,
                                        beta_MStoOE_out_param,
                                        expose_values, SNPs_self,
                                        SNPs_father_param, SNPs_mother_param,
                                        confounder_values, correlation_factor_values) {
    # 计算结局的确定性部分
    outcome_deterministic <- beta_exp_to_out_param * expose_values +
      beta_OStoOE_out_param * SNPs_self +
      beta_FStoOE_out_param * SNPs_father_param +
      beta_MStoOE_out_param * SNPs_mother_param +
      confounder_values + correlation_factor_values

    # 计算确定性部分的方差
    var_results <- var(outcome_deterministic, na.rm = TRUE)
    # 计算随机误差的标准差，以使最终结局的方差为1
    sd_err <- sqrt(max(0, 1 - var_results))
    # 添加随机误差，生成最终的结局值
    outcome_all <- outcome_deterministic + rnorm(length(SNPs_self),
      mean = 0, sd = sd_err
    ) # 注意：rnorm的N应与SNPs_self长度一致
    return(outcome_all)
  }

  # --- 3. 为所有 N_total 个体生成混杂因素和共享环境因素 ---
  confounder_exp_values <- rnorm(n,
    mean = 0,
    sd = sqrt(confounding_exp)
  )
  confounder_out_values <- rnorm(n,
    mean = 0,
    sd = sqrt(confounding_out)
  )
  correlation_factor_exp_values <- rnorm(n,
    mean = 0,
    sd = sqrt(correlation)
  )
  correlation_factor_out_values <- rnorm(n,
    mean = 0,
    sd = sqrt(correlation)
  )

  # --- 4. 为祖父母辈生成暴露数据 (用于选型婚配) ---
  Grandfather_Father_expose <- generate_expose_function(
    beta_OStoOE_exp, 0, 0,
    Grandfather_Father_SNPs, Grandfather_Father_SNPs, Grandfather_Father_SNPs,
    confounder_exp_values, correlation_factor_exp_values
  )
  Grandmother_Father_expose <- generate_expose_function(
    beta_OStoOE_exp, 0, 0,
    Grandmother_Father_SNPs, Grandmother_Father_SNPs, Grandmother_Father_SNPs,
    confounder_exp_values, correlation_factor_exp_values
  )
  Grandfather_Mother_expose <- generate_expose_function(
    beta_OStoOE_exp, 0, 0,
    Grandfather_Mother_SNPs, Grandfather_Mother_SNPs, Grandfather_Mother_SNPs,
    confounder_exp_values, correlation_factor_exp_values
  )
  Grandmother_Mother_expose <- generate_expose_function(
    beta_OStoOE_exp, 0, 0,
    Grandmother_Mother_SNPs, Grandmother_Mother_SNPs, Grandmother_Mother_SNPs,
    confounder_exp_values, correlation_factor_exp_values
  )


  Grandfather_Father_outcome <- generate_outcome_function(
    beta_exp_to_out, 0, 0, 0,
    Grandfather_Father_expose, Grandfather_Father_SNPs, Grandfather_Father_SNPs, Grandfather_Father_SNPs,
    confounder_out_values, correlation_factor_out_values
  )
  Grandmother_Father_outcome <- generate_outcome_function(
    beta_exp_to_out, 0, 0, 0,
    Grandmother_Father_expose, Grandmother_Father_SNPs, Grandmother_Father_SNPs, Grandmother_Father_SNPs,
    confounder_out_values, correlation_factor_out_values
  )
  Grandfather_Mother_outcome <- generate_outcome_function(
    beta_exp_to_out, 0, 0, 0,
    Grandfather_Mother_expose, Grandfather_Mother_SNPs, Grandfather_Mother_SNPs, Grandfather_Mother_SNPs,
    confounder_out_values, correlation_factor_out_values
  )
  Grandmother_Mother_outcome <- generate_outcome_function(
    beta_exp_to_out, 0, 0, 0,
    Grandmother_Mother_expose, Grandmother_Mother_SNPs, Grandmother_Mother_SNPs, Grandmother_Mother_SNPs,
    confounder_out_values, correlation_factor_out_values
  )


  # --- 5. 定义选型婚配函数 ---
  assortative_mating_function <- function(snps_1, expose_1, snps_2, expose_2,
                                          am_strength) {
    object_1 <- data.frame(snps = snps_1, expose = expose_1)
    object_2 <- data.frame(snps = snps_2, expose = expose_2)
    object_1 <- object_1 %>% dplyr::mutate(
      expose_star =
        expose + rnorm(dplyr::n(), mean = 0, sd = sqrt(am_strength))
    )
    object_2 <- object_2 %>% dplyr::mutate(
      expose_star =
        expose + rnorm(dplyr::n(), mean = 0, sd = sqrt(am_strength))
    )
    object_1 <- dplyr::arrange(object_1, expose_star)
    object_2 <- dplyr::arrange(object_2, expose_star)
    return(list(object_1, object_2))
  }

  # --- 6. 对祖父母辈进行选型婚配 ---
  grand_father_am_results <- assortative_mating_function(
    Grandfather_Father_SNPs, Grandfather_Father_expose,
    Grandmother_Father_SNPs, Grandmother_Father_outcome,
    assortative_mating_strength
  )
  grand_mother_am_results <- assortative_mating_function(
    Grandfather_Mother_SNPs, Grandfather_Mother_expose,
    Grandmother_Mother_SNPs, Grandmother_Mother_outcome,
    assortative_mating_strength
  )
  Grandfather_Father_SNPs <- grand_father_am_results[[1]]$snps
  Grandmother_Father_SNPs <- grand_father_am_results[[2]]$snps
  Grandfather_Mother_SNPs <- grand_mother_am_results[[1]]$snps
  Grandmother_Mother_SNPs <- grand_mother_am_results[[2]]$snps


  # --- 7. 生成父母代基因型 (通过模拟祖父母向父母传递等位基因) ---
  # 父亲的基因型 = 爷爷传递的等位基因 + 奶奶传递的等位基因
  Father_SNPs <- internal_get_transmitted_allele(Grandfather_Father_SNPs) +
    internal_get_transmitted_allele(Grandmother_Father_SNPs)
  # 母亲的基因型 = 外公传递的等位基因 + 外婆传递的等位基因
  Mother_SNPs <- internal_get_transmitted_allele(Grandfather_Mother_SNPs) +
    internal_get_transmitted_allele(Grandmother_Mother_SNPs)

  # --- 8. 为父母代生成暴露数据 (用于他们之间的选型婚配) ---
  # 父亲的暴露，受其自身SNP以及其父母(即子代的祖父母)SNP的遗传叠加效应影响
  Father_expose <- generate_expose_function(
    beta_OStoOE_exp, beta_FStoOE_exp, beta_MStoOE_exp, # 自身、其父(爷爷)、其母(奶奶)的SNP效应
    Father_SNPs, Grandfather_Father_SNPs, Grandmother_Father_SNPs,
    confounder_exp_values, correlation_factor_exp_values
  )
  # 母亲的暴露，受其自身SNP以及其父母(即子代的外祖父母)SNP的遗传叠加效应影响
  Mother_expose <- generate_expose_function(
    beta_OStoOE_exp, beta_FStoOE_exp, beta_MStoOE_exp, # 自身、其父(外公)、其母(外婆)的SNP效应
    Mother_SNPs, Grandfather_Mother_SNPs, Grandmother_Mother_SNPs, # 注意这里应该是祖父母的SNPs
    confounder_exp_values, correlation_factor_exp_values
  )

  # --- 12. 为父亲、母亲和子代生成结局数据 ---
  # 父亲的结局
  Father_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
    Father_expose, Father_SNPs, Grandfather_Father_SNPs, Grandmother_Father_SNPs, # 父亲的"父母"是爷爷奶奶
    confounder_out_values, correlation_factor_out_values
  )
  # 母亲的结局
  Mother_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
    Mother_expose, Mother_SNPs, Grandfather_Mother_SNPs, Grandmother_Mother_SNPs, # 母亲的"父母"是外公外婆
    confounder_out_values, correlation_factor_out_values
  )
  # --- 9. 对父母代进行选型婚配 ---
  parent_am_results <- assortative_mating_function(
    Father_SNPs, Father_expose, # 父亲的信息
    Mother_SNPs, Mother_outcome, # 母亲的信息
    assortative_mating_strength
  )
  # 更新父母的SNP和暴露向量，使其根据选型婚配的结果重新排序
  Father_SNPs <- parent_am_results[[1]]$snps # 排序后的父亲SNP
  Mother_SNPs <- parent_am_results[[2]]$snps # 排序后并与父亲配对的母亲SNP
  Mother_SNPs <- parent_am_results[[2]]$snps
  Father_expose <- generate_expose_function(
    beta_OStoOE_exp, beta_FStoOE_exp, beta_MStoOE_exp, # 自身、其父(爷爷)、其母(奶奶)的SNP效应
    Father_SNPs, Grandfather_Father_SNPs, Grandmother_Father_SNPs,
    confounder_exp_values, correlation_factor_exp_values
  )
  # 母亲的暴露，受其自身SNP以及其父母(即子代的外祖父母)SNP的遗传叠加效应影响
  Mother_expose <- generate_expose_function(
    beta_OStoOE_exp, beta_FStoOE_exp, beta_MStoOE_exp, # 自身、其父(外公)、其母(外婆)的SNP效应
    Mother_SNPs, Grandfather_Mother_SNPs, Grandmother_Mother_SNPs, # 注意这里应该是祖父母的SNPs
    confounder_exp_values, correlation_factor_exp_values
  )

  # --- 12. 为父亲、母亲和子代生成结局数据 ---
  # 父亲的结局
  Father_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
    Father_expose, Father_SNPs, Grandfather_Father_SNPs, Grandmother_Father_SNPs, # 父亲的"父母"是爷爷奶奶
    confounder_out_values, correlation_factor_out_values
  )
  # 母亲的结局
  Mother_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
    Mother_expose, Mother_SNPs, Grandfather_Mother_SNPs, Grandmother_Mother_SNPs, # 母亲的"父母"是外公外婆
    confounder_out_values, correlation_factor_out_values
  )
  # --- 9. 对父母代进行选型婚配 ---





  # --- 10. 生成子代基因型和暴露数据 ---
  # 子代的基因型 = 父亲传递的等位基因 + 母亲传递的等位基因
  Offspring_SNPs <- internal_get_transmitted_allele(Father_SNPs) +
    internal_get_transmitted_allele(Mother_SNPs)

  # 子代的暴露，受其自身SNP以及其父母SNP的遗传叠加效应影响
  Offspring_expose <- generate_expose_function(
    beta_OStoOE_exp, beta_FStoOE_exp, beta_MStoOE_exp, # 自身、父亲、母亲的SNP效应
    Offspring_SNPs, Father_SNPs, Mother_SNPs,
    confounder_exp_values, correlation_factor_exp_values
  )
  Offspring_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
    Offspring_expose, Offspring_SNPs, Father_SNPs, Mother_SNPs, # 子代的父母是Father和Mother
    confounder_out_values, correlation_factor_out_values
  )
  data <- data.frame(
    Father_SNPs = Father_SNPs,
    Mother_SNPs = Mother_SNPs,
    Offspring_SNPs = Offspring_SNPs,
    Father_expose = Father_expose,
    Mother_expose = Mother_expose,
    Offspring_expose = Offspring_expose,
    Father_outcome = Father_outcome,
    Mother_outcome = Mother_outcome,
    Offspring_outcome = Offspring_outcome
  )

  return(
    data = data
  )
}
generate_multiple_datasets_v4 <- function(
    n = 10, # 要生成的总数据集数量
    num_pleiotropic = 0,
    n_expose_heterogeneity = 0,
    N_exp = 1000, N_out = 1000,
    p_f = 0.3, p_m = 0.3,
    # --- 暴露效应 (当非零时的大小) ---
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3,
    # --- 暴露异质性 ---
    h_beta_FStoOE_exp = 0.1, h_beta_MStoOE_exp = 0.1,
    h_beta_OStoOE_exp = 0.3,
    # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
    # 子代基因型 -> 结局 (遗传叠加效应或基因多效性)
    mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0.05,
    # 父母基因型 -> 结局 (应为子代自身SNP对结局的效应，即基因多效性; 参数名可能需对应调整)
    mean_beta_OStoOE_out = 0.1, sd_beta_OStoOE_out = 0.05,
    prop_negative_pleiotropy = 0.5, # 在指定为多效性的SNP中，其效应为负的比例 (0 到 1)
    # 选型婚配
    assortative_mating_prob = 0,
    assortative_mating_strength = 0, # 选型婚配对结局的影响因子
    # 人群分层
    ## 定义人群分层的差异(次等位基因频率差异)
    crowd_stratification_differences = 0, # 用于模拟两个具有不同等位基因频率的亚群
    # --- 其他效应 ---
    beta_exp_to_out = 0, # 暴露对结局的真实因果效应
    beta_confounding_exp = 0.2, # 影响暴露的混杂因素的方差 (效应大小为1)
    beta_confounding_out = 0.2, # 影响结局的混杂因素的方差 (效应大小为1)
    correlation = 0.2, # 共享环境因素的方差
    seed = NULL) { # 随机数种子

  # --- 1. 输入验证与随机种子设置 ---
  if (!is.null(seed)) {
    set.seed(seed) # 设置随机种子以保证结果可重复性
  }

  # --- 2. 确定哪些数据集的SNP将具有潜在的多效性 ---
  n_pool <- 1:n # 创建一个从1到n的索引池
  # 从索引池中随机抽取 num_pleiotropic 个索引，这些索引对应的数据集将具有潜在的多效性SNP
  pleiotropic_indices <- sample(n_pool, num_pleiotropic)
  # 其余索引对应的数据集将作为有效的工具变量（无多效性）
  not_pleiotropic_indices <- setdiff(n_pool, pleiotropic_indices)

  expose_heterogeneity_indices <- sample(n_pool, n_expose_heterogeneity)
  not_expose_heterogeneity_indices <- setdiff(
    n_pool,
    expose_heterogeneity_indices
  )
  # --- 3. 循环生成 n 个独立的数据集 ---
  all_datasets <- vector("list", n) # 初始化一个列表以存储所有生成的数据集

  for (i in 1:n) { # 对每个数据集进行循环
    # --- 3a. 判断当前数据集 i 的特性 ---
    # 检查当前数据集的索引 i 是否在之前抽取的 "pleiotropic_indices" 中
    is_pleiotropic_candidate <- i %in%
      pleiotropic_indices # 如果为 TRUE, 则此数据集的SNP是潜在多效性的

    is_expose_heterogeneity_candidate <- i %in%
      expose_heterogeneity_indices

    # --- 3b. 为当前循环设置将要传递给底层数据生成函数的 beta 参数值 ---

    # --- 暴露效应参数 (在所有数据集中保持不变) ---
    if (is_expose_heterogeneity_candidate) {
      current_beta_FStoOE_exp <- beta_FStoOE_exp
      current_beta_MStoOE_exp <- beta_MStoOE_exp
      current_beta_OStoOE_exp <- beta_OStoOE_exp

      current_h_beta_FStoOE_exp <- h_beta_FStoOE_exp
      current_h_beta_MStoOE_exp <- h_beta_MStoOE_exp
      current_h_beta_OStoOE_exp <- h_beta_OStoOE_exp

      expose_heterogeneity_index <- "具有暴露异质性"
    } else {
      current_beta_FStoOE_exp <- beta_FStoOE_exp
      current_beta_MStoOE_exp <- beta_MStoOE_exp
      current_beta_OStoOE_exp <- beta_OStoOE_exp

      current_h_beta_FStoOE_exp <- beta_FStoOE_exp
      current_h_beta_MStoOE_exp <- beta_MStoOE_exp
      current_h_beta_OStoOE_exp <- beta_OStoOE_exp
      expose_heterogeneity_index <- "无暴露异质性"
    }



    # --- 结局效应 / 水平多效性参数 (核心修改部分) ---
    if (is_pleiotropic_candidate) {
      # 如果当前SNP是潜在多效性的:
      direction_multiplier <- if (runif(1) < prop_negative_pleiotropy) -1 else 1

      current_beta_FStoOE_out <- rnorm(1,
        mean = direction_multiplier *
          mean_beta_FStoOE_out, sd = sd_beta_FStoOE_out
      )
      current_beta_MStoOE_out <- rnorm(1,
        mean = direction_multiplier *
          mean_beta_MStoOE_out, sd = sd_beta_MStoOE_out
      )
      current_beta_OStoOE_out <- rnorm(1,
        mean = direction_multiplier *
          mean_beta_OStoOE_out, sd = sd_beta_OStoOE_out
      )

      snp_label <- "pleiotropic_variable"
    } else {
      # 如果当前SNP不是多效性的 (即，它是一个有效的工具变量):
      # 将其对结局的直接效应（多效性效应）设置为 0。
      current_beta_FStoOE_out <- 0
      current_beta_MStoOE_out <- 0
      current_beta_OStoOE_out <- 0

      snp_label <- "valid_instrument" # 标记此SNP为有效工具变量
    }

    # --- 选型婚配参数设置 ---
    # 根据 compatibility_selection_prop 概率决定当前数据集是否应用选型婚配
    is_compatibility_selection <- if (runif(1) <
      assortative_mating_prob) {
      1
    } else {
      0
    }
    if (is_compatibility_selection == 1) {
      # 如果应用选型婚配，则使用函数输入的选型婚配参数
      current_assortative_mating_strength <-
        assortative_mating_strength
    } else {
      current_assortative_mating_strength <- 1000 # 对结局的影响因子设为0
    }

    # --- 人群分层参数设置 ---
    # 为当前数据集随机选择等位基因频率，以模拟人群分层
    # current_p_f 从 (p_f, p_f + 差异) 中随机选一个
    current_p_f <- sample(c(p_f, p_f + crowd_stratification_differences), 1)
    # current_p_m 从 (p_m, p_m + 差异) 中随机选一个
    current_p_m <- sample(c(p_m, p_m + crowd_stratification_differences), 1)

    # --- 3c. 调用底层函数 (generate_mr_trio_data_2sample) 生成单个数据集 ---

    data_exp <- generate_mr_trio_data_sample(
      n = N_exp,
      p_f = current_p_f, p_m = current_p_m, # 使用当前循环的等位基因频率
      beta_FStoOE_exp = current_beta_FStoOE_exp,
      beta_MStoOE_exp = current_beta_MStoOE_exp,
      beta_OStoOE_exp = current_beta_OStoOE_exp,
      # 传递当前循环抽样生成的或设为0的结局效应值
      beta_FStoOE_out = current_beta_FStoOE_out,
      beta_MStoOE_out = current_beta_MStoOE_out,
      beta_OStoOE_out = current_beta_OStoOE_out,
      beta_exp_to_out = beta_exp_to_out,
      confounding_exp = beta_confounding_exp,
      confounding_out = beta_confounding_out,
      correlation = correlation, seed = NULL,
      assortative_mating_strength = current_assortative_mating_strength
    )


    data_out <- generate_mr_trio_data_sample(
      n = N_out,
      p_f = current_p_f, p_m = current_p_m, # 使用当前循环的等位基因频率
      beta_FStoOE_exp = current_h_beta_FStoOE_exp,
      beta_MStoOE_exp = current_h_beta_MStoOE_exp,
      beta_OStoOE_exp = current_h_beta_OStoOE_exp,
      # 传递当前循环抽样生成的或设为0的结局效应值
      beta_FStoOE_out = current_beta_FStoOE_out,
      beta_MStoOE_out = current_beta_MStoOE_out,
      beta_OStoOE_out = current_beta_OStoOE_out,
      beta_exp_to_out = beta_exp_to_out,
      confounding_exp = beta_confounding_exp,
      confounding_out = beta_confounding_out,
      correlation = correlation, seed = NULL,
      assortative_mating_strength = current_assortative_mating_strength
    )
    data_i <- list(data_exp, data_out)

    # --- 3d. 为生成的数据集添加唯一的 dataset_id 和 SNP 类型标签 ---
    # 检查底层函数返回结果的格式是否符合预期 (一个包含至少两个数据框的列表)
    if (length(data_i) >= 2 && is.data.frame(data_i[[1]]) && is.data.frame(data_i[[2]])) {
      data_i[[1]]$dataset_id <- i # 为暴露数据框添加 dataset_id
      data_i[[2]]$dataset_id <- i # 为结局数据框添加 dataset_id
      # 添加 SNP 类型标签 (例如 "pleiotropic_variable" 或 "valid_instrument")
      data_i[[1]]$snp_type <- snp_label
      data_i[[2]]$snp_type <- snp_label
      data_i[[1]]$expose_heterogeneity_index <- expose_heterogeneity_index
      data_i[[2]]$expose_heterogeneity_index <- expose_heterogeneity_index
    } else {
      warning(paste("数据集", i, "的生成结果格式不符合预期 (data_i 不是包含两个数据框的列表)，无法添加 ID 和标签。"))
    }
    all_datasets[[i]] <- data_i # 将生成的单个数据集 (data_i) 存储到总列表 all_datasets 中
  } # 结束 for 循环 (生成 n 个数据集)

  # --- 4. 合并所有独立生成的数据集 ---
  # 分别提取所有数据集中的暴露数据框和结局数据框
  all_exposure_dfs <- lapply(all_datasets, function(element) {
    if (length(element) >= 1 && is.data.frame(element[[1]])) {
      return(element[[1]]) # 返回暴露数据框
    } else {
      warning("在合并数据时，发现一个元素的结构不符合预期 (缺少暴露数据框)。")
      return(NULL) # 如果结构不符，返回NULL
    }
  })
  all_outcome_dfs <- lapply(all_datasets, function(element) {
    if (length(element) >= 2 && is.data.frame(element[[2]])) {
      return(element[[2]]) # 返回结局数据框
    } else {
      warning("在合并数据时，发现一个元素的结构不符合预期 (缺少结局数据框)。")
      return(NULL)
    }
  })

  # 移除可能的NULL元素，然后使用 rbind 合并
  all_exposure_dfs <- all_exposure_dfs[!sapply(all_exposure_dfs, is.null)]
  all_outcome_dfs <- all_outcome_dfs[!sapply(all_outcome_dfs, is.null)]

  combined_exposure_data <- if (length(all_exposure_dfs) > 0) do.call(rbind, all_exposure_dfs) else NULL
  combined_outcome_data <- if (length(all_outcome_dfs) > 0) do.call(rbind, all_outcome_dfs) else NULL



  # --- 5. 返回包含合并后的暴露和结局数据的列表 ---
  return(list(exposure_data = combined_exposure_data, outcome_data = combined_outcome_data))
}

# %% 函数测试区
if (FALSE) {
  test <- generate_multiple_datasets_v4(
    n = 10, # 要生成的总数据集数量
    num_pleiotropic = 0,
    n_expose_heterogeneity = 10,
    N_exp = 1000, N_out = 1000,
    p_f = 0.3, p_m = 0.3,
    # --- 暴露效应 (当非零时的大小) ---
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3,
    # --- 暴露异质性 ---
    h_beta_FStoOE_exp = 0, h_beta_MStoOE_exp = 0,
    h_beta_OStoOE_exp = 0.1
  )


  lm(test$outcome_data$Father_expose ~ test$outcome_data$Father_SNPs)
}

# %% 多SNPs同一数据

load_mr_trio_params <- function() {
  # 基本参数
  n <- 1000
  n_snps <<- 10
  p_f <<- 0.3
  p_m <<- 0.3

  # 暴露效应
  beta_fs_to_oe_exp <<- 0.3
  beta_ms_to_oe_exp <<- 0.3
  beta_os_to_oe_exp <<- 0.3

  # 结局效应 (直接多效性 / 遗传叠加效应)
  beta_fs_to_oe_out <<- 0
  beta_ms_to_oe_out <<- 0
  beta_os_to_oe_out <<- 0

  # 因果效应
  beta_exp_to_out <<- 0.4

  # 混杂效应
  confounding_exp <<- 0.2
  confounding_out <<- 0.2

  # 其他参数
  correlation <<- 0.2
  seed <<- NULL
  assortative_mating_strength <<- 1000


  n <<- 1000
  # 总样本量
  n_snps <<- 10
  # 要生成的总数据集数量
  n_pleiotropic <<- 0
  # 水平多效性
  n_expose_heterogeneity <<- 0
  # 暴露异质性
  p_f <<- 0.3
  p_m <<- 0.3

  # --- 暴露效应 (当非零时的大小) ---
  beta_fs_to_oe_exp <<- 0.1
  beta_ms_to_oe_exp <<- 0.1
  beta_os_to_oe_exp <<- 0.3
  # --- 暴露异质性 ---
  h_beta_fs_to_oe_exp <<- 0.1
  h_beta_ms_to_oe_exp <<- 0.1
  h_beta_os_to_oe_exp <<- 0.3
  # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
  mean_beta_fs_to_oe_out <<- 0.1
  sd_beta_fs_to_oe_out <<- 0.05
  mean_beta_ms_to_oe_out <<- 0.1
  sd_beta_ms_to_oe_out <<- 0.05
  mean_beta_os_to_oe_out <<- 0.1
  sd_beta_os_to_oe_out <<- 0.05
  p_neg_pleiotropy <<- 0.5
  # 选型婚配
  assortative_mating_strength <<- 0
  # 选型婚配对结局的影响因子
  # 人群分层
  crowd_differences <<- 0
  # 用于模拟两个具有不同等位基因频率的亚群
  # --- 其他效应 ---
  beta_exp_to_out <<- 0
  # 暴露对结局的真实因果效应
  beta_confounding_exp <<- 0.2
  # 影响暴露的混杂因素的方差 (效应大小为1)
  beta_confounding_out <<- 0.2
  # 影响结局的混杂因素的方差 (效应大小为1)
  correlation <<- 0.2
  # 共享环境因素的方差
  seed <<- NULL

  cat("所有MR Trio参数已加载到全局环境\n")
}
# load_mr_trio_params()
generate_snp_hwe_matrix <- function(n_samples, maf, n_snps = NULL) {
  # 参数验证和处理
  if (is.null(n_snps)) {
    # 如果没有指定n_snps，则根据maf向量长度确定
    if (length(maf) == 1) {
      stop("当maf为单个值时，必须指定n_snps参数")
    }
    n_snps <- length(maf)
  }

  # 如果maf是单个值，则扩展为向量
  if (length(maf) == 1) {
    maf <- rep(maf, n_snps)
  }

  # 检查maf向量长度是否与n_snps匹配
  if (length(maf) != n_snps) {
    stop("maf向量长度必须等于n_snps，或者maf为单个值")
  }

  # 检查MAF值是否在合理范围内
  if (any(maf < 0) || any(maf > 0.5)) {
    stop("MAF值必须在0到0.5之间")
  }

  # 初始化结果矩阵：行为样本，列为SNPs
  genotype_matrix <- matrix(0, nrow = n_samples, ncol = n_snps)

  # 为每个SNP生成基因型
  for (i in 1:n_snps) {
    # 计算Hardy-Weinberg平衡下的基因型频率
    # AA (0): (1-maf)^2, Aa (1): 2*maf*(1-maf), aa (2): maf^2
    genotype_freqs <- c((1 - maf[i])^2, 2 * maf[i] * (1 - maf[i]), maf[i]^2)

    # 生成基因型
    genotype_matrix[, i] <- sample(c(0, 1, 2),
      size = n_samples,
      replace = TRUE,
      prob = genotype_freqs
    )
  }

  # 添加列名和行名
  colnames(genotype_matrix) <- paste0("SNP_", 1:n_snps)
  rownames(genotype_matrix) <- paste0("Sample_", 1:n_samples)

  return(genotype_matrix)
}
generate_mr_trio_data_matrix <- function(
    n = 1000, n_snps = 10, # 每一个数据集的snps个数
    p_f = 0.3, p_m = 0.3, # p_m 当前未在SNP生成中使用
    # 暴露效应
    beta_fs_to_oe_exp = 0.3, beta_ms_to_oe_exp = 0.3,
    beta_os_to_oe_exp = 0.3,
    # 结局效应 (直接多效性 / 遗传叠加效应)
    beta_fs_to_oe_out = 0, beta_ms_to_oe_out = 0,
    beta_os_to_oe_out = 0,
    # 因果效应
    beta_exp_to_out = 0.4,
    # 混杂效应
    confounding_exp = 0.2, confounding_out = 0.2, # 混杂不用是系数
    # 其他参数
    correlation = 0.2, seed = NULL,
    # 选型婚配强度
    assortative_mating_strength = 1000) {
  set.seed(seed)


  # 输入父代的snps矩阵得到传递给子代的snps矩阵
  get_transmitted_allele <- function(parent_snps) {
    # 转换为矩阵格式（如果输入是数据框或向量）
    if (is.data.frame(parent_snps)) {
      parent_snps <- as.matrix(parent_snps)
    } else if (is.vector(parent_snps)) {
      parent_snps <- matrix(parent_snps, nrow = 1)
    }

    # 初始化结果矩阵
    transmitted_alleles <- array(NA, dim = dim(parent_snps))

    # 处理基因型为0的位点（aa -> a）
    transmitted_alleles[parent_snps == 0] <- 0

    # 处理基因型为2的位点（AA -> A）
    transmitted_alleles[parent_snps == 2] <- 1

    # 处理基因型为1的位点（Aa -> A 或 a，各50%概率）
    het_positions <- which(parent_snps == 1)
    if (length(het_positions) > 0) {
      transmitted_alleles[het_positions] <- rbinom(length(het_positions), 1, 0.5)
    }

    # 检查无效基因型
    invalid_positions <- which(!(parent_snps %in% c(0, 1, 2)) & !is.na(parent_snps))
    if (length(invalid_positions) > 0) {
      warning(paste("发现", length(invalid_positions), "个无效的基因型值"))
      transmitted_alleles[invalid_positions] <- NA
    }

    # 保持原始矩阵的行名和列名
    rownames(transmitted_alleles) <- rownames(parent_snps)
    colnames(transmitted_alleles) <- colnames(parent_snps)

    return(transmitted_alleles)
  }


  # --- 0. 生成爷爷奶奶外公外婆的snps ---

  # 父方祖父母 (爷爷、奶奶)
  grandfather_father_snps <- generate_snp_hwe_matrix(
    n_samples = n,
    maf = p_f, n_snps = n_snps
  )
  grandmother_father_snps <- generate_snp_hwe_matrix(
    n_samples = n,
    maf = p_f, n_snps = n_snps
  )
  # 母方祖父母 (外公、外婆)
  grandfather_mother_snps <- generate_snp_hwe_matrix(
    n_samples = n,
    maf = p_f, n_snps = n_snps
  )
  grandmother_mother_snps <- generate_snp_hwe_matrix(
    n_samples = n,
    maf = p_f, n_snps = n_snps
  )
  # 准备系数向量


  # --- 1. 准备系数向量，与混杂关联等变量的值 ---
  # ! 暴露组
  parameter_processing_function <- function(beta, n_snps) {
    if (is.matrix(beta)) {
      return(beta)
    } else {
      beta <- matrix(beta / n_snps, nrow = n_snps, ncol = 1)
      return(beta)
    }
  }


  beta_os_to_oe_exp_matrix <- parameter_processing_function(
    beta_os_to_oe_exp,
    n_snps
  )
  beta_fs_to_oe_exp_matrix <- parameter_processing_function(
    beta_fs_to_oe_exp,
    n_snps
  )
  beta_ms_to_oe_exp_matrix <- parameter_processing_function(
    beta_ms_to_oe_exp,
    n_snps
  )
  # ! 结局组
  beta_os_to_oe_out_matrix <- parameter_processing_function(
    beta_os_to_oe_out,
    n_snps
  )
  beta_fs_to_oe_out_matrix <- parameter_processing_function(
    beta_fs_to_oe_out,
    n_snps
  )
  beta_ms_to_oe_out_matrix <- parameter_processing_function(
    beta_ms_to_oe_out,
    n_snps
  )
  # ! 复杂与关联
  confounder_exp <- rnorm(n,
    mean = 0,
    sd = sqrt(confounding_exp)
  )
  confounder_out <- rnorm(n,
    mean = 0,
    sd = sqrt(confounding_out)
  )
  correlation_factor_exp <- rnorm(n,
    mean = 0,
    sd = sqrt(correlation)
  )
  correlation_factor_out <- rnorm(n,
    mean = 0,
    sd = sqrt(correlation)
  )


  # --- 2. 定义所有所需要的函数 ---

  generate_expose_function <- function(beta_os_to_oe_exp,
                                       beta_fs_to_oe_exp,
                                       beta_ms_to_oe_exp,
                                       snps,
                                       snps_father,
                                       snps_mother,
                                       confounder, correlation_factor) {
    # 计算暴露的确定性部分 (遗传效应 + 混杂 + 共享环境)
    results <- snps %*% beta_os_to_oe_exp +
      snps_father %*% beta_fs_to_oe_exp +
      snps_mother %*% beta_ms_to_oe_exp +
      confounder + correlation_factor
    # 计算确定性部分的方差
    var_results <- var(results, na.rm = TRUE)

    sd_err <- sqrt(max(0, 1 - var_results))

    # 添加随机误差，生成最终的暴露值
    expose <- results + rnorm(nrow(snps),
      mean = 0, sd = sd_err
    ) # 注意：rnorm的n应与snps_all长度一致
    return(expose)
  }

  generate_outcome_function <- function(beta_exp_to_out,
                                        beta_os_to_oe_out,
                                        beta_fs_to_oe_out,
                                        beta_ms_to_oe_out,
                                        expose, snps,
                                        snps_father, snps_mother,
                                        confounder, correlation_factor) {
    # 计算结局的确定性部分
    outcome_deterministic <- beta_exp_to_out * expose +
      snps %*% beta_os_to_oe_out +
      snps_father %*% beta_fs_to_oe_out +
      snps_mother %*% beta_ms_to_oe_out +
      confounder + correlation_factor

    # 计算确定性部分的方差
    var_results <- var(outcome_deterministic, na.rm = TRUE)
    # 计算随机误差的标准差，以使最终结局的方差为1
    sd_err <- sqrt(max(0, 1 - var_results))
    # 添加随机误差，生成最终的结局值
    outcome_all <- outcome_deterministic + rnorm(nrow(snps),
      mean = 0, sd = sd_err
    ) # 注意：rnorm的N应与SNPs_self长度一致
    return(outcome_all)
  }

  assortative_mating_function <- function(snps_1, expose_1, outcome_1, snps_2,
                                          expose_2, outcome_2, am_strength) {
    # 检查输入维度一致性
    if (is.matrix(snps_1) || is.data.frame(snps_1)) {
      n_individuals_1 <- nrow(snps_1)
    } else {
      n_individuals_1 <- length(snps_1)
    }

    if (is.matrix(snps_2) || is.data.frame(snps_2)) {
      n_individuals_2 <- nrow(snps_2)
    } else {
      n_individuals_2 <- length(snps_2)
    }

    # 检查样本数是否一致
    if (n_individuals_1 != length(expose_1) ||
      n_individuals_2 != length(expose_2) ||
      n_individuals_1 != n_individuals_2) {
      stop("SNPs矩阵的行数必须与暴露变量的长度一致，且两组样本数必须相等")
    }

    n <- n_individuals_1

    # 创建数据对象
    # 对于SNPs矩阵，我们需要保持矩阵结构
    if (is.matrix(snps_1) || is.data.frame(snps_1)) {
      # SNPs是矩阵的情况
      object_1 <- list(
        snps = snps_1,
        expose = expose_1,
        outcome = outcome_1,
        expose_star = expose_1 + rnorm(n, mean = 0, sd = sqrt(am_strength)),
        original_index = 1:n
      )

      object_2 <- list(
        snps = snps_2,
        expose = expose_2,
        outcome = outcome_2,
        outcome_star = outcome_2 + rnorm(n, mean = 0, sd = sqrt(am_strength)),
        original_index = 1:n
      )
    }
    # 根据expose_star排序
    order_1 <- order(object_1$expose_star)
    order_2 <- order(object_2$outcome_star)

    # 重新排列数据
    if (is.matrix(snps_1) || is.data.frame(snps_1)) {
      # 矩阵情况：需要重新排列每个组件
      object_1_sorted <- list(
        snps = object_1$snps[order_1, , drop = FALSE],
        expose = object_1$expose[order_1],
        outcome = object_1$outcome[order_1],
        expose_star = object_1$expose_star[order_1],
        original_index = object_1$original_index[order_1]
      )

      object_2_sorted <- list(
        snps = object_2$snps[order_2, , drop = FALSE],
        expose = object_2$expose[order_2],
        outcome = object_2$outcome[order_2],
        outcome_star = object_2$outcome_star[order_2],
        original_index = object_2$original_index[order_2]
      )

      # 保持原始的行名和列名
      if (!is.null(rownames(snps_1))) {
        rownames(object_1_sorted$snps) <- rownames(snps_1)[order_1]
        rownames(object_2_sorted$snps) <- rownames(snps_2)[order_2]
      }
    }

    return(list(object_1_sorted, object_2_sorted))
  }

  # 辅助函数：提取排序后的配对信息
  get_mating_pairs <- function(mating_result) {
    object_1 <- mating_result[[1]]
    object_2 <- mating_result[[2]]

    if (is.list(object_1) && !is.data.frame(object_1)) {
      # 矩阵情况
      pairs <- data.frame(
        pair_id = 1:length(object_1$expose),
        parent1_original_index = object_1$original_index,
        parent2_original_index = object_2$original_index,
        parent1_expose = object_1$expose,
        parent2_expose = object_2$expose,
        parent1_expose_star = object_1$expose_star,
        parent2_expose_star = object_2$expose_star
      )
    } else {
      # 向量情况
      pairs <- data.frame(
        pair_id = 1:nrow(object_1),
        parent1_original_index = object_1$original_index,
        parent2_original_index = object_2$original_index,
        parent1_expose = object_1$expose,
        parent2_expose = object_2$expose,
        parent1_expose_star = object_1$expose_star,
        parent2_expose_star = object_2$expose_star
      )
    }

    return(pairs)
  }

  # --- 3. 生成第1代的暴露和结局 ---
  beta_zero_matrix <- rep(0, n_snps)
  # 生成暴露
  grandfather_father_expose <- generate_expose_function(
    beta_os_to_oe_exp_matrix,
    beta_zero_matrix,
    beta_zero_matrix,
    grandfather_father_snps,
    grandfather_father_snps,
    grandfather_father_snps,
    confounder_exp,
    correlation_factor_exp
  )
  grandmother_father_expose <- generate_expose_function(
    beta_os_to_oe_exp_matrix,
    beta_zero_matrix,
    beta_zero_matrix,
    grandmother_father_snps,
    grandmother_father_snps,
    grandmother_father_snps,
    confounder_exp, correlation_factor_exp
  )
  grandfather_mother_expose <- generate_expose_function(
    beta_os_to_oe_exp_matrix, beta_zero_matrix, beta_zero_matrix,
    grandfather_mother_snps, grandfather_mother_snps, grandfather_mother_snps,
    confounder_exp, correlation_factor_exp
  )
  grandmother_mother_expose <- generate_expose_function(
    beta_os_to_oe_exp_matrix, beta_zero_matrix, beta_zero_matrix,
    grandmother_mother_snps, grandmother_mother_snps, grandmother_mother_snps,
    confounder_exp, correlation_factor_exp
  )

  # 生成结局
  grandfather_father_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_os_to_oe_out_matrix, beta_zero_matrix,
    beta_zero_matrix,
    grandfather_father_expose,
    grandfather_father_snps, grandfather_father_snps, grandfather_father_snps,
    confounder_out, correlation_factor_out
  )
  grandmother_father_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_os_to_oe_out_matrix, beta_zero_matrix,
    beta_zero_matrix,
    grandmother_father_expose,
    grandmother_father_snps, grandmother_father_snps, grandmother_father_snps,
    confounder_out, correlation_factor_out
  )
  grandfather_mother_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_os_to_oe_out_matrix,
    beta_zero_matrix,
    beta_zero_matrix,
    grandfather_mother_expose,
    grandfather_mother_snps, grandfather_mother_snps, grandfather_mother_snps,
    confounder_out, correlation_factor_out
  )
  grandmother_mother_outcome <- generate_outcome_function(
    beta_exp_to_out,
    beta_os_to_oe_out_matrix,
    beta_zero_matrix,
    beta_zero_matrix,
    grandmother_mother_expose,
    grandmother_mother_snps,
    grandmother_mother_snps,
    grandmother_mother_snps,
    confounder_out,
    correlation_factor_out
  )


  # --- 4. 对祖父母辈进行选型婚配 ---
  grand_father_am_results <- assortative_mating_function(
    grandfather_father_snps, grandfather_father_expose,
    grandfather_father_outcome,
    grandmother_father_snps, grandmother_father_expose,
    grandmother_father_outcome,
    assortative_mating_strength
  )
  grand_mother_am_results <- assortative_mating_function(
    grandfather_mother_snps, grandfather_mother_expose,
    grandfather_mother_outcome,
    grandmother_mother_snps, grandmother_mother_expose,
    grandmother_mother_outcome,
    assortative_mating_strength
  )
  grandfather_father_snps <- grand_father_am_results[[1]]$snps
  grandmother_father_snps <- grand_father_am_results[[2]]$snps
  grandfather_mother_snps <- grand_mother_am_results[[1]]$snps
  grandmother_mother_snps <- grand_mother_am_results[[2]]$snps


  # --- 7. 生成父母代基因型 (通过模拟祖父母向父母传递等位基因) ---

  father_snps <- get_transmitted_allele(grandfather_father_snps) +
    get_transmitted_allele(grandmother_father_snps)

  mother_snps <- get_transmitted_allele(grandfather_mother_snps) +
    get_transmitted_allele(grandmother_mother_snps)

  # --- 8. 为父母代生成暴露数据 (用于他们之间的选型婚配) ---
  # 父亲的暴露，受其自身SNP以及其父母(即子代的祖父母)SNP的遗传叠加效应影响
  father_expose <- generate_expose_function(
    beta_os_to_oe_exp_matrix, beta_fs_to_oe_exp_matrix,
    beta_ms_to_oe_exp_matrix, # 自身、其父(爷爷)、其母(奶奶)的snp效应
    father_snps, grandfather_father_snps, grandmother_father_snps,
    confounder_exp, correlation_factor_exp
  )
  # 母亲的暴露，受其自身snp以及其父母(即子代的外祖父母)snp的遗传叠加效应影响
  mother_expose <- generate_expose_function(
    beta_os_to_oe_exp_matrix, beta_fs_to_oe_exp_matrix,
    beta_ms_to_oe_exp_matrix, # 自身、其父(外公)、其母(外婆)的snp效应
    mother_snps, grandfather_mother_snps,
    grandmother_mother_snps, # 注意这里应该是祖父母的snps
    confounder_exp, correlation_factor_exp
  )

  # --- 12. 为父亲、母亲和子代生成结局数据 ---
  # 父亲的结局
  father_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_os_to_oe_out_matrix,
    beta_fs_to_oe_out_matrix, beta_ms_to_oe_out_matrix,
    father_expose, father_snps,
    grandfather_father_snps, grandmother_father_snps, # 父亲的"父母"是爷爷奶奶
    confounder_out, correlation_factor_out
  )
  # 母亲的结局
  mother_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_os_to_oe_out_matrix,
    beta_fs_to_oe_out_matrix, beta_ms_to_oe_out_matrix,
    mother_expose, mother_snps,
    grandfather_mother_snps, grandmother_mother_snps, # 母亲的"父母"是外公外婆
    confounder_out, correlation_factor_out
  )
  # --- 9. 对父母代进行选型婚配 ---
  parent_am_results <- assortative_mating_function(
    father_snps, father_expose, father_outcome, # 父亲的信息
    mother_snps, mother_expose, mother_outcome, # 母亲的信息
    assortative_mating_strength
  )
  # 更新父母的SNP和暴露向量，使其根据选型婚配的结果重新排序
  father_snps <- parent_am_results[[1]]$snps # 排序后的父亲snp
  father_expose <- parent_am_results[[1]]$expose
  father_outcome <- parent_am_results[[1]]$outcome
  mother_snps <- parent_am_results[[2]]$snps # 排序后并与父亲配对的母亲snp
  mother_expose <- parent_am_results[[2]]$expose
  mother_outcome <- parent_am_results[[2]]$outcome

  # --- 10. 生成子代基因型和暴露数据 ---
  # 子代的基因型 = 父亲传递的等位基因 + 母亲传递的等位基因
  offspring_snps <- get_transmitted_allele(father_snps) +
    get_transmitted_allele(mother_snps)

  # 子代的暴露，受其自身SNP以及其父母SNP的遗传叠加效应影响
  offspring_expose <- generate_expose_function(
    beta_os_to_oe_exp_matrix,
    beta_fs_to_oe_exp_matrix,
    beta_ms_to_oe_exp_matrix, # 自身、父亲、母亲的snp效应
    offspring_snps, father_snps, mother_snps,
    confounder_exp, correlation_factor_exp
  )
  offspring_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_os_to_oe_out_matrix,
    beta_fs_to_oe_out_matrix, beta_ms_to_oe_out_matrix,
    offspring_expose, offspring_snps,
    father_snps, mother_snps, # 子代的父母是father和mother
    confounder_out, correlation_factor_out
  )
  data <- data.frame(
    father_snps = father_snps,
    mother_snps = mother_snps,
    offspring_snps = offspring_snps,
    father_expose = father_expose,
    mother_expose = mother_expose,
    offspring_expose = offspring_expose,
    father_outcome = father_outcome,
    mother_outcome = mother_outcome,
    offspring_outcome = offspring_outcome
  )

  return(data)
}
if (FALSE) {
  beta_ms_to_oe_exp <- matrix(0, nrow = 10, ncol = 1)
  test1 <- generate_mr_trio_data_matrix( beta_fs_to_oe_exp = rep(0.5), beta_ms_to_oe_exp = 0.5, beta_os_to_oe_exp = 0.5)
  str(test1)
  lm(test1$father_expose ~ test1$father_snps.SNP_2)
  test_1 <- generate_mr_trio_data_matrix_ultra()

  fgwas_for_data_matrix(test_1[[1]])
}
generate_mr_trio_data_matrix_ultra <- function(
    n = 1000, # 总样本量
    n_snps = 10, # 要生成的总数据集数量
    n_pleiotropic = 0, # 水平多效性
    n_expose_heterogeneity = 0, # 暴露异质性
    p_f = 0.3, p_m = 0.3,
    # --- 暴露效应 (当非零时的大小) ---
    beta_fs_to_oe_exp = 0.1, beta_ms_to_oe_exp = 0.1,
    beta_os_to_oe_exp = 0.3,
    # --- 暴露异质性 ---
    h_beta_fs_to_oe_exp = 0.2, h_beta_ms_to_oe_exp = 0.1,
    h_beta_os_to_oe_exp = 0.3,
    # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
    mean_beta_fs_to_oe_out = 0.1, sd_beta_fs_to_oe_out = 0.05,
    mean_beta_ms_to_oe_out = 0.1, sd_beta_ms_to_oe_out = 0.05,
    mean_beta_os_to_oe_out = 0.1, sd_beta_os_to_oe_out = 0.05,
    p_neg_pleiotropy = 0.5,
    # 选型婚配
    assortative_mating_strength = 0, # 选型婚配对结局的影响因子
    # 人群分层
    crowd_differences = 0, # 用于模拟两个具有不同等位基因频率的亚群
    # --- 其他效应 ---
    beta_exp_to_out = 0, # 暴露对结局的真实因果效应
    confounding_exp = 0.2, # 影响暴露的混杂因素的方差 (效应大小为1)
    confounding_out = 0.2, # 影响结局的混杂因素的方差 (效应大小为1)
    correlation = 0.2, # 共享环境因素的方差
    seed = NULL) {
  # --- 1. 生成水平多效性的系数 ---
  beta_out_processing <- function(mean_beta_out, sd_beta_out, ...) {
    beta_out_matrix <- rnorm(n_snps,
      mean = mean_beta_out, sd = sd_beta_out
    )
    beta_out_nozero <- sample(1:n_snps, n_pleiotropic)
    if (length(beta_out_nozero) == 0) {
      beta_out_matrix[] <- 0
    } else {
      beta_out_matrix[-beta_out_nozero] <- 0
    }

    n_neg <- floor(n_pleiotropic * p_neg_pleiotropy)
    beta_out_neg <- sample(beta_out_nozero, n_neg)
    beta_out_matrix[beta_out_neg] <- -beta_out_matrix[beta_out_neg]
    beta_out_matrix <- matrix(beta_out_matrix, nrow = n_snps, ncol = 1)
    return(beta_out_matrix)
  }
  beta_fs_to_oe_out_matrix <- beta_out_processing(
    mean_beta_fs_to_oe_out, sd_beta_fs_to_oe_out,
    n_snps, n_pleiotropic, p_neg_pleiotropy
  )
  beta_ms_to_oe_out_matrix <- beta_out_processing(
    mean_beta_ms_to_oe_out, sd_beta_ms_to_oe_out,
    n_snps, n_pleiotropic, p_neg_pleiotropy
  )
  beta_os_to_oe_out_matrix <- beta_out_processing(
    mean_beta_os_to_oe_out, sd_beta_os_to_oe_out,
    n_snps, n_pleiotropic, p_neg_pleiotropy
  )
  # --- 2. 生成间接遗传异质性的系数 ---
  beta_exp_processing <- function(beta_exp,
                                  h_beta_exp, n_snps, n_pleiotropic) {
    beta_exp_matrix <- rep(beta_exp, n_snps)
    h_beta_exp_index <- sample(1:n_snps, n_pleiotropic)
    beta_exp_matrix[h_beta_exp_index] <- h_beta_exp
    beta_exp_matrix <- as.matrix(beta_exp_matrix)
    return(beta_exp_matrix)
  }
  beta_fs_to_oe_exp_matrix <- beta_exp_processing(
    beta_fs_to_oe_exp,
    h_beta_fs_to_oe_exp, n_snps, 0
  )
  beta_ms_to_oe_exp_matrix <- beta_exp_processing(
    beta_ms_to_oe_exp,
    h_beta_fs_to_oe_exp, n_snps, 0
  )
  beta_os_to_oe_exp_matrix <- beta_exp_processing(
    beta_os_to_oe_exp,
    h_beta_fs_to_oe_exp, n_snps, 0
  )
  h_beta_fs_to_oe_exp_matrix <- beta_exp_processing(
    beta_fs_to_oe_exp,
    h_beta_fs_to_oe_exp, n_snps, n_pleiotropic
  )
  h_beta_ms_to_oe_exp_matrix <- beta_exp_processing(
    beta_ms_to_oe_exp,
    h_beta_ms_to_oe_exp, n_snps, n_pleiotropic
  )
  h_beta_os_to_oe_exp_matrix <- beta_exp_processing(
    beta_os_to_oe_exp,
    h_beta_os_to_oe_exp, n_snps, n_pleiotropic
  )
  # --- 3. 生成两个数据 ---
  data_exp <- generate_mr_trio_data_matrix(
    n = n, n_snps = n_snps, # 每一个数据集的snps个数
    p_f = p_f, p_m = p_m, # p_m 当前未在SNP生成中使用
    # 暴露效应
    beta_fs_to_oe_exp = beta_fs_to_oe_exp_matrix,
    beta_ms_to_oe_exp = beta_ms_to_oe_exp_matrix,
    beta_os_to_oe_exp = beta_os_to_oe_exp_matrix,
    # 结局效应 (直接多效性 / 遗传叠加效应)
    beta_fs_to_oe_out = beta_fs_to_oe_out_matrix,
    beta_ms_to_oe_out = beta_ms_to_oe_out_matrix,
    beta_os_to_oe_out = beta_os_to_oe_out_matrix,
    # 因果效应
    beta_exp_to_out = beta_exp_to_out,
    # 混杂效应
    confounding_exp = confounding_exp,
    confounding_out = confounding_out, # 混杂不用是系数
    # 其他参数
    correlation = correlation, seed = seed,
    # 选型婚配强度
    assortative_mating_strength = assortative_mating_strength
  )

  data_out <- generate_mr_trio_data_matrix(
    n = n, n_snps = n_snps, # 每一个数据集的snps个数
    p_f = p_f, p_m = p_m, # p_m 当前未在SNP生成中使用
    # 暴露效应
    beta_fs_to_oe_exp = h_beta_fs_to_oe_exp_matrix,
    beta_ms_to_oe_exp = h_beta_ms_to_oe_exp_matrix,
    beta_os_to_oe_exp = h_beta_os_to_oe_exp_matrix,
    # 结局效应 (直接多效性 / 遗传叠加效应)
    beta_fs_to_oe_out = beta_fs_to_oe_out_matrix,
    beta_ms_to_oe_out = beta_ms_to_oe_out_matrix,
    beta_os_to_oe_out = beta_os_to_oe_out_matrix,
    # 因果效应
    beta_exp_to_out = beta_exp_to_out,
    # 混杂效应
    confounding_exp = confounding_exp,
    confounding_out = confounding_out, # 混杂不用是系数
    # 其他参数
    correlation = correlation, seed = seed,
    # 选型婚配强度
    assortative_mating_strength = assortative_mating_strength
  )
  results <- list(data_exp = data_exp, data_out = data_out)
  return(results)
}
