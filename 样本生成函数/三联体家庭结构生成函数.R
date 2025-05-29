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
generate_snp_hwe <- function(n_samples, maf) {
  q <- maf
  p <- 1 - q
  genotype_probs <- c(p^2, 2 * p * q, q^2)
  simulated_genotypes <- as.matrix(sample(c(0, 1, 2),
    size = n_samples,
    replace = TRUE,
    prob = genotype_probs
  ), ncol = 1)

  return(simulated_genotypes)
}

#' @title 模拟包含选型婚配的三代家系孟德尔随机化双样本数据
#' @description 为双样本孟德尔随机化(MR)研究生成模拟数据。
#'              该函数模拟了从祖父母到父母再到子代三代的遗传效应，
#'              并考虑了选型婚配、遗传叠加效应（dynastic effects）和混杂因素。
#'
#' @param N_exp 整数型。暴露GWAS的样本量。
#' @param N_out 整数型。结局GWAS的样本量。
#' @param overlap_prop 数值型。暴露和结局样本之间的重叠比例 (0 到 1)。
#' @param p_f 数值型。用于生成SNP的等位基因频率 (例如，次要等位基因频率MAF)。
#'            当前用于所有祖父母辈的SNP生成。
#' @param p_m 数值型。等位基因频率 (注意：在当前版本的SNP生成逻辑中，此参数 p_m 未被使用，
#'            所有祖父母的SNP均使用 p_f 生成)。
#'
#' @section 暴露的遗传效应参数:
#' @param beta_FStoOE_exp 数值型。父亲SNP对子代暴露的效应 (暴露的遗传叠加效应)。
#' @param beta_MStoOE_exp 数值型。母亲SNP对子代暴露的效应 (暴露的遗传叠加效应)。
#' @param beta_OStoOE_exp 数值型。子代自身SNP对子代暴露的效应。
#'
#' @section 结局的遗传效应参数 (直接多效性 / 结局的遗传叠加效应):
#' @param beta_FStoOE_out 数值型。父亲SNP对子代结局的效应。
#' @param beta_MStoOE_out 数值型。母亲SNP对子代结局的效应。
#' @param beta_OStoOE_out 数值型。子代自身SNP对子代结局的效应 (基因多效性)。
#'
#' @section 因果效应参数:
#' @param beta_exp_to_out 数值型。暴露对结局的真实因果效应。
#'
#' @section 混杂效应参数:
#' @param confounding_exp 数值型。影响暴露的混杂因素的方差。
#'                        (混杂因素的效应大小隐式设为1)。
#' @param confounding_out 数值型。影响结局的混杂因素的方差。
#'                        (混杂因素的效应大小隐式设为1)。
#'
#' @section 其他参数:
#' @param correlation 数值型。用于模拟共享环境 (例如家庭环境) 的方差组成部分，
#'                    该环境相似地影响暴露和结局。
#' @param seed 整数型。用于随机数生成的种子，以确保结果可重复。
#' @param assortative_mating_strength 数值型。选型婚配强度参数。具体而言，这是添加到真实暴露值上
#'                                    以创建用于选型婚配的代理表型的“噪音”的方差。
#'                                    如果此值较大，代理表型噪音的標準差 (`sqrt(assortative_mating_strength)`)
#'                                    相对于暴露本身的变异较大时，基于真实暴露的选型会较弱。
#' @param sample.outcome.betas.from.range 逻辑型。一个占位符或未来可能使用的功能，当前函数体中未使用。
#'
#' @return 一个包含两个数据框的列表:
#'         \describe{
#'           \item{exposure_data}{暴露GWAS样本的数据框，包含ID、父母代SNP、子代SNP，
#'                                以及这三代个体的暴露值。}
#'           \item{outcome_data}{结局GWAS样本的数据框，包含ID、父母代SNP、子代SNP，
#'                               以及这三代个体的结局值。}
#'         }
#' @md
generate_mr_trio_data_2sample <- function(
    N_exp = 1000, N_out = 1000, overlap_prop = 0,
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

  # 计算模拟所需的个体总数 (N_total 用于生成一个完整的人群库，后续从中抽样)
  N_total <- N_exp + N_out # 注意：这里的N_total是暴露和结局样本量之和，用于后续抽样。
  # 如果有重叠，实际独立个体数会少于N_exp + N_out。
  # 但为了简化选型婚配和家系结构模拟，我们先生成一个较大的人群。

  # 为所有 N_total 个体生成唯一ID或索引
  all_indices <- 1:N_total


  #' @title (内部辅助函数) 生成符合哈迪-温伯格平衡的SNP基因型
  #' @description 根据指定的等位基因频率 p 生成 N 个个体的SNP基因型。
  #'              基因型编码: 2 代表 AA (频率 p^2), 1 代表 Aa (频率 2pq), 0 代表 aa (频率 q^2),
  #'              其中 q = 1-p。参数 'p' 在这里指代的是编码为'2'的纯合子所对应的等位基因频率。
  #' @param N 整数型。需要生成的样本量。
  #' @param p 数值型。等位基因'A'的频率 (0到1之间)。
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

  #' @title (内部辅助函数) 模拟从单个亲本传递给子代的等位基因
  #' @description 根据亲本的SNP基因型，模拟孟德尔遗传中传递给子代的等位基因。
  #' @param parent_snps 数值向量。亲本的基因型 (0, 1, 或 2)。
  #'                    编码: 0=aa, 1=Aa, 2=AA (a是等位基因0, A是等位基因1)。
  #' @return 一个与 parent_snps 等长的数值向量，包含传递的等位基因 (0 或 1)。
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

  # --- 1. (外部SNP生成函数调用) 为祖父母辈生成SNP数据 ---
  # 假设 generate_snp_hwe 是一个外部定义的函数，它接受样本量和MAF作为参数
  # 并且其基因型编码与 internal_get_transmitted_allele 所期望的一致 (例如，0=纯合次要，1=杂合，2=纯合主要)
  # 或者 internal_get_transmitted_allele 的注释需要根据 generate_snp_hwe 的输出来调整。
  # 这里我们假设 generate_snp_hwe 返回的基因型编码是：0=aa, 1=Aa, 2=AA, maf是a的频率
  # 因此，p_f 代表的是等位基因'a'的频率(次要等位基因频率)

  # 父方祖父母 (爷爷、奶奶)
  Grandfather_Father_SNPs <- generate_snp_hwe(n_samples = N_total, maf = p_f)
  Grandmother_Father_SNPs <- generate_snp_hwe(n_samples = N_total, maf = p_f)

  # 母方祖父母 (外公、外婆)
  Grandfather_Mother_SNPs <- generate_snp_hwe(n_samples = N_total, maf = p_f)
  Grandmother_Mother_SNPs <- generate_snp_hwe(n_samples = N_total, maf = p_f)


  # --- 2. 定义暴露和结局的生成函数 ---

  #' @title (内部辅助函数) 生成暴露表型数据
  #' @description 根据个体自身SNP、父母SNP(遗传叠加效应)、混杂因素和共享环境因素生成暴露表型。
  #'              生成的暴露值会进行标准化，使得其总方差近似为1。
  #' @param beta_StoE_exp 数值型。个体自身SNP对暴露的效应系数。
  #' @param beta_FStoOE_exp 数值型。父亲SNP对个体(子代)暴露的效应系数。
  #' @param beta_MStoOE_exp 数值型。母亲SNP对个体(子代)暴露的效应系数。
  #' @param SNPs_all 数值向量。个体自身的SNP基因型。
  #' @param SNPs_Father 数值向量。个体父亲的SNP基因型。
  #' @param SNPs_Mother 数值向量。个体母亲的SNP基因型。
  #' @param confounder 数值向量。影响暴露的混杂因素。
  #' @param correlation_factor 数值向量。影响暴露的共享环境因素。
  #' @return 数值向量。生成的暴露表型值。
  generate_expose_function <- function(beta_StoE_exp, beta_FStoOE_exp, beta_MStoOE_exp,
                                       SNPs_all, SNPs_Father, SNPs_Mother,
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

  # --- 3. 为所有 N_total 个体生成混杂因素和共享环境因素 ---
  # 这些因素将在后续生成暴露和结局时使用
  # 影响暴露的混杂因素 (均值为0, 方差由 confounding_exp 定义)
  confounder_exp_values <- rnorm(N_total, mean = 0, sd = sqrt(confounding_exp))
  # 影响结局的混杂因素 (均值为0, 方差由 confounding_out 定义)
  confounder_out_values <- rnorm(N_total, mean = 0, sd = sqrt(confounding_out))

  # 共享环境因素 (例如，家庭环境，方差由 correlation 定义)
  # 假设暴露和结局受不同(但同分布的)共享环境因素影响，或这是同一个因素的不同实现
  correlation_factor_exp_values <- rnorm(N_total, mean = 0, sd = sqrt(correlation)) # 用于暴露
  correlation_factor_out_values <- rnorm(N_total, mean = 0, sd = sqrt(correlation)) # 用于结局

  # --- 4. 为祖父母辈生成暴露数据 (用于后续的选型婚配) ---
  # 注意: 对祖父母而言，他们的“父母”SNP效应参数(beta_FStoOE_exp, beta_MStoOE_exp)设为0，
  #       并且他们的“父亲”和“母亲”SNP都用自身SNP代替，因为我们不模拟曾祖父母。
  # 父方爷爷的暴露
  Grandfather_Father_expose <- generate_expose_function(
    beta_OStoOE_exp, 0, 0, # 自身SNP效应, 无父母SNP效应
    Grandfather_Father_SNPs, Grandfather_Father_SNPs, Grandfather_Father_SNPs, # 自身SNP作为父母SNP传入
    confounder_exp_values, correlation_factor_exp_values
  )
  # 父方奶奶的暴露
  Grandmother_Father_expose <- generate_expose_function(
    beta_OStoOE_exp, 0, 0,
    Grandmother_Father_SNPs, Grandmother_Father_SNPs, Grandmother_Father_SNPs,
    confounder_exp_values, correlation_factor_exp_values
  )
  # 母方外公的暴露
  Grandfather_Mother_expose <- generate_expose_function(
    beta_OStoOE_exp, 0, 0,
    Grandfather_Mother_SNPs, Grandfather_Mother_SNPs, Grandfather_Mother_SNPs,
    confounder_exp_values, correlation_factor_exp_values
  )
  # 母方外婆的暴露
  Grandmother_Mother_expose <- generate_expose_function(
    beta_OStoOE_exp, 0, 0,
    Grandmother_Mother_SNPs, Grandmother_Mother_SNPs, Grandmother_Mother_SNPs,
    confounder_exp_values, correlation_factor_exp_values
  )

  # --- 5. 定义选型婚配函数 ---
  #' @title (内部辅助函数) 执行选型婚配
  #' @description 根据个体的暴露代理表型对两组个体进行排序，以模拟选型婚配。
  #'              选型婚配强度通过添加到真实暴露值上的噪音方差来控制。
  #' @param snps_1 数值向量。第一组个体的SNP基因型。
  #' @param expose_1 数值向量。第一组个体的真实暴露值。
  #' @param snps_2 数值向量。第二组个体的SNP基因型。
  #' @param expose_2 数值向量。第二组个体的真实暴露值。
  #' @param am_strength 数值型。选型婚配强度参数 (噪音的方差)。
  #' @return 一个列表，包含两个排序后的数据框 (object_1, object_2)。
  #'         每个数据框包含排序后的 snps 和 expose，以及用于排序的 expose_star。
  assortative_mating_function <- function(snps_1, expose_1, snps_2, expose_2,
                                          am_strength) {
    # 创建包含SNP和暴露的数据框
    object_1 <- data.frame(snps = snps_1, expose = expose_1)
    object_2 <- data.frame(snps = snps_2, expose = expose_2)

    # 生成用于选型婚配的暴露代理表型 (真实暴露 + 随机噪音)
    # 噪音的标准差是 sqrt(am_strength)
    object_1 <- object_1 %>% mutate(expose_star = expose +
      rnorm(n(), mean = 0, sd = sqrt(am_strength)))
    object_2 <- object_2 %>% mutate(expose_star = expose +
      rnorm(n(), mean = 0, sd = sqrt(am_strength)))

    # 根据暴露代理表型对两组个体分别进行升序排序
    object_1 <- arrange(object_1, expose_star)
    object_2 <- arrange(object_2, expose_star)
    # 排序后，object_1的第一行将与object_2的第一行配对，以此类推。
    return(list(object_1, object_2))
  }

  # --- 6. 对祖父母辈进行选型婚配 ---
  # 父系祖父母 (爷爷奶奶) 进行选型婚配
  grand_father_am_results <- assortative_mating_function(
    Grandfather_Father_SNPs, Grandfather_Father_expose, # 爷爷的信息
    Grandmother_Father_SNPs, Grandmother_Father_expose, # 奶奶的信息
    assortative_mating_strength
  )
  # 母系祖父母 (外公外婆) 进行选型婚配
  grand_mother_am_results <- assortative_mating_function(
    Grandfather_Mother_SNPs, Grandfather_Mother_expose, # 外公的信息
    Grandmother_Mother_SNPs, Grandmother_Mother_expose, # 外婆的信息
    assortative_mating_strength
  )

  # 更新祖父母的SNP向量，使其根据选型婚配的结果重新排序
  # 这样，后续从这些排序后的祖父母传递等位基因时，就体现了选型婚配的效果
  Grandfather_Father_SNPs <- grand_father_am_results[[1]]$snps # 排序后的爷爷SNP
  Grandmother_Father_SNPs <- grand_father_am_results[[2]]$snps # 排序后并与爷爷配对的奶奶SNP

  Grandfather_Mother_SNPs <- grand_mother_am_results[[1]]$snps # 排序后的外公SNP
  Grandmother_Mother_SNPs <- grand_mother_am_results[[2]]$snps # 排序后并与外公配对的外婆SNP

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

  # --- 9. 对父母代进行选型婚配 ---
  parent_am_results <- assortative_mating_function(
    Father_SNPs, Father_expose, # 父亲的信息
    Mother_SNPs, Mother_expose, # 母亲的信息
    assortative_mating_strength
  )
  # 更新父母的SNP和暴露向量，使其根据选型婚配的结果重新排序
  Father_SNPs <- parent_am_results[[1]]$snps # 排序后的父亲SNP
  Mother_SNPs <- parent_am_results[[2]]$snps # 排序后并与父亲配对的母亲SNP
  Father_expose <- parent_am_results[[1]]$expose # 排序后的父亲暴露
  Mother_expose <- parent_am_results[[2]]$expose # 排序后并与父亲配对的母亲暴露

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

  # --- 11. 定义结局生成函数 ---
  #' @title (内部辅助函数) 生成结局表型数据
  #' @description 根据暴露值、个体自身SNP、父母SNP(遗传叠加效应/多效性)、混杂因素和共享环境因素生成结局表型。
  #'              生成的结局值会进行标准化，使得其总方差近似为1。
  #' @param beta_exp_to_out_param 数值型。暴露对结局的因果效应系数。 (参数名修改以避免与外部变量冲突)
  #' @param beta_OStoOE_out_param 数值型。个体自身SNP对结局的效应系数 (基因多效性)。
  #' @param beta_FStoOE_out_param 数值型。父亲SNP对个体(子代)结局的效应系数。
  #' @param beta_MStoOE_out_param 数值型。母亲SNP对个体(子代)结局的效应系数。
  #' @param expose_values 数值向量。个体的暴露值。 (参数名修改)
  #' @param SNPs_self 数值向量。个体自身的SNP基因型。 (参数名修改)
  #' @param SNPs_father_param 数值向量。个体父亲的SNP基因型。(参数名修改)
  #' @param SNPs_mother_param 数值向量。个体母亲的SNP基因型。(参数名修改)
  #' @param confounder_values 数值向量。影响结局的混杂因素。 (参数名修改)
  #' @param correlation_factor_values 数值向量。影响结局的共享环境因素。 (参数名修改)
  #' @return 数值向量。生成的结局表型值。
  generate_outcome_function <- function(beta_exp_to_out_param, beta_OStoOE_out_param,
                                        beta_FStoOE_out_param, beta_MStoOE_out_param,
                                        expose_values, SNPs_self, SNPs_father_param, SNPs_mother_param,
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
    outcome_all <- outcome_deterministic + rnorm(length(SNPs_self), mean = 0, sd = sd_err) # 注意：rnorm的N应与SNPs_self长度一致
    return(outcome_all)
  }

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
  # 子代的结局
  Offspring_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
    Offspring_expose, Offspring_SNPs, Father_SNPs, Mother_SNPs, # 子代的父母是Father和Mother
    confounder_out_values, correlation_factor_out_values
  )

  # --- 13. 创建包含所有模拟信息的完整数据集 ---
  # 这个数据集包含了 N_total 个家庭（父亲-母亲-子代）的完整信息
  data_all <- data.frame(
    ID = all_indices, # 个体/家庭的唯一ID，方便后续追踪和抽样
    Father_SNPs = Father_SNPs,
    Mother_SNPs = Mother_SNPs,
    Offspring_SNPs = Offspring_SNPs,
    Father_expose = Father_expose,
    Mother_expose = Mother_expose,
    Offspring_expose = Offspring_expose,
    Father_outcome = Father_outcome,
    Mother_outcome = Mother_outcome,
    Offspring_outcome = Offspring_outcome
    # 注意：这里可能还需要包含祖父母辈的SNP和暴露/结局数据，如果下游分析需要的话。
    # 当前版本主要集中在父母和子代。
  )

  # --- 14. 根据重叠比例为双样本MR抽取暴露组和结局组的样本索引 ---
  # 计算实际重叠的样本数量
  N_overlap_actual <- floor(min(N_exp, N_out) * overlap_prop)

  # 计算各组独有的样本数量
  N_exp_only_final <- N_exp - N_overlap_actual
  N_out_only_final <- N_out - N_overlap_actual

  # 从总个体索引池 (all_indices) 中进行抽样
  indices_pool <- all_indices # 初始化索引池

  # 1. 抽取重叠部分的样本索引
  indices_overlap <- sample(indices_pool, size = N_overlap_actual, replace = FALSE)
  indices_pool <- setdiff(indices_pool, indices_overlap) # 从池中移除已抽取的重叠样本

  # 2. 抽取暴露组独有的样本索引
  indices_exp_only <- sample(indices_pool, size = N_exp_only_final, replace = FALSE)
  indices_pool <- setdiff(indices_pool, indices_exp_only) # 从池中移除

  # 3. 抽取结局组独有的样本索引
  indices_out_only <- sample(indices_pool, size = N_out_only_final, replace = FALSE)
  # indices_pool <- setdiff(indices_pool, indices_out_only) # 池中剩余的未使用，如果需要可以记录

  # 组合最终的暴露组和结局组样本索引
  indices_exp_final <- c(indices_overlap, indices_exp_only)
  indices_out_final <- c(indices_overlap, indices_out_only)

  # --- 15. 根据抽取的索引创建最终的暴露组和结局组数据框 ---
  # 暴露组数据 (包含ID以及所有三代的SNP和暴露数据)
  exposure_data <- data_all %>%
    filter(ID %in% indices_exp_final) %>%
    # 选择研究所需的列，这里包含了所有三代的信息以供灵活分析
    dplyr::select(
      ID,
      Father_SNPs, Mother_SNPs, Offspring_SNPs,
      Father_expose, Mother_expose, Offspring_expose
      # 如果需要结局信息在暴露数据中（例如用于真实效应估计），可以添加
      # Father_outcome, Mother_outcome, Offspring_outcome
    )

  # 结局组数据 (包含ID以及所有三代的SNP和结局数据)
  outcome_data <- data_all %>%
    filter(ID %in% indices_out_final) %>%
    dplyr::select(
      ID,
      Father_SNPs, Mother_SNPs, Offspring_SNPs,
      # Father_expose, Mother_expose, Offspring_expose, # 通常结局GWAS不直接测量暴露
      Father_outcome, Mother_outcome, Offspring_outcome
    )

  # --- 16. 返回包含暴露数据和结局数据的列表 ---
  return(list(exposure_data = as.data.frame(exposure_data), outcome_data = as.data.frame(outcome_data)))
}

#' @title 生成多个模拟数据集用于孟德尔随机化分析 (版本3)
#' @description 该函数能够生成多个数据集，每个数据集模拟了包含潜在基因多效性、
#'              选型婚配和人群分层等复杂情况的三代家系遗传数据。
#'              它是对 generate_mr_trio_data_2sample 函数的封装调用。
#'
#' @param n 整数型。要生成的总数据集数量。
#' @param num_pleiotropic 整数型。指定在 n 个数据集中，有多少个数据集的SNP具有*潜在的*非零结局效应（即基因多效性）。
#' @param N_exp 整数型。每个数据集中暴露GWAS的样本量。
#' @param N_out 整数型。每个数据集中结局GWAS的样本量。
#' @param overlap_prop 数值型。暴露和结局样本之间的重叠比例 (0 到 1)。
#' @param p_f 数值型。用于生成SNP的基础等位基因频率 (例如，父系来源的等位基因频率或一个基准频率)。
#' @param p_m 数值型。用于生成SNP的另一个基础等位基因频率 (例如，母系来源的等位基因频率或另一个基准频率)。
#'
#' @section 暴露效应参数:
#' @param beta_FStoOE_exp 数值型。当存在时，父亲SNP对子代暴露的效应大小。
#' @param beta_MStoOE_exp 数值型。当存在时，母亲SNP对子代暴露的效应大小。
#' @param beta_OStoOE_exp 数值型。当存在时，子代自身SNP对子代暴露的效应大小。
#'
#' @section 结局效应/水平多效性参数 (当SNP具有多效性时):
#' @param mean_beta_FStoOE_out 数值型。父亲SNP对子代结局的多效性效应大小的均值 (绝对值)。
#' @param sd_beta_FStoOE_out 数值型。父亲SNP对子代结局的多效性效应大小的标准差。
#' @param mean_beta_MStoOE_out 数值型。母亲SNP对子代结局的多效性效应大小的均值 (绝对值)。
#' @param sd_beta_MStoOE_out 数值型。母亲SNP对子代结局的多效性效应大小的标准差。
#' @param mean_beta_OStoOE_out 数值型。子代自身SNP对子代结局的多效性效应大小的均值 (绝对值)。
#' @param sd_beta_OStoOE_out 数值型。子代自身SNP对子代结局的多效性效应大小的标准差。
#' @param prop_negative_pleiotropy 数值型。在被指定为具有多效性的SNP中，其多效性效应为负向的比例 (0 到 1)。
#'
#' @section 选型婚配参数:

#' @section 人群分层参数:
#' @param crowd_stratification_differences 数值型。用于模拟人群分层的等位基因频率差异。
#'                                         例如，如果p_f是0.3，此参数是0.05，则实际使用的p_f可能是0.3或0.35。
#'
#' @section 其他效应参数:
#' @param beta_exp_to_out 数值型。暴露对结局的真实因果效应。
#' @param beta_confounding_exp 数值型。影响暴露的混杂因素的方差。
#' @param beta_confounding_out 数值型。影响结局的混杂因素的方差。
#' @param correlation 数值型。用于模拟共享环境的方差组成部分。
#' @param seed 整数型。用于随机数生成的种子，以确保结果可重复。
#'
#' @return 一个列表，包含两个合并后的数据框：
#'         \code{[[1]]} 是所有数据集中暴露组的信息 (\code{combined_exposure_data})。
#'         \code{[[2]]} 是所有数据集中结局组的信息 (\code{combined_outcome_data})。
#'         每个数据框都额外添加了 \code{dataset_id} 和 \code{snp_type} 列。
#' @md
generate_multiple_datasets_v3 <- function(
    n = 10, # 要生成的总数据集数量
    num_pleiotropic = 1, # 指定多少个数据集的SNP具有 *潜在的* 非零结局效应(即基因多效性)
    N_exp = 1000, N_out = 1000, overlap_prop = 0,
    p_f = 0.3, p_m = 0.3,
    # --- 暴露效应 (当非零时的大小) ---
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3,
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
    beta_exp_to_out = 0.4, # 暴露对结局的真实因果效应
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


  # --- 3. 循环生成 n 个独立的数据集 ---
  all_datasets <- vector("list", n) # 初始化一个列表以存储所有生成的数据集

  for (i in 1:n) { # 对每个数据集进行循环
    # --- 3a. 判断当前数据集 i 的特性 ---
    # 检查当前数据集的索引 i 是否在之前抽取的 "pleiotropic_indices" 中
    is_pleiotropic_candidate <- i %in% pleiotropic_indices # 如果为 TRUE, 则此数据集的SNP是潜在多效性的

    # --- 3b. 为当前循环设置将要传递给底层数据生成函数的 beta 参数值 ---

    # --- 暴露效应参数 (在所有数据集中保持不变) ---
    current_beta_FStoOE_exp <- beta_FStoOE_exp
    current_beta_MStoOE_exp <- beta_MStoOE_exp
    current_beta_OStoOE_exp <- beta_OStoOE_exp


    # --- 结局效应 / 水平多效性参数 (核心修改部分) ---
    if (is_pleiotropic_candidate) {
      # 如果当前SNP是潜在多效性的:
      # 1. 决定这个多效性SNP的效应方向 (正或负)
      #    rbinom(1, 1, prob) 返回0或1。这里用 runif(1) < prop_negative_pleiotropy 判断是否为负向。
      direction_multiplier <- if (runif(1) < prop_negative_pleiotropy) -1 else 1 # -1代表负向, 1代表正向

      # 2. 为这个特定的SNP i (数据集 i) 抽样生成多效性效应的大小
      #    效应大小从一个正态分布中抽样，该分布的均值为 (方向 * 输入的均值参数)，标准差为输入的标准差参数。
      current_beta_FStoOE_out <- rnorm(1, mean = direction_multiplier * mean_beta_FStoOE_out, sd = sd_beta_FStoOE_out)
      current_beta_MStoOE_out <- rnorm(1, mean = direction_multiplier * mean_beta_MStoOE_out, sd = sd_beta_MStoOE_out)
      current_beta_OStoOE_out <- rnorm(1, mean = direction_multiplier * mean_beta_OStoOE_out, sd = sd_beta_OStoOE_out)

      snp_label <- "pleiotropic_variable" # 标记此SNP为具有（随机生成的）多效性效应
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
    is_compatibility_selection <- if (runif(1) < assortative_mating_prob) 1 else 0
    if (is_compatibility_selection == 1) {
      # 如果应用选型婚配，则使用函数输入的选型婚配参数
      current_assortative_mating_strength <- assortative_mating_strength
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
    # 将当前循环中确定的参数传递给底层函数
    # 注意: 底层函数 generate_mr_trio_data_2sample 需要能够接收并处理这里传递的所有参数，
    data_i <- generate_mr_trio_data_2sample(
      N_exp = N_exp, N_out = N_out,
      overlap_prop = overlap_prop,
      p_f = current_p_f, p_m = current_p_m, # 使用当前循环的等位基因频率
      beta_FStoOE_exp = current_beta_FStoOE_exp,
      beta_MStoOE_exp = current_beta_MStoOE_exp,
      beta_OStoOE_exp = current_beta_OStoOE_exp,
      # 传递当前循环抽样生成的或设为0的结局效应值
      beta_FStoOE_out = current_beta_FStoOE_out,
      beta_MStoOE_out = current_beta_MStoOE_out,
      beta_OStoOE_out = current_beta_OStoOE_out,
      beta_exp_to_out = beta_exp_to_out,
      # beta_confounding_exp 和 beta_confounding_out 可能是笔误，应为 confounding_exp, confounding_out (方差参数)
      # 或者底层函数期望的是效应大小，这里传递的是方差，需确认。假设底层函数处理这些。
      confounding_exp = beta_confounding_exp, # 传递混杂方差 (或效应，取决于底层函数定义)
      confounding_out = beta_confounding_out, # 传递混杂方差 (或效应，取决于底层函数定义)
      correlation = correlation, seed = NULL, # 在循环内部不应重设主种子，除非有意为每个数据集设独立子种子
      assortative_mating_strength = current_assortative_mating_strength
    )

    # --- 3d. 为生成的数据集添加唯一的 dataset_id 和 SNP 类型标签 ---
    # 检查底层函数返回结果的格式是否符合预期 (一个包含至少两个数据框的列表)
    if (length(data_i) >= 2 && is.data.frame(data_i[[1]]) && is.data.frame(data_i[[2]])) {
      data_i[[1]]$dataset_id <- i # 为暴露数据框添加 dataset_id
      data_i[[2]]$dataset_id <- i # 为结局数据框添加 dataset_id
      # 添加 SNP 类型标签 (例如 "pleiotropic_variable" 或 "valid_instrument")
      data_i[[1]]$snp_type <- snp_label
      data_i[[2]]$snp_type <- snp_label
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


  # --- [可选步骤] 重命名结局数据框中的列名 ---
  # 这个步骤的目的是，如果结局数据框中的列名与暴露数据框中的某些列名因来源（如父母、子代）而相似，
  # 仅以 "_outcome" 后缀区分，这里尝试将其重命名，例如，统一为类似暴露的命名方式以便后续处理。
  # 但这种操作具有特定性，需要确保符合实际需求。
  if (!is.null(combined_outcome_data) && nrow(combined_outcome_data) > 0) {
    current_names <- names(combined_outcome_data)
    # 尝试将以 "_outcome" 结尾的列名替换为以 "_expose" 结尾
    # 例如 "Father_outcome" -> "Father_expose"
    # 注意: 这个替换逻辑可能需要根据底层函数 generate_mr_trio_data_2sample 的实际输出列名进行调整。
    # 原注释提示: # 注意: 这里可能需要调整以匹配你的底层函数输出
    # 更安全的做法是只重命名确实存在的、且符合模式的列
    cols_to_rename <- grep("_outcome$", current_names, value = TRUE)
    if (length(cols_to_rename) > 0) {
      new_col_names <- sub("_outcome$", "_expose", cols_to_rename)
      names(combined_outcome_data)[match(cols_to_rename, names(combined_outcome_data))] <- new_col_names
    }
  }


  # --- 5. 返回包含合并后的暴露和结局数据的列表 ---
  return(list(exposure_data = combined_exposure_data, outcome_data = combined_outcome_data))
}

# %% 新升级的数据生成函数
#' @param n 代表整数
#' @param p 代表率
#' @param beta 代表回归系数
#' @param var 代表方差
#' @param r 代表相关系数
#' @param bool 代表逻辑
generate_mr_trio_data_ultra <- function(
    n_snps = 3,
    n_independent = 500, n_trio = 500, p_overlap = 0,
    p_f = 0.3, p_m = 0.3, # p_m 当前未在SNP生成中使用
    # 暴露效应
    beta_fs_oe_exp = 0.3, beta_ms_oe_exp = 0.3,
    beta_os_oe_exp = 0.3,
    # 结局效应 (直接多效性 / 遗传叠加效应)
    beta_fs_oe_out = 0, beta_ms_oe_out = 0,
    beta_os_oe_out = 0,
    # 因果效应
    beta_exp_to_out = 0.4,
    # 混杂效应
    var_confounding_exp = 0.2, var_confounding_out = 0.2,
    # 其他参数
    r_correlation = 0.2, n_seed = NULL,
    # 选型婚配强度(跨性状)
    n_am = 1000) {
  # --- 0. 设置随机种子 (确保结果可重复性) ---
  set.seed(seed)

  # --- 系数弄成矩阵
  matrix_beta_os_oe_exp <- matrix(beta_os_oe_exp / n_snps,
    ncol = 1, nrow = n_snps
  )
  matrix_beta_fs_oe_exp <- matrix(beta_fs_oe_exp / n_snps,
    ncol = 1, nrow = n_snps
  )
  matrix_beta_ms_oe_exp <- matrix(beta_ms_oe_exp / n_snps,
    ncol = 1, nrow = n_snps
  )
  matrix_beta_os_oe_out <- matrix(beta_os_oe_out / n_snps,
    ncol = 1, nrow = n_snps
  )
  matrix_beta_fs_oe_out <- matrix(beta_fs_oe_out / n_snps,
    ncol = 1, nrow = n_snps
  )
  matrix_beta_ms_oe_out <- matrix(beta_ms_oe_out / n_snps,
    ncol = 1, nrow = n_snps
  )
  # ---


  #' @title (内部辅助函数) 模拟从单个亲本传递给子代的等位基因
  #' @description 根据亲本的SNP基因型，模拟孟德尔遗传中传递给子代的等位基因。
  #' @param parent_snps 数值向量。亲本的基因型 (0, 1, 或 2)。
  #'                    编码: 0=aa, 1=Aa, 2=AA (a是等位基因0, A是等位基因1)。
  #' @return 一个与 parent_snps 等长的数值向量，包含传递的等位基因 (0 或 1)。
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

  # 生成个体数据的
  ## 生成爷爷奶奶外公外婆的数据
  grandfather_father_snps <- matrix(0, nrow = n_independent, ncol = n_snps)
  grandmother_father_snps <- matrix(0, nrow = n_independent, ncol = n_snps)
  grandfather_mother_snps <- matrix(0, nrow = n_independent, ncol = n_snps)
  grandmother_mother_snps <- matrix(0, nrow = n_independent, ncol = n_snps)
  for (i in 1:n_snps) {
    grandfather_father_snps[, i] <- generate_snp_hwe(n_samples = n_independent, maf = p_f)
    grandmother_father_snps[, i] <- generate_snp_hwe(n_samples = n_independent, maf = p_f)
    grandfather_mother_snps[, i] <- generate_snp_hwe(n_samples = n_independent, maf = p_f)
    grandmother_mother_snps[, i] <- generate_snp_hwe(n_samples = n_independent, maf = p_f)
  }



  # --- 2. 定义暴露和结局的生成函数 ---
  generate_expose_function <- function(matrix_beta_os_oe_exp,
                                       matrix_beta_fs_oe_exp,
                                       matrix_beta_ms_oe_exp,
                                       snps, snps_father, snps_mother,
                                       values_confounder,
                                       values_correlation_factor) {
    values_expose <- snps %*% matrix_beta_os_oe_exp +
      snps_father %*% matrix_beta_fs_oe_exp +
      snps_mother %*% matrix_beta_ms_oe_exp +
      values_confounder + values_correlation_factor

    # 计算确定性部分的方差
    var_results <- var(values_expose)

    # 计算随机误差的标准差，以使最终暴露的方差为1
    # 如果 var_results >= 1, 则误差标准差为0
    sd_err <- sqrt(1 - var_results)

    # 添加随机误差，生成最终的暴露值
    value_expose <- value_expose + rnorm(n_independent, mean = 0, sd = sd_err)
    return(value_expose)
  }

  generate_outcome_function <- function(beta_exp_to_out,
                                        matrix_beta_os_oe_out,
                                        matrix_beta_fs_oe_out,
                                        matrix_beta_ms_oe_out,
                                        value_expose, snps, snps_father, snps_mother,
                                        values_confounder, values_correlation_factor) {
    # 计算结局的确定性部分
    outcome_deterministic <- beta_exp_to_out * value_expose +
      snps %*% matrix_beta_os_oe_out +
      snps_father %*% matrix_beta_fs_oe_out +
      snps_mother %*% matrix_beta_ms_oe_out +
      values_confounder + values_correlation_factor

    # 计算确定性部分的方差
    var_results <- var(outcome_deterministic, na.rm = TRUE)
    # 计算随机误差的标准差，以使最终结局的方差为1
    sd_err <- sqrt(max(0, 1 - var_results))
    # 添加随机误差，生成最终的结局值
    outcome_all <- outcome_deterministic + rnorm(n_independent, mean = 0, sd = sd_err) # 注意：rnorm的N应与SNPs_self长度一致
    return(outcome_all)
  }

  # --- 3. 为所有 N_total 个体生成混杂因素和共享环境因素 ---
  # 这些因素将在后续生成暴露和结局时使用
  # 影响暴露的混杂因素 (均值为0, 方差由 confounding_exp 定义)
  values_confounder_exp <- rnorm(n_independent,
    mean = 0,
    sd = sqrt(var_confounding_exp)
  )
  values_confounder_out <- rnorm(n_independent,
    mean = 0,
    sd = sqrt(var_confounding_exp)
  )

  values_correlation_factor_exp <- rnorm(n_independent,
    mean = 0,
    sd = sqrt(r_correlation)
  ) # 用于暴露
  values_correlation_factor_out <- rnorm(n_independent,
    mean = 0,
    sd = sqrt(r_correlation)
  ) # 用于结局


  # 生成暴露
  grandfather_father_expose <- generate_expose_function(
    matrix_beta_os_oe_exp, matrix(0, nrow = n_snps, ncol = 1), matrix(0, nrow = n_snps, ncol = 1), # 自身SNP效应, 无父母SNP效应
    grandfather_father_snps, grandfather_father_snps, grandfather_father_snps, # 自身SNP作为父母SNP传入
    values_confounder_exp, values_correlation_factor_exp
  )
  grandmother_father_expose <- generate_expose_function(
    matrix_beta_os_oe_exp, matrix(0, nrow = n_snps, ncol = 1), matrix(0, nrow = n_snps, ncol = 1), # 自身SNP效应, 无父母SNP效应
    grandmother_father_snps, grandfather_father_snps, grandfather_father_snps, # 自身SNP作为父母SNP传入
    values_confounder_exp, values_correlation_factor_exp
  )
  grandfather_mother_expose <- generate_expose_function(
    matrix_beta_os_oe_exp, matrix(0, nrow = n_snps, ncol = 1), matrix(0, nrow = n_snps, ncol = 1), # 自身SNP效应, 无父母SNP效应
    grandfather_mother_snps, grandfather_father_snps, grandfather_father_snps, # 自身SNP作为父母SNP传入
    values_confounder_exp, values_correlation_factor_exp
  )
  grandmother_mother_expose <- generate_expose_function(
    matrix_beta_os_oe_exp, matrix(0, nrow = n_snps, ncol = 1), matrix(0, nrow = n_snps, ncol = 1), # 自身SNP效应, 无父母SNP效应
    grandmother_mother_snps, grandfather_father_snps, grandfather_father_snps, # 自身SNP作为父母SNP传入
    values_confounder_exp, values_correlation_factor_exp
  )


  # 生成结局
  grandfather_father_outcome <- generate_outcome_function(
    beta_exp_to_out, matrix_beta_os_oe_out,
    matrix(0, nrow = n_snps, ncol = 1), matrix(0, nrow = n_snps, ncol = 1), # 自身SNP效应, 无父母SNP效应
    grandfather_father_expose, grandfather_father_snps,
    grandfather_father_snps, grandfather_father_snps, # 自身SNP作为父母SNP传入
    values_confounder_out, values_correlation_factor_out
  )
  grandmother_father_outcome <- generate_outcome_function(
    beta_exp_to_out, matrix_beta_os_oe_out,
    matrix(0, nrow = n_snps, ncol = 1), matrix(0, nrow = n_snps, ncol = 1), # 自身SNP效应, 无父母SNP效应
    grandmother_father_expose, grandmother_father_snps,
    grandfather_father_snps, grandfather_father_snps, # 自身SNP作为父母SNP传入
    values_confounder_out, values_correlation_factor_out
  )
  grandfather_mother_outcome <- generate_outcome_function(
    beta_exp_to_out, matrix_beta_os_oe_out,
    matrix(0, nrow = n_snps, ncol = 1), matrix(0, nrow = n_snps, ncol = 1), # 自身SNP效应, 无父母SNP效应
    grandfather_mother_expose, grandfather_mother_snps,
    grandfather_father_snps, grandfather_father_snps, # 自身SNP作为父母SNP传入
    values_confounder_out, values_correlation_factor_out
  )
  grandmother_mother_outcome <- generate_outcome_function(
    beta_exp_to_out, matrix_beta_os_oe_out,
    matrix(0, nrow = n_snps, ncol = 1), matrix(0, nrow = n_snps, ncol = 1), # 自身SNP效应, 无父母SNP效应
    grandmother_mother_expose, grandmother_mother_snps,
    grandfather_father_snps, grandfather_father_snps, # 自身SNP作为父母SNP传入
    values_confounder_out, values_correlation_factor_out
  )


  assortative_mating_function <- function(snps_1, expose_1, outcome_1,
                                          snps_2, expose_2, outcome_2,
                                          am_strength) {
    # 创建包含SNP和暴露的数据框
    object_1 <- data.frame(
      snps = snps_1,
      expose = expose_1, outcome = outcome_1
    )
    object_2 <- data.frame(
      snps = snps_2,
      expose = expose_2, outcome = outcome_2
    )

    # 生成用于选型婚配的暴露代理表型 (真实暴露 + 随机噪音)
    # 噪音的标准差是 sqrt(am_strength)
    object_1 <- object_1 %>% mutate(expose_star = expose +
      rnorm(n(), mean = 0, sd = sqrt(am_strength)))
    object_2 <- object_2 %>% mutate(outcome_star = outcome +
      rnorm(n(), mean = 0, sd = sqrt(am_strength)))

    # 根据暴露代理表型对两组个体分别进行升序排序
    object_1 <- arrange(object_1, expose_star)
    object_2 <- arrange(object_2, outcome_star)
    # 排序后，object_1的第一行将与object_2的第一行配对，以此类推。
    return(list(object_1, object_2))
  }

  # --- 6. 对祖父母辈进行选型婚配 ---
  # 父系祖父母 (爷爷奶奶) 进行选型婚配
  grand_father_am_results <- assortative_mating_function(
    grandfather_father_snps, grandfather_father_expose,
    grandfather_father_outcome, # 爷爷的信息
    grandmother_father_snps, grandmother_father_expose,
    grandmother_father_outcome, # 奶奶的信息
    assortative_mating_strength
  )
  # 母系祖父母 (外公外婆) 进行选型婚配
  grand_mother_am_results <- assortative_mating_function(
    grandfather_mother_snps, grandfather_mother_expose,
    grandfather_mother_outcome, # 外公的信息
    grandmother_mother_snps, grandmother_mother_expose,
    grandmother_mother_outcome, # 外婆的信息
    assortative_mating_strength
  )

  # 更新祖父母的SNP向量，使其根据选型婚配的结果重新排序
  # 这样，后续从这些排序后的祖父母传递等位基因时，就体现了选型婚配的效果
  grandfather_father_snps <- as.matrix(grand_father_am_results[[1]][, 1:n_snps]) # 排序后的爷爷SNP
  grandmother_father_snps <- as.matrix(grand_father_am_results[[2]][, 1:n_snps]) # 排序后并与爷爷配对的奶奶SNP

  grandfather_mother_snps <- as.matrix(grand_mother_am_results[[1]][, 1:n_snps]) # 排序后的外公SNP
  grandmother_mother_snps <- as.matrix(grand_mother_am_results[[2]][, 1:n_snps]) # 排序后并与外公配对的外婆SNP

  # --- 7. 生成父母代基因型 (通过模拟祖父母向父母传递等位基因) ---
  # 父亲的基因型 = 爷爷传递的等位基因 + 奶奶传递的等位基因
  father_snps <- matrix(0, nrow = n_independent, ncol = n_snps)
  mother_snps <- matrix(0, nrow = n_independent, ncol = n_snps)
  for (i in 1:n_snps) {
    father_snps[, i] <- internal_get_transmitted_allele(grandfather_father_snps[, i]) +
      internal_get_transmitted_allele(grandmother_father_snps[, i])
    # 母亲的基因型 = 外公传递的等位基因 + 外婆传递的等位基因
    mother_snps[, i] <- internal_get_transmitted_allele(grandfather_mother_snps[, i]) +
      internal_get_transmitted_allele(grandmother_mother_snps[, i])
  }


  # --- 8. 为父母代生成暴露数据 (用于他们之间的选型婚配) ---


  father_expose <- generate_expose_function(
    matrix_beta_os_oe_exp, matrix_beta_fs_oe_exp, matrix_beta_ms_oe_exp, # 自身SNP效应, 无父母SNP效应
    father_snps, grandfather_father_snps, grandmother_father_snps, # 自身SNP作为父母SNP传入
    values_confounder_exp, values_correlation_factor_exp
  )
  mother_expose <- generate_expose_function(
    matrix_beta_os_oe_exp, matrix_beta_fs_oe_exp, matrix_beta_ms_oe_exp, # 自身SNP效应, 无父母SNP效应
    mother_snps, grandfather_mother_snps, grandmother_mother_snps, # 自身SNP作为父母SNP传入
    values_confounder_exp, values_correlation_factor_exp
  )
  father_outcome <- generate_outcome_function(
    beta_exp_to_out, matrix_beta_os_oe_out,
    matrix_beta_fs_oe_out, matrix_beta_ms_oe_out, # 自身SNP效应, 无父母SNP效应
    father_expose, father_snps,
    grandfather_father_snps, grandmother_father_snps, # 自身SNP作为父母SNP传入
    values_confounder_out, values_correlation_factor_out
  )
  mother_outcome <- generate_outcome_function(
    beta_exp_to_out, matrix_beta_os_oe_out,
    matrix_beta_fs_oe_out, matrix_beta_ms_oe_out, # 自身SNP效应, 无父母SNP效应
    mother_expose, mother_snps,
    grandfather_mother_snps, grandmother_mother_snps, # 自身SNP作为父母SNP传入
    values_confounder_out, values_correlation_factor_out
  )


  # --- 9. 对父母代进行选型婚配 ---
  parent_am_results <- assortative_mating_function(
    father_snps, father_expose, father_outcome, # 父亲的信息
    mother_snps, mother_expose, mother_outcome, # 母亲的信息
    assortative_mating_strength
  )
  # 更新父母的SNP和暴露向量，使其根据选型婚配的结果重新排序
  Father_SNPs <- parent_am_results[[1]]$snps # 排序后的父亲SNP
  Mother_SNPs <- parent_am_results[[2]]$snps # 排序后并与父亲配对的母亲SNP
  Father_expose <- parent_am_results[[1]]$expose # 排序后的父亲暴露
  Mother_expose <- parent_am_results[[2]]$expose # 排序后并与父亲配对的母亲暴露
  Father_outcome <- parent_am_results[[1]]$outcome # 排序后的父亲暴露
  Mother_outcome <- parent_am_results[[2]]$outcome # 排序后并与父亲配对的母亲暴露

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



  # --- 12. 为父亲、母亲和子代生成结局数据 ---
  # 父亲的结局
  # 子代的结局
  Offspring_outcome <- generate_outcome_function(
    beta_exp_to_out, beta_OStoOE_out, beta_FStoOE_out, beta_MStoOE_out,
    Offspring_expose, Offspring_SNPs, Father_SNPs, Mother_SNPs, # 子代的父母是Father和Mother
    confounder_out_values, correlation_factor_out_values
  )

  # --- 13. 创建包含所有模拟信息的完整数据集 ---
  # 这个数据集包含了 N_total 个家庭（父亲-母亲-子代）的完整信息
  data_all <- data.frame(
    ID = all_indices, # 个体/家庭的唯一ID，方便后续追踪和抽样
    Father_SNPs = Father_SNPs,
    Mother_SNPs = Mother_SNPs,
    Offspring_SNPs = Offspring_SNPs,
    Father_expose = Father_expose,
    Mother_expose = Mother_expose,
    Offspring_expose = Offspring_expose,
    Father_outcome = Father_outcome,
    Mother_outcome = Mother_outcome,
    Offspring_outcome = Offspring_outcome
    # 注意：这里可能还需要包含祖父母辈的SNP和暴露/结局数据，如果下游分析需要的话。
    # 当前版本主要集中在父母和子代。
  )

  # --- 14. 根据重叠比例为双样本MR抽取暴露组和结局组的样本索引 ---
  # 计算实际重叠的样本数量
  N_overlap_actual <- floor(min(N_exp, N_out) * overlap_prop)

  # 计算各组独有的样本数量
  N_exp_only_final <- N_exp - N_overlap_actual
  N_out_only_final <- N_out - N_overlap_actual

  # 从总个体索引池 (all_indices) 中进行抽样
  indices_pool <- all_indices # 初始化索引池

  # 1. 抽取重叠部分的样本索引
  indices_overlap <- sample(indices_pool, size = N_overlap_actual, replace = FALSE)
  indices_pool <- setdiff(indices_pool, indices_overlap) # 从池中移除已抽取的重叠样本

  # 2. 抽取暴露组独有的样本索引
  indices_exp_only <- sample(indices_pool, size = N_exp_only_final, replace = FALSE)
  indices_pool <- setdiff(indices_pool, indices_exp_only) # 从池中移除

  # 3. 抽取结局组独有的样本索引
  indices_out_only <- sample(indices_pool, size = N_out_only_final, replace = FALSE)
  # indices_pool <- setdiff(indices_pool, indices_out_only) # 池中剩余的未使用，如果需要可以记录

  # 组合最终的暴露组和结局组样本索引
  indices_exp_final <- c(indices_overlap, indices_exp_only)
  indices_out_final <- c(indices_overlap, indices_out_only)

  # --- 15. 根据抽取的索引创建最终的暴露组和结局组数据框 ---
  # 暴露组数据 (包含ID以及所有三代的SNP和暴露数据)
  exposure_data <- data_all %>%
    filter(ID %in% indices_exp_final) %>%
    # 选择研究所需的列，这里包含了所有三代的信息以供灵活分析
    dplyr::select(
      ID,
      Father_SNPs, Mother_SNPs, Offspring_SNPs,
      Father_expose, Mother_expose, Offspring_expose
      # 如果需要结局信息在暴露数据中（例如用于真实效应估计），可以添加
      # Father_outcome, Mother_outcome, Offspring_outcome
    )

  # 结局组数据 (包含ID以及所有三代的SNP和结局数据)
  outcome_data <- data_all %>%
    filter(ID %in% indices_out_final) %>%
    dplyr::select(
      ID,
      Father_SNPs, Mother_SNPs, Offspring_SNPs,
      # Father_expose, Mother_expose, Offspring_expose, # 通常结局GWAS不直接测量暴露
      Father_outcome, Mother_outcome, Offspring_outcome
    )

  # --- 16. 返回包含暴露数据和结局数据的列表 ---
  return(list(exposure_data = as.data.frame(exposure_data), outcome_data = as.data.frame(outcome_data)))
}
generate_multiple_datasets_v4_ultra <- function(
    n = 10, # 要生成的总数据集数量
    num_pleiotropic = 1, # 指定多少个数据集的SNP具有 *潜在的* 非零结局效应(即基因多效性)
    n_independent_individual = 5,
    N_exp = 1000, N_out = 1000, overlap_prop = 0,
    p_f = 0.3, p_m = 0.3,
    # --- 暴露效应 (当非零时的大小) ---
    beta_FStoOE_exp = 0.1, beta_MStoOE_exp = 0.1,
    beta_OStoOE_exp = 0.3,
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
    beta_exp_to_out = 0.4, # 暴露对结局的真实因果效应
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
  # 独立个体数据集指示变量
  independent_indices <- 1:n_independent_individual

  # --- 3. 循环生成 n 个独立的数据集 ---
  all_datasets <- vector("list", n) # 初始化一个列表以存储所有生成的数据集

  for (i in 1:n) { # 对每个数据集进行循环
    # --- 3a. 判断当前数据集 i 的特性 ---
    # 检查当前数据集的索引 i 是否在之前抽取的 "pleiotropic_indices" 中
    is_pleiotropic_candidate <- i %in% pleiotropic_indices # 如果为 TRUE, 则此数据集的SNP是潜在多效性的

    # --- 3b. 为当前循环设置将要传递给底层数据生成函数的 beta 参数值 ---

    # --- 暴露效应参数 (在所有数据集中保持不变) ---
    current_beta_FStoOE_exp <- beta_FStoOE_exp
    current_beta_MStoOE_exp <- beta_MStoOE_exp
    current_beta_OStoOE_exp <- beta_OStoOE_exp


    # --- 结局效应 / 水平多效性参数 (核心修改部分) ---
    if (is_pleiotropic_candidate) {
      # 如果当前SNP是潜在多效性的:
      # 1. 决定这个多效性SNP的效应方向 (正或负)
      #    rbinom(1, 1, prob) 返回0或1。这里用 runif(1) < prop_negative_pleiotropy 判断是否为负向。
      direction_multiplier <- if (runif(1) < prop_negative_pleiotropy) -1 else 1 # -1代表负向, 1代表正向

      # 2. 为这个特定的SNP i (数据集 i) 抽样生成多效性效应的大小
      #    效应大小从一个正态分布中抽样，该分布的均值为 (方向 * 输入的均值参数)，标准差为输入的标准差参数。
      current_beta_FStoOE_out <- rnorm(1, mean = direction_multiplier * mean_beta_FStoOE_out, sd = sd_beta_FStoOE_out)
      current_beta_MStoOE_out <- rnorm(1, mean = direction_multiplier * mean_beta_MStoOE_out, sd = sd_beta_MStoOE_out)
      current_beta_OStoOE_out <- rnorm(1, mean = direction_multiplier * mean_beta_OStoOE_out, sd = sd_beta_OStoOE_out)

      snp_label <- "pleiotropic_variable" # 标记此SNP为具有（随机生成的）多效性效应
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
    is_compatibility_selection <- if (runif(1) < assortative_mating_prob) 1 else 0
    if (is_compatibility_selection == 1) {
      # 如果应用选型婚配，则使用函数输入的选型婚配参数
      current_assortative_mating_strength <- assortative_mating_strength
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
    # 将当前循环中确定的参数传递给底层函数
    # 注意: 底层函数 generate_mr_trio_data_2sample 需要能够接收并处理这里传递的所有参数，
    data_i <- generate_mr_trio_data_ultra(
      N_exp = N_exp, N_out = N_out,
      overlap_prop = overlap_prop,
      p_f = current_p_f, p_m = current_p_m, # 使用当前循环的等位基因频率
      beta_FStoOE_exp = current_beta_FStoOE_exp,
      beta_MStoOE_exp = current_beta_MStoOE_exp,
      beta_OStoOE_exp = current_beta_OStoOE_exp,
      # 传递当前循环抽样生成的或设为0的结局效应值
      beta_FStoOE_out = current_beta_FStoOE_out,
      beta_MStoOE_out = current_beta_MStoOE_out,
      beta_OStoOE_out = current_beta_OStoOE_out,
      beta_exp_to_out = beta_exp_to_out,
      # beta_confounding_exp 和 beta_confounding_out 可能是笔误，应为 confounding_exp, confounding_out (方差参数)
      # 或者底层函数期望的是效应大小，这里传递的是方差，需确认。假设底层函数处理这些。
      confounding_exp = beta_confounding_exp, # 传递混杂方差 (或效应，取决于底层函数定义)
      confounding_out = beta_confounding_out, # 传递混杂方差 (或效应，取决于底层函数定义)
      correlation = correlation, seed = NULL, # 在循环内部不应重设主种子，除非有意为每个数据集设独立子种子
      assortative_mating_strength = current_assortative_mating_strength
    )
    # --- 独立个体参数设置 ---
    if (i %in% independent_indices) {
      data_i[[1]] <- data_i[[1]] %>% dplyr::select(ID, Offspring_SNPs, Offspring_expose)
      data_i[[2]] <- data_i[[2]] %>% dplyr::select(ID, Offspring_SNPs, Offspring_outcome)
    }

    # --- 3d. 为生成的数据集添加唯一的 dataset_id 和 SNP 类型标签 ---
    # 检查底层函数返回结果的格式是否符合预期 (一个包含至少两个数据框的列表)
    if (length(data_i) >= 2 && is.data.frame(data_i[[1]]) && is.data.frame(data_i[[2]])) {
      data_i[[1]]$dataset_id <- i # 为暴露数据框添加 dataset_id
      data_i[[2]]$dataset_id <- i # 为结局数据框添加 dataset_id
      # 添加 SNP 类型标签 (例如 "pleiotropic_variable" 或 "valid_instrument")
      data_i[[1]]$snp_type <- snp_label
      data_i[[2]]$snp_type <- snp_label
    } else {
      warning(paste("数据集", i, "的生成结果格式不符合预期 (data_i 不是包含两个数据框的列表)，无法添加 ID 和标签。"))
    }

    all_datasets[[i]] <- data_i # 将生成的单个数据集 (data_i) 存储到总列表 all_datasets 中
  } # 结束 for 循环 (生成 n 个数据集)


  # --- 5. 返回包含合并后的暴露和结局数据的列表 ---
  return(all_datasets)
}
