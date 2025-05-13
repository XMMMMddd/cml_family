# %%
library(dplyr)
library(doParallel)
library(purrr)
source("GWAS 模拟/GWAS模拟（单次）.R")

# %%
# GWAS 多次模拟函数

gwas_multiple_simulation <- function(
    n_simulations, param) {
    # 检测可用核数
    num_cores <- detectCores() - 1

    # 设置集群
    cluster <- makeCluster(num_cores)
    clusterEvalQ(cluster, {
        source("GWAS 模拟/GWAS模拟（单次）.R")
    })


    # 并行化计算
    results <- parLapply(
        cl = cluster, # 使用命名参数 cl
        X = 1:n_simulations, # 迭代序列
        # --- 匿名函数仍然需要接收 sim_index ---
        fun = function(sim_index, sim_params_list) {
            # 第一个参数 sim_index 接收迭代序号，但后面不用它
            # 调用内部模拟函数
            gwas_simulation_once(
                # --- 从列表解包参数 ---
                n = sim_params_list$n,
                num_pleiotropic = sim_params_list$num_pleiotropic,
                N_out = sim_params_list$N_out,
                # 样本重叠
                overlap_prop = sim_params_list$overlap_prop,
                # 直接效应
                beta_FStoOE_exp = sim_params_list$beta_FStoOE_exp,
                beta_MStoOE_exp = sim_params_list$beta_MStoOE_exp,
                beta_OStoOE_exp = sim_params_list$beta_OStoOE_exp,
                # 水平多效性
                prop_negative_pleiotropy =
                    sim_params_list$prop_negative_pleiotropy, # 不满足inside假设
                mean_beta_FStoOE_out = sim_params_list$mean_beta_FStoOE_out,
                sd_beta_FStoOE_out = sim_params_list$sd_beta_FStoOE_out,
                mean_beta_MStoOE_out = sim_params_list$mean_beta_MStoOE_out,
                sd_beta_MStoOE_out = sim_params_list$sd_beta_MStoOE_out,
                mean_beta_OStoOE_out = sim_params_list$mean_beta_OStoOE_out,
                sd_beta_OStoOE_out = sim_params_list$sd_beta_OStoOE_out,
                # 选型婚配
                ## 选型婚配基因型
                compatibility_selection_prop =
                    sim_params_list$compatibility_selection_prop,
                compatibility_selection_geno =
                    sim_params_list$compatibility_selection_geno,
                correlation_param = sim_params_list$correlation_param,
                ## 选型婚配暴露相关
                compatibility_selection_factor_exp =
                    sim_params_list$compatibility_selection_factor_exp,
                compatibility_selection_factor_out =
                    sim_params_list$compatibility_selection_factor_out,
                # 人群分层（双人群差异）
                crowd_stratification_differences =
                    sim_params_list$crowd_stratification_differences,
                # 其他参数设置
                beta_exp_to_out = sim_params_list$beta_exp_to_out,
                beta_confounding_exp = sim_params_list$beta_confounding_exp,
                beta_confounding_out = sim_params_list$beta_confounding_out,
                correlation = sim_params_list$correlation,
                seed = NULL
            )
        },
        # --- 将 'param' 列表传递给匿名函数的 'sim_params_list' 参数 ---
        sim_params_list = param
    )

    # 关闭集群
    stopCluster(cluster)


    return(results)
}
# %%
# 小模拟函数
sim_params_list <- list( # --- 从列表解包参数 ---
    n = 10, num_pleiotropic = 0, N_out = 1000,
    # 样本重叠
    overlap_prop = 0,
    # 直接效应
    beta_FStoOE_exp = 0.3,
    beta_MStoOE_exp = 0.3,
    beta_OStoOE_exp = 0.3,
    # 水平多效性
    prop_negative_pleiotropy =
        0, # 不满足inside假设
    mean_beta_FStoOE_out = 0,
    sd_beta_FStoOE_out = 0.05,
    mean_beta_MStoOE_out = 0,
    sd_beta_MStoOE_out = 0.05,
    mean_beta_OStoOE_out = 0,
    sd_beta_OStoOE_out = 0.05,
    # 选型婚配
    ## 选型婚配基因型
    compatibility_selection_prop = 0,
    compatibility_selection_geno = "independent", correlation_param = 0.5,
    ## 选型婚配暴露相关
    compatibility_selection_factor_exp = 0,
    compatibility_selection_factor_out = 0,
    # 人群分层（双人群差异）
    crowd_stratification_differences = 0,
    # 其他参数设置
    beta_exp_to_out = 0,
    beta_confounding_exp = 0.2,
    beta_confounding_out = 0.2,
    correlation = 0.2,
    seed = NULL
)
small_simulation <- gwas_multiple_simulation(500, sim_params_list)
small_simulation_dataframe <- simulation_ListtoTibble(small_simulation)$tibble
small_simulation_dataframe <- small_simulation_dataframe %>%
    mutate(beta_exp_z = beta_exp / beta_se_exp)
small_simulation_dataframe <- small_simulation_dataframe %>%
    mutate(beta_out_z = beta_out / beta_se_out)
small_simulation_dataframe <- small_simulation_dataframe %>%
    mutate(beta_exp_pvalue = 2 * pnorm(abs(beta_exp_z), lower.tail = FALSE))
small_simulation_dataframe <- small_simulation_dataframe %>%
    mutate(beta_out_pvalue = 2 * pnorm(abs(beta_out_z), lower.tail = FALSE))
# %%
mean(small_simulation_dataframe$beta_out_pvalue < 0.05)
