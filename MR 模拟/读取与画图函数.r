library(data.table)
basic_simulation_result <- as.data.frame(fread("MR 模拟/MR 模拟结果/sim_IVA-0_propNeg-0_overlap-0_indirect-0_beta-0.csv"))
mean(basic_simulation_result$B_p_value < 0.05)
