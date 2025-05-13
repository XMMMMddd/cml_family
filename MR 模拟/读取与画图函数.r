library(data.table)
basic_simulation_result <- as.data.frame(fread("MR 模拟/MR 模拟结果/基础测试.csv"))
mean(basic_simulation_result$G_p_value < 0.05)
