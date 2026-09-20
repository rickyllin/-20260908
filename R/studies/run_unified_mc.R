###############################################################################
# 工程一：統一所有蒙地卡羅的模擬規模
#
#   所有實驗一律 reps = 100、固定亂數種子，取代先前 reps = 8~50 的試跑結果。
#   輸出覆寫 mc_results.csv / lambda_selection.csv / lambda_grid3.csv，
#   並另存一份 *_reps100.csv 以保留來源可追溯。
###############################################################################

source("R/studies/mc_exposure.R")
source("R/studies/lambda_selection.R")
source("R/studies/lambda3_and_resampling.R")

REPS <- 100
SEED <- 20260921
truth <- mc_truth("Female", 2001, 2024)
cat(sprintf("真值：%d 齡組 x %d 年，reps = %d，seed = %d\n",
            length(truth$a), length(truth$k), REPS, SEED))

t_all <- Sys.time()

## ---- (1) 暴露數 x 估計方法 ----
cat("\n[1/3] mc_run\n")
res <- mc_run(truth, reps = REPS, seed = SEED)
write.csv(res, "output/tables/mc_results.csv", row.names = FALSE)
write.csv(res, "output/tables/mc_results_reps100.csv", row.names = FALSE)
attr(res, "drift_true") <- (truth$k[length(truth$k)] - truth$k[1]) / (length(truth$k) - 1)
mc_plots(res, file = "output/figures/mc_plots.png")
cat(sprintf("  完成，累計 %.1f 分鐘\n", as.numeric(Sys.time()-t_all, units="mins")))

## ---- (3) 三方向 lambda 格點（先跑較快者）----
cat("\n[2/3] lambda_grid3\n")
g <- lambda_grid3(truth, 2e5, reps = REPS, seed = SEED, verbose = FALSE)
g <- g[order(g$sse_b), ]
write.csv(g, "output/tables/lambda_grid3.csv", row.names = FALSE)
print(head(g, 5), digits = 4, row.names = FALSE)
cat(sprintf("  完成，累計 %.1f 分鐘\n", as.numeric(Sys.time()-t_all, units="mins")))

## ---- (2) lambda 選法 vs oracle（最耗時）----
cat("\n[3/3] lambda_study\n")
st <- lambda_study(truth, reps = REPS, B = 15, seed = SEED, verbose = FALSE)
write.csv(st, "output/tables/lambda_selection.csv", row.names = FALSE)
write.csv(st, "output/tables/lambda_selection_reps100.csv", row.names = FALSE)
print(st, digits = 4, row.names = FALSE)

cat(sprintf("\n全部完成，共 %.1f 分鐘\n", as.numeric(Sys.time()-t_all, units="mins")))
