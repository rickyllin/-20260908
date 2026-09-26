###############################################################################
# 四張 EDA 診斷圖：從資料分布決定該用哪個估計量
#
# 對應 0922 指導的第 (1) 項——須先自資料的分布與偏誤特性出發，
# 方能說明為何需要所提方法。本腳本把「模擬才得出的門檻」翻譯成
# 「拿到資料即可計算的量」。
#
# 輸出：output/figures/figEDA_panel.png
#       output/tables/tableEDA_diagnostics.csv
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/fig_axis_utils.R")
source("R/core/eda_diagnostics.R")
source("R/core/fh_firth_kalman.R")

SEED <- 20260926
YEARS <- 2001:2024

dat  <- load_data(sex = "Female")
keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]

## 真值：以全國女性 2001-2024 的卜瓦松 LC 配適為真值，年齡權重取最後一年
fit   <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(fit$a, rep(1, ncol(D0))) + outer(fit$b, fit$k))
wage  <- E0[, ncol(E0)] / sum(E0[, ncol(E0)])

sim_one <- function(N, seed) {
  set.seed(seed)
  E <- outer(N * wage, rep(1, ncol(D0)))
  dimnames(E) <- dimnames(D0)
  D <- matrix(rpois(length(E), E * mtrue), nrow(E), dimnames = dimnames(E))
  list(D = D, E = E)
}

cat("=== 計算四項診斷（sigma 由 Firth 預備配適的 mu_hat 取得）===\n")
firth_fit <- function(D, E) lc_poisson_firth(D, E, firth = TRUE)

s20 <- sim_one(2e5, SEED)
s05 <- sim_one(5e4, SEED + 1)

dg_all <- lc_eda(D0,     E0,     fitter = firth_fit)
dg_20  <- lc_eda(s20$D,  s20$E,  fitter = firth_fit)
dg_05  <- lc_eda(s05$D,  s05$E,  fitter = firth_fit)

labs <- c(sprintf("全國女性（%.0f 萬）", sum(E0[, ncol(E0)]) / 1e4),
          "模擬 20 萬", "模擬 5 萬")

for (i in seq_along(labs)) {
  dg <- list(dg_all, dg_20, dg_05)[[i]]
  cat("\n---- ", labs[i], " ----\n", sep = "")
  print(dg$summary, row.names = FALSE)
  cat("建議：\n"); cat(paste0("  ", dg$advice, collapse = "\n"), "\n")
}

## ---- 匯出 ----
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)

tab <- do.call(rbind, lapply(seq_along(labs), function(i) {
  dg <- list(dg_all, dg_20, dg_05)[[i]]
  data.frame(資料 = labs[i], dg$summary, stringsAsFactors = FALSE)
}))
write.csv(tab, "output/tables/tableEDA_diagnostics.csv", row.names = FALSE)

png_cjk("output/figures/figEDA_panel.png", width = 2000, height = 560, res = 150)
plot_lc_eda(dg_all, dg_20, dg_05, labels = labs)
dev.off()

cat("\n已輸出 output/figures/figEDA_panel.png 與 output/tables/tableEDA_diagnostics.csv\n")
