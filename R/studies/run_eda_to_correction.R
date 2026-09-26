###############################################################################
# EDA 診斷 -> 估計量選擇 -> 偏誤改善：把三者接成一條鏈
#
# 回答的問題：0926 那批探索裡，哪些部分真的改善了小樣本 alpha / beta 的偏誤，
#             而且改善的位置可由 EDA 事前指出？
#
# 設計：以全國女性 2001-2024 的卜瓦松 LC 配適為真值，
#       D ~ Poisson(N * w_x * m_true)，reps = 100、固定種子。
#
#   表 1（逐年齡）：零格比例（EDA 診斷②）× 各法的 alpha 偏誤
#                   -> 檢驗「零格比例可否事前指出 alpha 偏誤的位置」
#   表 2（彙總）  ：五種估計量 × 四項指標 × 三個人口規模
#
# 輸出：output/tables/tableH_alpha_by_age.csv
#       output/tables/tableI_estimator_summary.csv
#       output/figures/figH_eda_vs_alpha_bias.png
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/rabbi_mazzuco_replication.R")
source("R/core/heteropca_lc.R")
source("R/core/fh_firth_kalman.R")

SEED  <- 20260926
REPS  <- 100
NS    <- c(1e4, 5e4, 2e5)
YEARS <- 2001:2024

dat  <- load_data(sex = "Female")
keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A  <- nrow(D0); Tn <- ncol(D0); ages <- rownames(D0)

truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn])
drift_true <- (truth$k[Tn] - truth$k[1]) / (Tn - 1)
cat(sprintf("真值：drift = %.4f，beta 範圍 %.4f - %.4f\n",
            drift_true, min(truth$b), max(truth$b)))

## ---- 五種估計量 -----------------------------------------------------------
##  A 標準 LC（現行作法）      等權重 SVD、零格以 0.5 替代
##  B 加權 SVD（w = mu_hat）   Wilmoth (1993) 的修正版：權重用配適值而非觀測 D
##  C 卜瓦松 MLE               Brouhns et al. (2002)
##  D Firth                    Firth (1993) 偏誤減少
##  E Firth + FH               再加 Fay-Herriot 經驗貝氏收縮（僅動 beta）
ESTS <- c("A 標準 LC", "B 加權SVD(w=mu)", "C 卜瓦松 MLE", "D Firth", "E Firth+FH")

fit_one <- function(D, E, est) {
  switch(est,
    "A 標準 LC"        = lc_svd_fit(log(pmax(D, 0.5) / E)),
    "B 加權SVD(w=mu)"  = lc_wsvd(D, E, wmode = "mu"),
    "C 卜瓦松 MLE"     = lc_poisson_firth(D, E, firth = FALSE),
    "D Firth"          = lc_poisson_firth(D, E, firth = TRUE),
    "E Firth+FH"       = lc_fh(D, E, fitter = function(D, E) lc_poisson_firth(D, E, firth = TRUE)))
}

## 發散判準與真值的中心化相關
ok_fit   <- function(f) all(is.finite(f$a)) && all(is.finite(f$b)) &&
                        all(is.finite(f$k)) && max(abs(f$k)) < 1e3
cor_b    <- function(b) { if (sd(b) < 1e-12) NA_real_ else cor(b, truth$b) }

res_age <- list(); res_sum <- list()

for (N in NS) {
  E <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu_exp <- E * mtrue

  ## 儲存：逐重複 × 逐年齡的 alpha；逐重複的彙總量
  aM <- array(NA_real_, c(REPS, A, length(ESTS)))
  bS <- dS <- cS <- matrix(NA_real_, REPS, length(ESTS))
  div <- integer(length(ESTS)); zero_ct <- numeric(REPS)

  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000 * which(NS == N) + r)
    D <- matrix(rpois(length(E), mu_exp), A, dimnames = dimnames(E))
    zero_ct[r] <- sum(D == 0)
    if (r == 1) zero_by_age <- rowMeans(D == 0) else
      zero_by_age <- zero_by_age + rowMeans(D == 0)

    for (j in seq_along(ESTS)) {
      f <- tryCatch(fit_one(D, E, ESTS[j]), error = function(e) NULL)
      if (is.null(f) || !ok_fit(f)) { div[j] <- div[j] + 1L; next }
      aM[r, , j] <- f$a - truth$a
      bS[r, j]   <- sum((f$b - truth$b)^2) / sum(truth$b^2)   # 相對 SSE
      cS[r, j]   <- cor_b(f$b)
      dS[r, j]   <- (f$k[Tn] - f$k[1]) / (Tn - 1) - drift_true
    }
    if (r %% 20 == 0) cat(sprintf("  N=%.0e  rep %d/%d\n", N, r, REPS))
  }
  zero_by_age <- zero_by_age / REPS

  ## ---- 表 1：逐年齡 ----
  res_age[[as.character(N)]] <- data.frame(
    N = N, 年齡 = ages,
    期望死亡中位 = round(apply(mu_exp, 1, median), 1),
    零格比例 = round(zero_by_age, 3),
    setNames(as.data.frame(lapply(seq_along(ESTS), function(j)
      round(apply(aM[, , j], 2, median, na.rm = TRUE), 4))), ESTS),
    check.names = FALSE, stringsAsFactors = FALSE)

  ## ---- 表 2：彙總（穩健統計；alpha 取逐年齡偏誤絕對值的中位數）----
  res_sum[[as.character(N)]] <- data.frame(
    N = N, 估計量 = ESTS,
    平均零格 = round(mean(zero_ct), 1),
    alpha偏誤中位 = round(sapply(seq_along(ESTS), function(j)
      median(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    alpha偏誤最大 = round(sapply(seq_along(ESTS), function(j)
      max(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    SSE_beta = round(apply(bS, 2, median, na.rm = TRUE), 4),
    cor_beta = round(apply(cS, 2, median, na.rm = TRUE), 3),
    漂移偏誤 = round(apply(dS, 2, median, na.rm = TRUE), 4),
    發散率 = round(div / REPS, 3),
    stringsAsFactors = FALSE)

  cat(sprintf("\n===== N = %.0e（平均零格 %.1f / %d）=====\n",
              N, mean(zero_ct), A * Tn))
  print(res_sum[[as.character(N)]][, -1], row.names = FALSE)
}

tabH <- do.call(rbind, res_age); tabI <- do.call(rbind, res_sum)
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
write.csv(tabH, "output/tables/tableH_alpha_by_age.csv", row.names = FALSE)
write.csv(tabI, "output/tables/tableI_estimator_summary.csv", row.names = FALSE)

## ---- 診斷②（零格比例）與 alpha 偏誤的關係 ----
cat("\n=== EDA 診斷② 是否能事前指出 alpha 偏誤的位置 ===\n")
for (N in NS) {
  d <- res_age[[as.character(N)]]
  cat(sprintf("N=%.0e：零格比例 vs |標準 LC 的 alpha 偏誤| 的 Spearman 相關 = %.3f",
              N, cor(d$零格比例, abs(d[["A 標準 LC"]]), method = "spearman")))
  cat(sprintf("；Firth 為 %.3f\n",
              cor(d$零格比例, abs(d[["D Firth"]]), method = "spearman")))
}

## ---- 圖：零格比例（EDA）與各法 alpha 偏誤並排 ----
png("output/figures/figH_eda_vs_alpha_bias.png", width = 1900, height = 700, res = 150)
op <- par(mfrow = c(1, 3), mar = c(4.4, 4.6, 3.4, 1.2), bg = "#fcfcfb",
          col.axis = "#52514e", col.lab = "#52514e", fg = "#d6d5d0",
          cex.main = 1.05, font.main = 1)
n_ <- c(1, 4, rep(5, A - 2)); amid <- c(0, cumsum(n_)[-A]) + c(0.5, 2, rep(2.5, A - 2))
pal <- c("#c2410c", "#0d366b", "#2a78d6")
for (N in NS) {
  d <- res_age[[as.character(N)]]
  plot(NA, xlim = range(amid), ylim = range(unlist(d[, ESTS[c(1, 3, 4)]]), 0),
       xlab = "年齡", ylab = expression(alpha ~ "的偏誤（中位數）"),
       main = sprintf("N = %s（平均零格 %.0f/%d）",
                      format(N, big.mark = ",", scientific = FALSE),
                      res_sum[[as.character(N)]]$平均零格[1], A * Tn),
       col.main = "#0b0b0b")
  grid(NA, NULL, col = "#e8e7e2", lty = 1, lwd = .7); abline(h = 0, col = "#b8b7b2")
  ## 背景：零格比例（右軸概念，以灰色柱表示相對高度）
  ry <- par("usr")[3] + d$零格比例 * (par("usr")[4] - par("usr")[3])
  rect(amid - 1.6, par("usr")[3], amid + 1.6, ry, col = "#eb68341f", border = NA)
  for (i in seq_along(c(1, 3, 4)))
    lines(amid, d[[ESTS[c(1, 3, 4)][i]]], col = pal[i], lwd = 2, type = "b", pch = 16, cex = .55)
  legend("topright", c(ESTS[c(1, 3, 4)], "零格比例（相對高度）"),
         col = c(pal, "#eb683455"), lwd = c(2, 2, 2, 8), pch = c(16, 16, 16, NA),
         bty = "n", cex = .8)
}
par(op); dev.off()

cat("\n已輸出 tableH / tableI 與 figH_eda_vs_alpha_bias.png\n")
