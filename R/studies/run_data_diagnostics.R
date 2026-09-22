###############################################################################
# 資料診斷：暴露數對 Lee-Carter 參數估計的影響（逐年齡）
#
# 圖示規格比照 Wang, Yue & Chong (2018) IME 78: 351-359 第 354~355 頁：
#   Fig.2  alpha_x / beta_x 的「偏誤」與「標準誤」，橫軸年齡，按人口規模分層
#   Fig.3  alpha_x / beta_x 的「T-ratio = 偏誤 / 標準誤」，附 +-2 參考線
#
# T-ratio 的意義：偏誤相對於單次估計之標準誤的倍數。
#   |T| > 2 表示該年齡的系統性偏誤已超過抽樣變異的兩倍，
#   亦即「偏誤」而非「變異數」才是該處的主要誤差來源。
#
# 本檔另加一組 SVD vs 卜瓦松 MLE 的對照，用以指認偏誤的來源。
#
# 用法（工作目錄為 repo 根目錄）： Rscript R/studies/run_data_diagnostics.R
###############################################################################

source("R/core/penalized_lc.R")
source("R/studies/mc_exposure.R")
source("R/core/fig_axis_utils.R")

FIGDIR <- "output/figures"; TABDIR <- "output/tables"
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

SEX   <- "Female"; Y0 <- 2001; Y1 <- 2024
REPS  <- 1000
REPS_MLE <- 300                      # MLE 較慢，另設較小的重複次數
SEED  <- 20260922
NS    <- c(1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6, 2e6, 5e6)

truth <- mc_truth(SEX, Y0, Y1)
A <- length(truth$b); ages <- truth$ages; mid <- make_grid(A)$mid
cat(sprintf("真值：%s %d-%d，%d 齡組 x %d 年\n", SEX, Y0, Y1, A, length(truth$k)))
cat(sprintf("重複次數：SVD %d 次、卜瓦松 MLE %d 次；人口規模 %d 個水準\n\n",
            REPS, REPS_MLE, length(NS)))

## ===================== 1. 逐年齡的偏誤與標準誤 ============================
#' @param est  "svd" 或 "mle"
diag_one <- function(N, reps, est = "svd") {
  E <- N * truth$w
  Am <- Bm <- matrix(NA_real_, reps, A)
  for (r in seq_len(reps)) {
    D <- matrix(rpois(length(E), E * truth$m), nrow(E), ncol(E),
                dimnames = dimnames(truth$m))
    f <- if (est == "svd")
      tryCatch(lc_svd_fit(log(pmax(D, 0.5) / E)), error = function(e) NULL)
    else
      tryCatch({ g <- lc_penalized(D, E, lambda = 0, init_beta = "svd")
                 list(a = g$alpha, b = g$beta) }, error = function(e) NULL)
    if (is.null(f)) next
    Am[r, ] <- f$a; Bm[r, ] <- f$b
  }
  ok <- stats::complete.cases(Am) & stats::complete.cases(Bm)
  Am <- Am[ok, , drop = FALSE]; Bm <- Bm[ok, , drop = FALSE]
  data.frame(
    N = N, estimator = est, age = ages, mid = mid, nrep = nrow(Am),
    bias_a = colMeans(Am) - truth$a, se_a = apply(Am, 2, sd),
    bias_b = colMeans(Bm) - truth$b, se_b = apply(Bm, 2, sd),
    true_a = truth$a, true_b = truth$b)
}

set.seed(SEED)
cat("[1/2] SVD 估計\n")
res <- do.call(rbind, lapply(NS, function(N) {
  cat(sprintf("   N = %-9.0f\n", N)); diag_one(N, REPS, "svd") }))
res$t_a <- res$bias_a / res$se_a
res$t_b <- res$bias_b / res$se_b
write.csv(res, file.path(TABDIR, "tableE_diagnostics_svd.csv"), row.names = FALSE)

cat("\n[2/2] 卜瓦松最大概似估計（對照）\n")
res_m <- do.call(rbind, lapply(NS, function(N) {
  cat(sprintf("   N = %-9.0f\n", N)); diag_one(N, REPS_MLE, "mle") }))
res_m$t_a <- res_m$bias_a / res_m$se_a
res_m$t_b <- res_m$bias_b / res_m$se_b
write.csv(res_m, file.path(TABDIR, "tableE_diagnostics_mle.csv"), row.names = FALSE)

## ===================== 2. 繪圖 ============================================
cols <- colorRampPalette(c("grey15", "red", "darkgreen"))(length(NS))
ltys <- rep(c(2, 1, 1), length.out = length(NS))
leg  <- function(pos = "topright", cex = .62)
  legend(pos, legend = format(NS, big.mark = ",", scientific = FALSE),
         col = cols, lty = ltys, lwd = 1.4, bty = "n", cex = cex, ncol = 1)

panel <- function(df, col, ylab, main, href = NULL, ylim = NULL) {
  yy <- df[[col]]
  if (is.null(ylim)) ylim <- range(c(yy, href), na.rm = TRUE, finite = TRUE)
  plot(NA, xlim = range(mid), ylim = ylim, xlab = "Age", ylab = ylab, main = main)
  if (!is.null(href)) abline(h = href, lty = 3, col = "goldenrod3")
  abline(h = 0, lty = 2, col = "blue")
  for (i in seq_along(NS)) {
    s <- df[df$N == NS[i], ]
    lines(s$mid, s[[col]], col = cols[i], lty = ltys[i], lwd = 1.4)
  }
}

## --- 圖 E1：alpha 的偏誤與標準誤 ---
png(file.path(FIGDIR, "figE1_alpha_bias_se.png"), 1200, 520, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4.4, 3, 1))
panel(res, "bias_a", expression(Bias~of~alpha[x]), "Bias of alpha"); leg()
panel(res, "se_a",   expression(s.e.~of~alpha[x]), "s.e. of alpha"); leg()
par(op); dev.off()

## --- 圖 E2：beta 的偏誤與標準誤（對應 WYC Fig.2）---
png(file.path(FIGDIR, "figE2_beta_bias_se.png"), 1200, 520, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4.4, 3, 1))
panel(res, "bias_b", expression(Bias~of~beta[x]), "Bias of beta"); leg()
panel(res, "se_b",   expression(s.e.~of~beta[x]), "s.e. of beta"); leg()
par(op); dev.off()

## --- 圖 E3：T-ratio（對應 WYC Fig.3）---
png(file.path(FIGDIR, "figE3_tratio.png"), 1200, 520, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4.4, 3, 1))
ylim_t <- range(c(res$t_a, res$t_b, -2.5, 2.5), na.rm = TRUE, finite = TRUE)
ylim_t <- c(max(ylim_t[1], -8), min(ylim_t[2], 8))
panel(res, "t_a", expression(Bias/s.e.~of~alpha[x]), "T-ratio of alpha",
      href = c(-2, 2), ylim = ylim_t); leg("bottomright")
panel(res, "t_b", expression(Bias/s.e.~of~beta[x]),  "T-ratio of beta",
      href = c(-2, 2), ylim = ylim_t); leg("bottomright")
par(op); dev.off()

## --- 圖 E4：SVD vs 卜瓦松 MLE 的 T-ratio 對照 ---
png(file.path(FIGDIR, "figE4_tratio_svd_vs_mle.png"), 1200, 520, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4.4, 3, 1))
for (nm in c("t_a", "t_b")) {
  lab <- if (nm == "t_a") "alpha" else "beta"
  yy <- c(res$t_a, res$t_b, res_m$t_a, res_m$t_b)
  yl <- c(max(min(yy, na.rm=TRUE), -8), min(max(yy, na.rm=TRUE), 8))
  plot(NA, xlim = range(mid), ylim = yl, xlab = "Age",
       ylab = sprintf("Bias/s.e. of %s", lab),
       main = sprintf("T-ratio of %s: SVD vs Poisson MLE", lab))
  abline(h = c(-2, 2), lty = 3, col = "goldenrod3"); abline(h = 0, lty = 2, col = "blue")
  for (i in seq_along(NS)) {
    s1 <- res  [res  $N == NS[i], ]; s2 <- res_m[res_m$N == NS[i], ]
    lines(s1$mid, s1[[nm]], col = "grey35", lwd = 1.2)
    lines(s2$mid, s2[[nm]], col = "red",    lwd = 1.2)
  }
  legend("bottomright", c("SVD", "Poisson MLE"), col = c("grey35","red"),
         lwd = 1.6, bty = "n", cex = .75)
}
par(op); dev.off()

## ===================== 3. 摘要表 ==========================================
smry <- do.call(rbind, lapply(NS, function(N) {
  s  <- res  [res  $N == N, ]; m <- res_m[res_m$N == N, ]
  data.frame(N = N,
             sig_a_svd = sum(abs(s$t_a) > 2), sig_b_svd = sum(abs(s$t_b) > 2),
             sig_a_mle = sum(abs(m$t_a) > 2), sig_b_mle = sum(abs(m$t_b) > 2),
             max_abs_t_a_svd = max(abs(s$t_a)), max_abs_t_b_svd = max(abs(s$t_b)),
             max_abs_t_a_mle = max(abs(m$t_a)), max_abs_t_b_mle = max(abs(m$t_b)))
}))
cat(sprintf("\n=== |T-ratio| > 2 的年齡組數（共 %d 組）===\n", A))
print(smry[, c("N","sig_a_svd","sig_b_svd","sig_a_mle","sig_b_mle")],
      row.names = FALSE)
write.csv(smry, file.path(TABDIR, "tableE_tratio_summary.csv"), row.names = FALSE)
cat(sprintf("\n圖已輸出至 %s/（figE1~figE4），表已輸出至 %s/\n", FIGDIR, TABDIR))
