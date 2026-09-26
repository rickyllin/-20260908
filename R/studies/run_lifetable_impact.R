###############################################################################
# 把估計流程接到生命表：參數的改善能否轉為平均餘命的改善？
#
#   本文第伍節顯示 Firth 使 alpha 的最大偏誤降低 12.7 倍，但該偏誤集中於
#   幼年組，而幼年組的死亡率極低；漂移項的偏誤則直接乘上預測期長度。
#   因此「參數估得準」未必等於「平均餘命估得準」，須直接檢驗。
#
#   第一部分：各估計量的 e0 偏誤（配適末年與預測 10 年後）
#   第二部分：e0 誤差的年齡別歸因——逐一把某年齡組的估計值換回真值，
#             看 e0 誤差減少多少，藉以判斷哪些年齡的偏誤真正影響 e0
#
# 輸出：output/tables/tableO_e0_bias.csv
#       output/tables/tableP_e0_attribution.csv
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/rabbi_mazzuco_replication.R")
source("R/core/heteropca_lc.R")
source("R/core/fh_firth_kalman.R")

SEED <- 20260926; REPS <- 100; H <- 10
NS <- c(1e4, 5e4, 2e5)          # 與 run_eda_to_correction.R 一致
YEARS <- 2001:2024

dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0); ages <- rownames(D0)

truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn])
drift_true <- (truth$k[Tn] - truth$k[1]) / (Tn - 1)

gn <- make_grid(A)$n
e0 <- function(m) life_table(m, gn)$e0
e0_true_fit  <- e0(mtrue[, Tn])
e0_true_fore <- e0(exp(truth$a + truth$b * (truth$k[Tn] + H * drift_true)))
cat(sprintf("真值：配適末年 e0 = %.3f；預測 %d 年後 e0 = %.3f\n",
            e0_true_fit, H, e0_true_fore))

ESTS <- c("標準 LC", "加權 SVD", "卜瓦松 MLE", "Firth", "Firth+FH")
fit_one <- function(D, E, est) switch(est,
  "標準 LC"    = lc_svd_fit(log(pmax(D, 0.5) / E)),
  "加權 SVD"   = lc_wsvd(D, E, wmode = "mu"),
  "卜瓦松 MLE" = lc_poisson_firth(D, E, firth = FALSE),
  "Firth"      = lc_poisson_firth(D, E, firth = TRUE),
  "Firth+FH"   = lc_fh(D, E, fitter = function(D, E) lc_poisson_firth(D, E, firth = TRUE)))
ok <- function(f) all(is.finite(f$a)) && all(is.finite(f$b)) &&
                  all(is.finite(f$k)) && max(abs(f$k)) < 1e3

## ---------------- 第一部分：e0 的偏誤 --------------------------------------
res <- list()
for (N in NS) {
  E <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu_exp <- E * mtrue
  eF <- eP <- matrix(NA_real_, REPS, length(ESTS))
  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000 * which(NS == N) + r)
    D <- matrix(rpois(length(E), mu_exp), A, dimnames = dimnames(E))
    for (j in seq_along(ESTS)) {
      f <- tryCatch(fit_one(D, E, ESTS[j]), error = function(e) NULL)
      if (is.null(f) || !ok(f)) next
      mfit <- exp(f$a + f$b * f$k[Tn])
      dh   <- (f$k[Tn] - f$k[1]) / (Tn - 1)
      mfor <- exp(f$a + f$b * (f$k[Tn] + H * dh))
      eF[r, j] <- e0(mfit) - e0_true_fit
      eP[r, j] <- e0(mfor) - e0_true_fore
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }
  res[[as.character(N)]] <- data.frame(
    N = N, 估計量 = ESTS,
    e0配適偏誤 = round(apply(eF, 2, median, na.rm = TRUE), 3),
    e0配適RMSE = round(sqrt(apply(eF^2, 2, median, na.rm = TRUE)), 3),
    e0預測偏誤 = round(apply(eP, 2, median, na.rm = TRUE), 3),
    e0預測RMSE = round(sqrt(apply(eP^2, 2, median, na.rm = TRUE)), 3),
    stringsAsFactors = FALSE)
  cat(sprintf("\n===== N = %.0e =====\n", N))
  print(res[[as.character(N)]][, -1], row.names = FALSE)
}
tabO <- do.call(rbind, res)

## ---------------- 第二部分：e0 誤差的年齡別歸因 ----------------------------
## 逐一把某年齡組的估計死亡率換回真值，e0 誤差減少的幅度即為該年齡的貢獻。
cat("\n=== e0 誤差的年齡別歸因（N=5e4，標準 LC 與 Firth）===\n")
Nx <- 5e4; RA <- 100
E <- outer(Nx * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0); mu_exp <- E * mtrue
contrib <- array(NA_real_, c(RA, A, 2)); base_err <- matrix(NA_real_, RA, 2)
for (r in seq_len(RA)) {
  set.seed(SEED + 2000 + r)
  D <- matrix(rpois(length(E), mu_exp), A, dimnames = dimnames(E))
  for (j in 1:2) {
    f <- tryCatch(fit_one(D, E, c("標準 LC", "Firth")[j]), error = function(e) NULL)
    if (is.null(f) || !ok(f)) next
    mfit <- exp(f$a + f$b * f$k[Tn]); mt <- mtrue[, Tn]
    base_err[r, j] <- e0(mfit) - e0_true_fit
    for (x in seq_len(A)) {
      mm <- mfit; mm[x] <- mt[x]                  # 只把第 x 組換回真值
      contrib[r, x, j] <- base_err[r, j] - (e0(mm) - e0_true_fit)
    }
  }
}
tabP <- data.frame(
  年齡 = ages,
  真值死亡率 = signif(mtrue[, Tn], 3),
  標準LC貢獻 = round(apply(contrib[, , 1], 2, median, na.rm = TRUE), 4),
  Firth貢獻  = round(apply(contrib[, , 2], 2, median, na.rm = TRUE), 4),
  stringsAsFactors = FALSE)
tabP$標準LC占比 <- round(100 * abs(tabP$標準LC貢獻) / sum(abs(tabP$標準LC貢獻)), 1)
print(tabP, row.names = FALSE)
cat(sprintf("\n基準誤差中位數：標準 LC %.3f 歲、Firth %.3f 歲\n",
            median(base_err[, 1], na.rm = TRUE), median(base_err[, 2], na.rm = TRUE)))
cat(sprintf("30 歲以下各組的貢獻合計：標準 LC %.4f 歲（占 %.1f%%）\n",
            sum(tabP$標準LC貢獻[1:7]),
            100 * abs(sum(tabP$標準LC貢獻[1:7])) / abs(median(base_err[,1], na.rm=TRUE))))

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(tabO, "output/tables/tableO_e0_bias.csv", row.names = FALSE)
write.csv(tabP, "output/tables/tableP_e0_attribution.csv", row.names = FALSE)
cat("\n已輸出 tableO / tableP\n")
