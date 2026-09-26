###############################################################################
# 變體 F／G：檢查函數能否改善 alpha，以及加權是否為必要條件
#
#   可檢驗的預測（見 研究筆記/0926_穩健估計/轉向穩健估計_素材增減與緒論鋪陳.md）：
#     檢查函數 rho_tau 是等權重的，與 SVD 的等權重最小平方同構，因此
#       F（中位數 LC，等權重）    應改善 alpha，但不應改善漂移項
#       G（加權中位數，w = mu_hat）應同時改善兩者
#   不論成立與否都是可寫的結果：成立則精確劃出穩健損失的作用範圍；
#   不成立則表示機制還有一層尚未辨識。
#
#   種子與 run_eda_to_correction.R 完全相同，故本表可直接併入該表比較。
#   標準 LC 列為內部對照，其數值應與 tableI 逐格相同。
#
# 輸出：output/tables/tableN_quantile_lc.csv
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
source("R/core/quantile_lc.R")

SEED <- 20260926; REPS <- 100
NS <- c(1e4, 5e4, 2e5)           # 與 run_eda_to_correction.R 一致，不可更動
YEARS <- 2001:2024

dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0)

truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn])
drift_true <- (truth$k[Tn] - truth$k[1]) / (Tn - 1)
cat(sprintf("真值：drift = %.4f\n", drift_true))

ESTS <- c("A 標準 LC（對照）", "F 中位數 LC", "G 加權中位數 LC")
fit_one <- function(D, E, est) switch(est,
  "A 標準 LC（對照）" = lc_svd_fit(log(pmax(D, 0.5) / E)),
  "F 中位數 LC"       = lc_quantile(D, E, wmode = "none"),
  "G 加權中位數 LC"   = lc_quantile(D, E, wmode = "mu"))
ok_fit <- function(f) all(is.finite(f$a)) && all(is.finite(f$b)) &&
                      all(is.finite(f$k)) && max(abs(f$k)) < 1e3

res <- list(); res_age <- list()
for (N in NS) {
  E <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu_exp <- E * mtrue
  aM <- array(NA_real_, c(REPS, A, length(ESTS)))
  bS <- dS <- cS <- matrix(NA_real_, REPS, length(ESTS))
  div <- integer(length(ESTS)); zc <- numeric(REPS)

  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000 * which(NS == N) + r)
    D <- matrix(rpois(length(E), mu_exp), A, dimnames = dimnames(E))
    zc[r] <- sum(D == 0)
    for (j in seq_along(ESTS)) {
      f <- tryCatch(fit_one(D, E, ESTS[j]), error = function(e) NULL)
      if (is.null(f) || !ok_fit(f)) { div[j] <- div[j] + 1L; next }
      aM[r, , j] <- f$a - truth$a
      bS[r, j] <- sum((f$b - truth$b)^2) / sum(truth$b^2)
      cS[r, j] <- if (sd(f$b) < 1e-12) NA_real_ else cor(f$b, truth$b)
      dS[r, j] <- (f$k[Tn] - f$k[1]) / (Tn - 1) - drift_true
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }
  res[[as.character(N)]] <- data.frame(
    N = N, 估計量 = ESTS, 平均零格 = round(mean(zc), 1),
    alpha偏誤中位 = round(sapply(seq_along(ESTS), function(j)
      median(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    alpha偏誤最大 = round(sapply(seq_along(ESTS), function(j)
      max(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    SSE_beta = round(apply(bS, 2, median, na.rm = TRUE), 4),
    cor_beta = round(apply(cS, 2, median, na.rm = TRUE), 3),
    漂移偏誤 = round(apply(dS, 2, median, na.rm = TRUE), 4),
    發散率 = round(div / REPS, 3), stringsAsFactors = FALSE)
  res_age[[as.character(N)]] <- data.frame(
    N = N, 年齡 = rownames(D0), 零格比例 = round(rowMeans(mu_exp < 0.7), 2),
    setNames(as.data.frame(lapply(seq_along(ESTS), function(j)
      round(apply(aM[, , j], 2, median, na.rm = TRUE), 4))), ESTS),
    check.names = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("\n===== N = %.0e（平均零格 %.1f/%d）=====\n", N, mean(zc), A * Tn))
  print(res[[as.character(N)]][, -1], row.names = FALSE)
}

tab <- do.call(rbind, res)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(tab, "output/tables/tableN_quantile_lc.csv", row.names = FALSE)
write.csv(do.call(rbind, res_age), "output/tables/tableN_quantile_by_age.csv",
          row.names = FALSE)

cat("\n=== 幼年組（零格集中處）的 alpha 偏誤 ===\n")
for (N in NS) {
  d <- res_age[[as.character(N)]]
  print(d[d$年齡 %in% c("1-4","5-9","10-14","50-54","85-89"), -1], row.names = FALSE)
  cat("\n")
}
cat("已輸出 tableN_quantile_lc.csv\n")
