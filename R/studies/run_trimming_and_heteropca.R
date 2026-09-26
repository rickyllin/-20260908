###############################################################################
# 穩健估計的兩側：截尾類（負面）與 HeteroPCA（正面）
#
#   截尾類（MTL、LTS）以「對抗少數污染觀測」為設計目標。本文變體 E 已證明
#   無離群值亦有系統性偏誤，故預期無效——本腳本驗證之，並指出失效的機制：
#   兩種截尾規則丟掉的都是最具資訊的格，只是丟法不同。
#
#   HeteroPCA（Zhang, Cai & Wu 2022）則直接針對異質變異：樣本 Gram 矩陣的
#   對角線被噪音變異數灌水，該法以秩一近似反覆填補對角線。
#
# 輸出：output/tables/tableJ_robust_ladder.csv    估計量階梯（含截尾與 HeteroPCA）
#       output/tables/tableK_trim_retention.csv   逐年齡保留率（失效機制的證據）
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/rabbi_mazzuco_replication.R")
source("R/core/heteropca_lc.R")
source("R/core/fh_firth_kalman.R")
source("R/core/reml_kappa.R")

SEED <- 20260926; REPS <- 100; NS <- c(5e4, 2e5)
YEARS <- 2001:2024

dat  <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0); ages <- rownames(D0)

truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn])
drift_true <- (truth$k[Tn] - truth$k[1]) / (Tn - 1)

firth_fit <- function(D, E) lc_poisson_firth(D, E, firth = TRUE)

ESTS <- c("標準 LC", "卜瓦松 MLE", "Firth", "HeteroPCA", "HeteroPCA+FH",
          "MTL（保留 90%）", "LTS（保留 90%）")

fit_one <- function(D, E, est) {
  switch(est,
    "標準 LC"         = lc_svd_fit(log(pmax(D, 0.5) / E)),
    "卜瓦松 MLE"      = lc_poisson_firth(D, E, firth = FALSE),
    "Firth"           = firth_fit(D, E),
    "HeteroPCA"       = lc_heteropca(D, E),
    "HeteroPCA+FH"    = { f <- lc_heteropca(D, E)
                          v <- 1 / info_beta(f$a, f$b, f$k, E)
                          s <- fh_shrink(f$b, v)
                          list(a = f$a, b = s$beta, k = f$k) },
    "MTL（保留 90%）" = trimmed_lc(D, E, frac = 0.90),
    "LTS（保留 90%）" = lts_lc(D, E, frac = 0.90))
}
ok_fit <- function(f) all(is.finite(f$a)) && all(is.finite(f$b)) &&
                      all(is.finite(f$k)) && max(abs(f$k)) < 1e3

## ---------------- 第一部分：估計量階梯 -------------------------------------
res <- list()
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
      bS[r, j]   <- sum((f$b - truth$b)^2) / sum(truth$b^2)
      cS[r, j]   <- if (sd(f$b) < 1e-12) NA_real_ else cor(f$b, truth$b)
      dS[r, j]   <- (f$k[Tn] - f$k[1]) / (Tn - 1) - drift_true
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
  cat(sprintf("\n===== N = %.0e（平均零格 %.1f/%d）=====\n", N, mean(zc), A * Tn))
  print(res[[as.character(N)]][, -1], row.names = FALSE)
}
tabJ <- do.call(rbind, res)

## ---------------- 第二部分：截尾掉的是哪些格 -------------------------------
## 於 N = 5e6（模型正確、無零格、無離群值）檢查逐年齡保留率。
## 若截尾是中性的，各年齡的保留率應大致相同。
cat("\n=== 逐年齡保留率（N=5e6，模型正確、無零格）===\n")
NBIG <- 5e6; RB <- 20
Eb <- outer(NBIG * wage, rep(1, Tn)); dimnames(Eb) <- dimnames(D0)
keepM <- keepL <- matrix(NA_real_, RB, A)
for (r in seq_len(RB)) {
  set.seed(SEED + 9000 + r)
  D <- matrix(rpois(length(Eb), Eb * mtrue), A, dimnames = dimnames(Eb))
  fm <- tryCatch(trimmed_lc(D, Eb, 0.90), error = function(e) NULL)
  fl <- tryCatch(lts_lc(D, Eb, 0.90),     error = function(e) NULL)
  if (!is.null(fm)) keepM[r, ] <- rowMeans(fm$W)
  if (!is.null(fl)) keepL[r, ] <- rowMeans(fl$W)
}
tabK <- data.frame(
  年齡 = ages,
  beta_true = round(truth$b, 4),
  MTL保留率 = round(colMeans(keepM, na.rm = TRUE), 2),
  LTS保留率 = round(colMeans(keepL, na.rm = TRUE), 2),
  期望死亡中位 = round(apply(Eb * mtrue, 1, median), 0),
  stringsAsFactors = FALSE)
print(tabK, row.names = FALSE)

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(tabJ, "output/tables/tableJ_robust_ladder.csv", row.names = FALSE)
write.csv(tabK, "output/tables/tableK_trim_retention.csv", row.names = FALSE)
cat("\n已輸出 tableJ / tableK\n")
