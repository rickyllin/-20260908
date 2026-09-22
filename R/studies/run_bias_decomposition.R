###############################################################################
# 偏誤來源的階梯分解，與估計後的校正
#
# 動機（2026-09-22 會議）：老師指出「原始資料本身會有對數的偏誤，
# 但是否還有其他因素影響」。本檔以階梯式設計逐一移除候選因素，
# 檢驗對數轉換是否為唯一來源。
#
# 候選因素：
#   (1) 零死亡格的替代值（D=0 以 0.5 取代）
#   (2) 對數的二階偏誤（Jensen）：E[log D] ~= log(mu) - 1/(2 mu)
#   (3) SVD 的等權重最小平方（Wilmoth 1993 指出應以死亡數加權）
#   (4) 奇異向量的訊噪比擾動（與對數無關，任何含噪矩陣皆有）
#
# 估計量階梯：
#   A 標準 LC        SVD on log(max(D,0.5)/E)        (1)(2)(3)(4) 全在
#   B ＋Jensen 校正   SVD on log(D/E) + 1/(2D)        移除 (2)
#   C 加權 SVD        以 D 為權重的加權最小平方        移除 (3)
#   D 卜瓦松 MLE      不取對數                        移除 (1)(2)(3)
#   E 高斯對照        對真值 log m 加上均值為零、
#                     變異數與卜瓦松相當的噪音後做 SVD  只剩 (4)
#
#   -> 若 E 仍有偏誤，即證明對數並非唯一來源。
#
# 第二部分：估計後的拔靴法偏誤校正
#   theta_corrected = 2*theta_hat - mean(theta_bootstrap)
#
# 用法： Rscript R/studies/run_bias_decomposition.R
###############################################################################

source("R/core/penalized_lc.R")
source("R/studies/mc_exposure.R")

TABDIR <- "output/tables"; FIGDIR <- "output/figures"
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

REPS <- 400; SEED <- 20260922
NS   <- c(2e4, 5e4, 2e5, 1e6)
truth <- mc_truth("Female", 2001, 2024)
A <- length(truth$b); Tn <- length(truth$k)
drift_true <- (truth$k[Tn] - truth$k[1]) / (Tn - 1)
logm_true <- log(truth$m)

## 以死亡數為權重的加權 SVD（Wilmoth 1993）：以交替最小平方求秩一近似
wsvd_fit <- function(Z, W, iter = 50) {
  a <- rowSums(W * Z) / pmax(rowSums(W), 1e-12)
  R <- Z - a
  b <- rep(1 / nrow(Z), nrow(Z)); k <- as.vector(crossprod(b, R))
  for (i in seq_len(iter)) {
    b <- rowSums(W * R * matrix(k, nrow(R), ncol(R), byrow = TRUE)) /
         pmax(rowSums(W * matrix(k^2, nrow(R), ncol(R), byrow = TRUE)), 1e-12)
    k <- colSums(W * R * b) / pmax(colSums(W * b^2), 1e-12)
    k <- k - mean(k)
    s <- sum(b); if (abs(s) > 1e-12) { b <- b / s; k <- k * s }
  }
  list(a = a, b = b, k = k)
}

fit_one <- function(D, E, est) {
  switch(est,
    A = lc_svd_fit(log(pmax(D, 0.5) / E)),
    B = { Z <- log(pmax(D, 0.5) / E) + 1 / (2 * pmax(D, 0.5))   # Jensen 一階校正
          lc_svd_fit(Z) },
    C = wsvd_fit(log(pmax(D, 0.5) / E), W = pmax(D, 0.5)),
    D = { g <- lc_penalized(D, E, lambda = 0, init_beta = "svd")
          list(a = g$alpha, b = g$beta, k = g$kappa) })
}

## ===================== 1. 階梯分解 ========================================
set.seed(SEED)
rows <- list()
for (N in NS) {
  E <- N * truth$w
  sd_gauss <- sqrt(1 / pmax(E * truth$m, 1e-8))   # var(log(D/E)) ~= 1/mu
  acc <- setNames(lapply(c("A","B","C","D","E"), function(z)
           list(a = matrix(NA_real_, REPS, A), b = matrix(NA_real_, REPS, A),
                d = numeric(REPS))), c("A","B","C","D","E"))
  zero <- numeric(REPS)
  for (r in seq_len(REPS)) {
    D <- matrix(rpois(length(E), E * truth$m), nrow(E), ncol(E),
                dimnames = dimnames(truth$m))
    zero[r] <- sum(D == 0)
    for (est in c("A","B","C","D")) {
      f <- tryCatch(fit_one(D, E, est), error = function(e) NULL)
      if (is.null(f)) next
      acc[[est]]$a[r, ] <- f$a; acc[[est]]$b[r, ] <- f$b
      acc[[est]]$d[r]   <- (f$k[Tn] - f$k[1]) / (Tn - 1)
    }
    # E：高斯對照——無對數轉換、無零格、變異數與卜瓦松相當
    Zg <- logm_true + matrix(rnorm(length(E), 0, sd_gauss), nrow(E), ncol(E))
    fg <- lc_svd_fit(Zg)
    acc$E$a[r, ] <- fg$a; acc$E$b[r, ] <- fg$b
    acc$E$d[r]   <- (fg$k[Tn] - fg$k[1]) / (Tn - 1)
  }
  # 少數模擬會發散（SSE 可達中位數的 10^6 倍），故一律採用穩健統計：
  # 逐年齡取中位數估偏誤、逐次模擬的 SSE 取中位數。
  for (est in c("A","B","C","D","E")) {
    S <- acc[[est]]; ok <- stats::complete.cases(S$a) & stats::complete.cases(S$b)
    Am <- S$a[ok,,drop=FALSE]; Bm <- S$b[ok,,drop=FALSE]
    med_a <- apply(Am, 2, median); med_b <- apply(Bm, 2, median)
    sse_r <- rowSums(sweep(Bm, 2, truth$b)^2)
    rows[[length(rows)+1]] <- data.frame(
      N = N, zero_cells = mean(zero), estimator = est, nrep = sum(ok),
      bias_a_med  = median(med_a - truth$a),          # 逐年齡中位數偏誤的中位數
      max_abs_bias_a = max(abs(med_a - truth$a)),
      sse_b_med   = median(sse_r),                    # 逐次 SSE 的中位數
      sse_b_bias  = sum((med_b - truth$b)^2),         # 中位數估計的偏誤平方和
      diverge_rate = mean(sse_r > 100 * median(sse_r)),
      cor_b       = median(apply(Bm, 1, function(v)
                      if (sd(v) < 1e-12) NA else cor(v-mean(v), truth$b-mean(truth$b))),
                      na.rm = TRUE),
      drift_bias  = median(S$d[ok]) - drift_true)
  }
  cat(sprintf("N = %-9.0f 完成（平均零格 %.1f）\n", N, mean(zero)))
}
dec <- do.call(rbind, rows)
lab <- c(A="A 標準 LC", B="B +Jensen 校正", C="C 加權 SVD",
         D="D 卜瓦松 MLE", E="E 高斯對照（只剩訊噪比）")
dec$label <- lab[dec$estimator]
write.csv(dec, file.path(TABDIR, "tableF_bias_decomposition.csv"), row.names = FALSE)

cat("\n=== 偏誤來源的階梯分解 ===\n")
for (N in NS) {
  cat(sprintf("\n-- N = %.0e（平均零格 %.1f / %d）--\n",
              N, dec$zero_cells[dec$N==N][1], A*Tn))
  s <- dec[dec$N == N, c("label","bias_a_med","sse_b_med","cor_b","drift_bias","diverge_rate")]
  print(s, digits = 3, row.names = FALSE)
}

## ===================== 2. 估計後的拔靴法偏誤校正 ==========================
cat("\n\n=== 估計後校正：參數式拔靴法的偏誤校正 ===\n")
B_BOOT <- 60; REPS2 <- 150
set.seed(SEED + 1)
corr <- do.call(rbind, lapply(c(5e4, 2e5), function(N) {
  E <- N * truth$w
  raw <- cor_ <- matrix(NA_real_, REPS2, A); draw <- dcor <- numeric(REPS2)
  for (r in seq_len(REPS2)) {
    D <- matrix(rpois(length(E), E*truth$m), nrow(E), ncol(E), dimnames=dimnames(truth$m))
    f <- lc_svd_fit(log(pmax(D,0.5)/E))
    mhat <- exp(outer(f$a, rep(1,Tn)) + outer(f$b, f$k))
    Bst <- matrix(NA_real_, B_BOOT, A); dst <- numeric(B_BOOT)
    for (bb in seq_len(B_BOOT)) {
      Db <- matrix(rpois(length(E), E*mhat), nrow(E), ncol(E))
      fb <- lc_svd_fit(log(pmax(Db,0.5)/E))
      Bst[bb, ] <- fb$b; dst[bb] <- (fb$k[Tn]-fb$k[1])/(Tn-1)
    }
    raw[r, ]  <- f$b
    cor_[r, ] <- 2*f$b - colMeans(Bst)                  # 偏誤校正
    draw[r]   <- (f$k[Tn]-f$k[1])/(Tn-1)
    dcor[r]   <- 2*draw[r] - mean(dst)
  }
  data.frame(N = N,
    bias2_raw = sum((apply(raw,2,median)-truth$b)^2),
    bias2_cor = sum((apply(cor_,2,median)-truth$b)^2),
    sse_med_raw = median(rowSums(sweep(raw,2,truth$b)^2)),
    sse_med_cor = median(rowSums(sweep(cor_,2,truth$b)^2)),
    drift_bias_raw = median(draw)-drift_true, drift_bias_cor = median(dcor)-drift_true)
}))
print(corr, digits = 3, row.names = FALSE)
write.csv(corr, file.path(TABDIR, "tableG_bootstrap_bias_correction.csv"), row.names = FALSE)
cat(sprintf("\n已輸出 %s/tableF_*.csv 與 tableG_*.csv\n", TABDIR))
