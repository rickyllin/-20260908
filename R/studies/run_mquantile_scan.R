###############################################################################
# 變體 H-K：把「穩健程度」由開關改成連續刻度，並檢驗兩項可事先寫下的預測
#
#   背景。第伍節之六以變體 A（標準 LC）與變體 F（中位數 LC）兩點劃出穩健
#   損失的作用範圍，結論是中位數的崩潰點 50% 即其上界。但這兩點其實是
#   同一個損失函數族的兩個極限（見 R/core/mquantile_lc.R 的說明）：
#       M-分位數 rho_{tau,c}，c -> 0 為分位數迴歸、c -> Inf 為非對稱最小平方
#   因此可以把 (tau, c) 當成連續刻度掃描，檢驗：
#
#   預測一（穩健程度）。零格比例低於崩潰點時（N = 2e5，幼年組 30-36%），
#     最適的 c 應偏小（偏向中位數）；零格比例遠高於崩潰點時（N = 1e4，
#     幼年組 93-95%），中位數本身已被扭曲，最適的 c 應偏大（偏向最小平方）。
#     亦即最適穩健程度應隨人口規模單調移動。
#
#   預測二（非對稱方向）。零格以 0.5 替代必然「高估」該格死亡率，故被扭曲
#     的是殘差的右尾。若如此，壓抑右尾（tau < 0.5，或 c_pos < c_neg）應優於
#     壓抑左尾。注意第柒節之三（二）原先的推測方向相反，本研究以此檢驗之。
#
#   預測三（Hyndman-Ullah）。Hyndman and Ullah (2007) 的穩健化以「年」為
#     單位。小人口的噪音在每一格獨立、每一年同樣嘈雜，故不存在離群年，
#     逐年穩健化應毫無作用。
#
#   種子與 run_eda_to_correction.R / run_quantile_lc.R 完全相同，
#   故變體 A 一列可與 tableI、tableN 逐格對照，作為內部對照。
#
# 輸出：output/tables/tableQ_mquantile_scan.csv     （總表）
#       output/tables/tableQ_mquantile_by_age.csv   （逐年齡 alpha 偏誤）
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
source("R/core/quantile_lc.R")
source("R/core/mquantile_lc.R")
source("R/core/robust_fpca_lc.R")

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

## --- 估計量清單 -------------------------------------------------------------
## 組別 1：穩健程度的掃描（tau = 0.5，等權重）——檢驗預測一
## 組別 2：非對稱參數的掃描（c = 1.345，等權重）——檢驗預測二
## 組別 3：非對稱截點（Xu and Chen 2018）——檢驗預測二的另一種實作
## 組別 4：與資訊加權疊加
## 組別 5：Hyndman-Ullah 逐年穩健化——檢驗預測三
SPEC <- list(
  list(id = "A  標準 LC（對照）",            g = 1, f = function(D, E) lc_svd_fit(log(pmax(D, 0.5) / E))),
  list(id = "H1 M-q c=3.0",                  g = 1, f = function(D, E) lc_mquantile(D, E, 0.5, k_c = 3.0)),
  list(id = "H2 M-q c=1.345（Huber 標準）",  g = 1, f = function(D, E) lc_mquantile(D, E, 0.5, k_c = 1.345)),
  list(id = "H3 M-q c=0.7",                  g = 1, f = function(D, E) lc_mquantile(D, E, 0.5, k_c = 0.7)),
  list(id = "H4 M-q c=0（分位數極限）",      g = 1, f = function(D, E) lc_mquantile(D, E, 0.5, k_c = 0)),
  list(id = "I1 M-q tau=0.3",                g = 2, f = function(D, E) lc_mquantile(D, E, 0.3, k_c = 1.345)),
  list(id = "I2 M-q tau=0.4",                g = 2, f = function(D, E) lc_mquantile(D, E, 0.4, k_c = 1.345)),
  list(id = "I3 M-q tau=0.6",                g = 2, f = function(D, E) lc_mquantile(D, E, 0.6, k_c = 1.345)),
  list(id = "I4 M-q tau=0.7",                g = 2, f = function(D, E) lc_mquantile(D, E, 0.7, k_c = 1.345)),
  list(id = "J1 只壓右尾 c+=1.345",          g = 3, f = function(D, E) lc_mquantile(D, E, 0.5, k_cn = Inf,   k_cp = 1.345)),
  list(id = "J2 只壓左尾 c-=1.345",          g = 3, f = function(D, E) lc_mquantile(D, E, 0.5, k_cn = 1.345, k_cp = Inf)),
  list(id = "K1 M-q c=1.345 + 資訊加權",     g = 4, f = function(D, E) lc_mquantile(D, E, 0.5, k_c = 1.345, wmode = "mu")),
  list(id = "K2 M-q c=0 + 資訊加權（≈G）",   g = 4, f = function(D, E) lc_mquantile(D, E, 0.5, k_c = 0,     wmode = "mu")),
  list(id = "L1 HU 逐年降權 lambda=3",       g = 5, f = function(D, E) lc_robust_year(D, E, rule = "huber", lambda = 3)),
  list(id = "L2 HU 逐年剔除 2 年",           g = 5, f = function(D, E) lc_robust_year(D, E, rule = "trim", n_trim = 2)),
  list(id = "L3 HU 逐年剔除 4 年",           g = 5, f = function(D, E) lc_robust_year(D, E, rule = "trim", n_trim = 4))
)
ESTS <- vapply(SPEC, `[[`, character(1), "id")
GRP  <- vapply(SPEC, `[[`, numeric(1),  "g")
J <- length(SPEC)

ok_fit <- function(f) all(is.finite(f$a)) && all(is.finite(f$b)) &&
                      all(is.finite(f$k)) && max(abs(f$k)) < 1e3

res <- list(); res_age <- list()
for (N in NS) {
  E <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu_exp <- E * mtrue
  aM <- array(NA_real_, c(REPS, A, J))
  bS <- dS <- cS <- matrix(NA_real_, REPS, J)
  div <- integer(J); zc <- numeric(REPS); ndown <- numeric(REPS)

  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000 * which(NS == N) + r)
    D <- matrix(rpois(length(E), mu_exp), A, dimnames = dimnames(E))
    zc[r] <- sum(D == 0)
    for (j in seq_len(J)) {
      f <- tryCatch(SPEC[[j]]$f(D, E), error = function(e) NULL)
      if (is.null(f) || !ok_fit(f)) { div[j] <- div[j] + 1L; next }
      if (j == which(ESTS == "L1 HU 逐年降權 lambda=3")) ndown[r] <- f$n_down
      aM[r, , j] <- f$a - truth$a
      bS[r, j] <- sum((f$b - truth$b)^2) / sum(truth$b^2)
      cS[r, j] <- if (sd(f$b) < 1e-12) NA_real_ else cor(f$b, truth$b)
      dS[r, j] <- (f$k[Tn] - f$k[1]) / (Tn - 1) - drift_true
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }
  res[[as.character(N)]] <- data.frame(
    N = N, 組 = GRP, 估計量 = ESTS, 平均零格 = round(mean(zc), 1),
    alpha偏誤中位 = round(sapply(seq_len(J), function(j)
      median(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    alpha偏誤最大 = round(sapply(seq_len(J), function(j)
      max(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    SSE_beta = round(apply(bS, 2, median, na.rm = TRUE), 4),
    cor_beta = round(apply(cS, 2, median, na.rm = TRUE), 3),
    漂移偏誤 = round(apply(dS, 2, median, na.rm = TRUE), 4),
    發散率 = round(div / REPS, 3), stringsAsFactors = FALSE)
  res_age[[as.character(N)]] <- data.frame(
    N = N, 年齡 = rownames(D0), 零格比例 = round(rowMeans(mu_exp < 0.7), 2),
    期望死亡 = round(rowMeans(mu_exp), 1),
    setNames(as.data.frame(lapply(seq_len(J), function(j)
      round(apply(aM[, , j], 2, median, na.rm = TRUE), 4))), ESTS),
    check.names = FALSE, stringsAsFactors = FALSE)
  cat(sprintf("\n===== N = %.0e（平均零格 %.1f/%d；HU 平均降權年數 %.2f）=====\n",
              N, mean(zc), A * Tn, mean(ndown)))
  print(res[[as.character(N)]][, -1], row.names = FALSE)
}

tab <- do.call(rbind, res)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(tab, "output/tables/tableQ_mquantile_scan.csv", row.names = FALSE)
write.csv(do.call(rbind, res_age), "output/tables/tableQ_mquantile_by_age.csv",
          row.names = FALSE)

## --- 預測二的逐年齡檢驗：最適 tau 是否隨零格比例變號 --------------------
cat("\n\n=== 逐年齡 alpha 偏誤：tau 的方向是否隨零格比例改變 ===\n")
TCOL <- c("I1 M-q tau=0.3", "I2 M-q tau=0.4", "H2 M-q c=1.345（Huber 標準）",
          "I3 M-q tau=0.6", "I4 M-q tau=0.7")
for (N in NS) {
  d <- res_age[[as.character(N)]]
  sub <- d[, c("年齡", "零格比例", "期望死亡", TCOL)]
  names(sub)[3 + seq_along(TCOL)] <- c("t.30", "t.40", "t.50", "t.60", "t.70")
  sub$最適tau <- c(0.3, 0.4, 0.5, 0.6, 0.7)[apply(abs(as.matrix(
    sub[, c("t.30", "t.40", "t.50", "t.60", "t.70")])), 1, which.min)]
  cat(sprintf("\n--- N = %.0e ---\n", N))
  print(sub, row.names = FALSE)
}
cat("\n已輸出 tableQ_mquantile_scan.csv / tableQ_mquantile_by_age.csv\n")
