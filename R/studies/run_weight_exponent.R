###############################################################################
# 檢驗 Koenker (2005) 定理 5.1：分位數迴歸的最適權重是局部密度，不是 1/sigma
#
#   背景。本文變體 G 與 K1 以 w = mu_hat 加權，那是 Var(log m_hat) ~ 1/mu
#   的倒數，即「最小平方」的最適權重。但 Koenker (2005, sec. 5.3) 證明
#   加權分位數迴歸的最適權重為該分位處的局部密度 f_i(xi_i)：
#
#     Rather than weighting by the reciprocals of the standard deviations of
#     the observations, quantile regression weights should be proportional to
#     the local density evaluated at the quantile of interest.
#
#   其定理 5.1 給出加權估計量的漸近變異數 tau(1-tau) D2^{-1}(tau)，
#   D2(tau) = lim n^{-1} sum f_i^2(xi_i) x_i x_i'，並證明該估計量相對於
#   未加權版本為無條件更有效率。
#
#   移到本問題。log(D/E) 在 mu 不太小時近似 N(log m, 1/mu)，其中位數處的
#   密度為 sqrt(mu / 2pi)，故檢查函數的最適權重應為 sqrt(mu)，而非 mu。
#
#   可檢驗的預測：以 w = mu^p 掃描冪次 p，最適的 p 應隨 Huber 截點 c
#   由 0.5（c -> 0，檢查函數）移動到 1.0（c -> Inf，最小平方）。
#   若成立，則權重的選擇可由理論而非試誤決定；若不成立（例如最適 p 一律
#   接近 1），則表示本問題的離散性使常態近似失效，該定理不適用——
#   兩種結果都可寫。
#
#   種子與 run_eda_to_correction.R / run_mquantile_scan.R 完全相同。
#
# 輸出：output/tables/tableR_weight_exponent.csv
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
source("R/core/mquantile_lc.R")

SEED <- 20260926; REPS <- 100
NS <- c(1e4, 5e4, 2e5)
YEARS <- 2001:2024

dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0)
truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn])
drift_true <- (truth$k[Tn] - truth$k[1]) / (Tn - 1)
cat(sprintf("真值：drift = %.4f\n", drift_true))

CS   <- c(0, 0.7, 1.345, Inf)           # 穩健程度（0 = 檢查函數，Inf = 最小平方）
POWS <- c(0, 0.25, 0.5, 0.75, 1.0)      # 權重冪次 w = mu^p
GRID <- expand.grid(wpow = POWS, k_c = CS)
GRID$id <- sprintf("c=%-5s p=%.2f",
                   ifelse(is.infinite(GRID$k_c), "Inf", as.character(GRID$k_c)),
                   GRID$wpow)
J <- nrow(GRID)
cat(sprintf("共 %d 個組合 x %d 規模 x %d 次\n", J, length(NS), REPS))

ok_fit <- function(f) all(is.finite(f$a)) && all(is.finite(f$b)) &&
                      all(is.finite(f$k)) && max(abs(f$k)) < 1e3

res <- list()
for (N in NS) {
  E <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu_exp <- E * mtrue
  aM <- array(NA_real_, c(REPS, A, J))
  bS <- dS <- cS <- matrix(NA_real_, REPS, J)
  div <- integer(J)

  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000 * which(NS == N) + r)
    D <- matrix(rpois(length(E), mu_exp), A, dimnames = dimnames(E))
    for (j in seq_len(J)) {
      f <- tryCatch(lc_mquantile(D, E, tau = 0.5, k_c = GRID$k_c[j],
                                 wmode = "pow", wpow = GRID$wpow[j]),
                    error = function(e) NULL)
      if (is.null(f) || !ok_fit(f)) { div[j] <- div[j] + 1L; next }
      aM[r, , j] <- f$a - truth$a
      bS[r, j] <- sum((f$b - truth$b)^2) / sum(truth$b^2)
      cS[r, j] <- if (sd(f$b) < 1e-12) NA_real_ else cor(f$b, truth$b)
      dS[r, j] <- (f$k[Tn] - f$k[1]) / (Tn - 1) - drift_true
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }
  d <- data.frame(
    N = N, c = GRID$k_c, wpow = GRID$wpow, 組合 = GRID$id,
    alpha中位 = round(sapply(seq_len(J), function(j)
      median(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    alpha最大 = round(sapply(seq_len(J), function(j)
      max(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    SSE_beta = round(apply(bS, 2, median, na.rm = TRUE), 4),
    cor_beta = round(apply(cS, 2, median, na.rm = TRUE), 3),
    漂移偏誤 = round(apply(dS, 2, median, na.rm = TRUE), 4),
    發散率 = round(div / REPS, 3), stringsAsFactors = FALSE)
  res[[as.character(N)]] <- d

  cat(sprintf("\n===== N = %.0e =====\n", N))
  for (cc in CS) {
    sub <- d[d$c == cc | (is.infinite(cc) & is.infinite(d$c)), ]
    if (!nrow(sub)) next
    lab <- if (is.infinite(cc)) "Inf" else as.character(cc)
    cat(sprintf("  c = %-5s  ", lab))
    cat(sprintf("p=%.2f:%.4f  ", sub$wpow, sub$alpha中位), sep = "")
    best_a <- sub$wpow[which.min(sub$alpha中位)]
    best_d <- sub$wpow[which.min(abs(sub$漂移偏誤))]
    best_b <- sub$wpow[which.min(sub$SSE_beta)]
    cat(sprintf("| 最適 p：alpha %.2f、漂移 %.2f、beta %.2f\n",
                best_a, best_d, best_b))
  }
}

tab <- do.call(rbind, res)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(tab, "output/tables/tableR_weight_exponent.csv", row.names = FALSE)

cat("\n\n=== 預測的檢驗：最適冪次是否隨 c 由 0.5 移向 1.0 ===\n")
cat(sprintf("%-9s %-8s %-8s %-8s %-8s\n", "N", "c", "最適p(a)", "最適p(漂移)", "最適p(beta)"))
for (N in NS) {
  d <- res[[as.character(N)]]
  for (cc in CS) {
    sub <- d[(d$c == cc) | (is.infinite(cc) & is.infinite(d$c)), ]
    if (!nrow(sub)) next
    cat(sprintf("%-9.0e %-8s %-8.2f %-8.2f %-8.2f\n", N,
                if (is.infinite(cc)) "Inf" else as.character(cc),
                sub$wpow[which.min(sub$alpha中位)],
                sub$wpow[which.min(abs(sub$漂移偏誤))],
                sub$wpow[which.min(sub$SSE_beta)]))
  }
}
cat("\n已輸出 tableR_weight_exponent.csv\n")
