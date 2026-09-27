###############################################################################
# 聯合估計迴歸與尺度：Zoubir et al. (2018, sec. 3.5) 的 M-Lasso 估計方程
#
#   問題。M-分位數的截點 c 必須以殘差的尺度為單位（c_actual = c * sigma），
#   故 sigma 的估計方式直接決定截點落在哪裡。本文原先的作法是每次迭代以
#   MAD 重估 sigma，這是一個外插的權宜作法：MAD 在常態下的效率僅約 37%，
#   且與 c 的取值無關。
#
#   Zoubir et al. (2018, sec. 3.5) 依 Ollila (2016) 的 M-Lasso 估計方程，
#   把迴歸與尺度定義為同一組零次梯度方程的解，尺度方程為
#       (1/N) sum_i psi_c( r_i / sigma )^2 = beta_c,
#       beta_c = E[psi_c(Z)^2], Z ~ N(0,1)
#   即 Huber 的 Proposal 2。此式使 (beta, sigma) 成為同一個目標函數的
#   聯合 M-估計量，且尺度與 c 相容。
#
#   可檢驗的預測：聯合估計應使結果對起始值較不敏感、收斂較快，
#   並在效率上略優於 MAD 版；但因本文的殘差分布並非常態（幼年組為離散、
#   有原子），beta_c 的常態校正未必正確，故改善幅度可能有限甚至為負。
#   兩種結果都可寫：若無改善，即說明「尺度估計不是本問題的瓶頸」。
#
#   種子與 run_mquantile_scan.R 完全相同。
#
# 輸出：output/tables/tableS_joint_scale.csv
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

SPEC <- list(
  list(id = "c=0.7   MAD   等權重", c = 0.7,   sm = "mad",   wm = "none"),
  list(id = "c=0.7   聯合  等權重", c = 0.7,   sm = "joint", wm = "none"),
  list(id = "c=1.345 MAD   等權重", c = 1.345, sm = "mad",   wm = "none"),
  list(id = "c=1.345 聯合  等權重", c = 1.345, sm = "joint", wm = "none"),
  list(id = "c=1.345 MAD   資訊加權", c = 1.345, sm = "mad",   wm = "mu"),
  list(id = "c=1.345 聯合  資訊加權", c = 1.345, sm = "joint", wm = "mu")
)
ESTS <- vapply(SPEC, `[[`, character(1), "id"); J <- length(SPEC)
ok_fit <- function(f) all(is.finite(f$a)) && all(is.finite(f$b)) &&
                      all(is.finite(f$k)) && max(abs(f$k)) < 1e3

res <- list()
for (N in NS) {
  E <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu_exp <- E * mtrue
  aM <- array(NA_real_, c(REPS, A, J))
  bS <- dS <- cS <- itM <- scM <- matrix(NA_real_, REPS, J)
  div <- integer(J)
  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000 * which(NS == N) + r)
    D <- matrix(rpois(length(E), mu_exp), A, dimnames = dimnames(E))
    for (j in seq_len(J)) {
      f <- tryCatch(lc_mquantile(D, E, 0.5, k_c = SPEC[[j]]$c,
                                 wmode = SPEC[[j]]$wm,
                                 scale_mode = SPEC[[j]]$sm),
                    error = function(e) NULL)
      if (is.null(f) || !ok_fit(f)) { div[j] <- div[j] + 1L; next }
      aM[r, , j] <- f$a - truth$a
      bS[r, j] <- sum((f$b - truth$b)^2) / sum(truth$b^2)
      cS[r, j] <- if (sd(f$b) < 1e-12) NA_real_ else cor(f$b, truth$b)
      dS[r, j] <- (f$k[Tn] - f$k[1]) / (Tn - 1) - drift_true
      itM[r, j] <- f$iter; scM[r, j] <- f$scale
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }
  res[[as.character(N)]] <- data.frame(
    N = N, 估計量 = ESTS,
    alpha中位 = round(sapply(seq_len(J), function(j)
      median(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    alpha最大 = round(sapply(seq_len(J), function(j)
      max(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    SSE_beta = round(apply(bS, 2, median, na.rm = TRUE), 4),
    cor_beta = round(apply(cS, 2, median, na.rm = TRUE), 3),
    漂移偏誤 = round(apply(dS, 2, median, na.rm = TRUE), 4),
    尺度中位 = round(apply(scM, 2, median, na.rm = TRUE), 4),
    迭代中位 = round(apply(itM, 2, median, na.rm = TRUE), 1),
    未收斂率 = round(colMeans(itM >= 60, na.rm = TRUE), 2),
    發散率 = round(div / REPS, 3), stringsAsFactors = FALSE)
  cat(sprintf("\n===== N = %.0e =====\n", N))
  print(res[[as.character(N)]][, -1], row.names = FALSE)
}

tab <- do.call(rbind, res)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(tab, "output/tables/tableS_joint_scale.csv", row.names = FALSE)
cat("\n已輸出 tableS_joint_scale.csv\n")
