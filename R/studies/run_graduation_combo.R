###############################################################################
# 修勻與中心化校正的組合（分支研究：修勻 × 貝氏混合，第一問）
#
#   問題。表 tab:smr 顯示兩件事。其一，Whittaker 比值修勻把 alpha 偏誤的
#   「中位數」由 0.0568 壓到 0.0095（優於本文中心化加權的 0.0223），
#   但「最大值」只由 0.8658 降到 0.7622、SSE(beta) 只由 6.71 降到 6.41、
#   e0 偏誤完全不動（-1.601 對 -1.672）。其二，本文的中心化加權把最大值
#   壓到 0.0721、e0 偏誤壓到 +0.003，但 SSE(beta) 停在 1.12、
#   e0 標準差 0.564 仍不如 Firth 的 0.422。
#
#   亦即兩者處理的是不相交的病灶：修勻處理非零格的跨年齡震盪，
#   中心化處理零格替代與 Jensen 偏誤。若此解讀正確，組合起來應
#   同時改善兩側，且改善量可事先預測。本程式即檢驗此一預測。
#
#   三個插入位置。
#     (甲) 修勻（僅）：對 log 死亡率逐年修勻後配適，不使用參考母體。
#          這是把 Yue et al. (2019) 的「先修勻後配模型」改成自含式版本，
#          以區隔修勻本身的貢獻與參考母體的貢獻。
#     (乙) 修勻 -> 中心化＋加權：先修勻，再對修勻後的資料施以中心化。
#     (丙) 中心化＋加權（mu_hat 先修勻）：資料不動，只把代入 b(.) 的
#          mu_hat 先修勻。這一項的理論動機最明確——b(mu_hat) 是插入式
#          估計，mu_hat 的噪音逐格注入是 beta 變差的已知機制，
#          而 b 在小 mu 處為凸，故 E[b(mu_hat)] - b(mu) 另有一項
#          約 b''(mu) Var(mu_hat) / 2。修勻 mu_hat 同時壓低兩者。
#
#   修勻的三項實作選擇，皆須明寫。
#     1. 權重取 mu_hat（由前置的標準 Lee-Carter 配適取得）而非人口數 E。
#        log(D/E) 的變異約為 1/mu，故資訊權重是 mu 而非 E。
#     2. 差分懲罰只施於 5-9 歲以上的二十個年齡組。年齡組為 0、1-4、
#        5-9、…、100+，前兩組的組距與其餘不同，且嬰幼兒死亡率本身
#        在年齡上陡降、並非平滑。對前兩組施以二階差分懲罰會扭曲
#        正是校正最需要處理的兩個年齡。
#     3. h = lambda * mean(w)，掃描 lambda 以呈現對平滑量的敏感度。
#        lambda 的選取是分析者的選擇，與 c0 同性質；原則上的解法是
#        把 h 寫成變異成分再以 REML 定出（見分支計畫第二問）。
#
# 輸出：output/tables/tableG1_graduation_combo.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/partial_smr.R")

SEED <- 20261006; REPS <- 100; NS <- c(1e4, 5e4); YEARS <- 2001:2024
LAM  <- c(0.3, 1, 3, 10); PEN_FROM <- 3L     # 由 5-9 歲起施以懲罰

dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0)
truth  <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue  <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage   <- E0[, Tn] / sum(E0[, Tn])
e0_true <- life_table(mtrue[, Tn])$e0
BF <- make_logbias(0.5)

#' 自含式 Whittaker 修勻：對 A x T 的對數尺度矩陣逐年修勻
grad_mat <- function(L, W, lam, z = 2, from = PEN_FROM) {
  idx <- from:nrow(L); out <- L
  for (j in seq_len(ncol(L))) {
    w <- W[idx, j]
    out[idx, j] <- whittaker_1d(L[idx, j], w = w, h = lam * mean(w), z = z)
  }
  out
}

#' lc_alpha_only 的推廣版：資料端可換成修勻後的對數死亡率，
#' 且可選擇在代入 b(.) 之前先把 mu_hat 修勻。
lc_alpha_gen <- function(lm_, E, bfun, weight = TRUE, lam_mu = NULL,
                         maxit = 20, tol = 1e-9) {
  Tn <- ncol(lm_)
  f <- lc_svd_fit(lm_); a <- f$a; b <- f$b; k <- f$k
  for (it in seq_len(maxit)) {
    mu_hat <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    mu_b   <- if (is.null(lam_mu)) mu_hat else
                E * exp(grad_mat(log(mu_hat / E), mu_hat, lam_mu))
    bbar   <- rowMeans(bfun(mu_b))
    if (weight) {
      W  <- mu_hat / mean(mu_hat)
      am <- rowSums(W * lm_) / rowSums(W)
      Z  <- (lm_ - am) * sqrt(W)
      sv <- svd(Z); u1 <- sv$u[, 1]; v1 <- sv$v[, 1]
      if (sum(u1) < 0) { u1 <- -u1; v1 <- -v1 }
      b2 <- u1 / sum(u1)
      k2 <- sv$d[1] * v1 * sum(u1) / sqrt(pmax(colMeans(W), 1e-12))
      k2 <- k2 - mean(k2); sb <- sum(b2); b2 <- b2 / sb; k2 <- k2 * sb
    } else {
      f2 <- lc_svd_fit(lm_); am <- f2$a; b2 <- f2$b; k2 <- f2$k
    }
    a2 <- am - bbar
    dif <- max(abs(a2 - a), abs(b2 - b), abs(k2 - k) / max(1, max(abs(k))))
    a <- a2; b <- b2; k <- k2
    if (dif < tol) break
  }
  names(a) <- names(b) <- rownames(lm_); names(k) <- colnames(lm_)
  list(a = a, b = b, k = k, iter = it)
}

lab <- c("標準 Lee-Carter", "Firth", "中心化＋加權",
         sprintf("甲 修勻（僅）λ=%g", LAM),
         sprintf("乙 修勻→中心化＋加權 λ=%g", LAM),
         sprintf("丙 中心化＋加權（μ̂ 修勻）λ=%g", LAM))
K <- length(lab)

res <- list()
for (ni in seq_along(NS)) {
  N <- NS[ni]
  E  <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu <- E * mtrue
  aM <- array(NA_real_, c(REPS, A, K)); bS <- eS <- matrix(NA_real_, REPS, K)

  for (r in seq_len(REPS)) {
    set.seed(SEED + 13000 * ni + r)
    D   <- matrix(rpois(length(E), mu), A, dimnames = dimnames(E))
    lm_ <- log(pmax(D, 0.5) / E)
    f0  <- lc_svd_fit(lm_)                       # 前置配適，供修勻的權重用
    Wmu <- E * exp(outer(f0$a, rep(1, Tn)) + outer(f0$b, f0$k))

    fits <- vector("list", K)
    fits[[1]] <- f0
    fits[[2]] <- tryCatch(lc_poisson_firth(D, E, firth = TRUE), error = function(e) NULL)
    fits[[3]] <- tryCatch(lc_alpha_only(D, E, bfun = BF, weight = TRUE),
                          error = function(e) NULL)
    j <- 3L
    for (lam in LAM) { j <- j + 1L
      fits[[j]] <- tryCatch(lc_svd_fit(grad_mat(lm_, Wmu, lam)), error = function(e) NULL) }
    for (lam in LAM) { j <- j + 1L
      fits[[j]] <- tryCatch(lc_alpha_gen(grad_mat(lm_, Wmu, lam), E, BF, TRUE),
                            error = function(e) NULL) }
    for (lam in LAM) { j <- j + 1L
      fits[[j]] <- tryCatch(lc_alpha_gen(lm_, E, BF, TRUE, lam_mu = lam),
                            error = function(e) NULL) }

    for (j in seq_len(K)) {
      f <- fits[[j]]
      if (is.null(f) || !all(is.finite(f$a)) || !all(is.finite(f$b))) next
      aM[r, , j] <- f$a - truth$a
      bS[r, j]   <- sum((f$b - truth$b)^2) / sum(truth$b^2)
      mh <- exp(outer(f$a, rep(1, Tn)) + outer(f$b, f$k))
      eS[r, j] <- tryCatch(life_table(mh[, Tn])$e0, error = function(e) NA_real_)
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }

  res[[ni]] <- data.frame(
    N = N, 估計量 = lab,
    alpha中位 = round(sapply(seq_len(K), function(j)
      median(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    alpha最大 = round(sapply(seq_len(K), function(j)
      max(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    SSE_beta  = round(apply(bS, 2, median, na.rm = TRUE), 4),
    e0偏誤    = round(apply(eS, 2, median, na.rm = TRUE) - e0_true, 3),
    e0標準差  = round(apply(eS, 2, sd, na.rm = TRUE), 3),
    e0均方根  = round(sqrt(apply(eS, 2, function(v)
      mean((v - e0_true)^2, na.rm = TRUE))), 3),
    stringsAsFactors = FALSE)
  cat(sprintf("\n===== N = %.0e （真值 e0 = %.3f）=====\n", N, e0_true))
  print(res[[ni]][, -1], row.names = FALSE)
}
tab <- do.call(rbind, res)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(tab, "output/tables/tableG1_graduation_combo.csv", row.names = FALSE)
cat("\n已輸出 tableG1_graduation_combo.csv\n")
