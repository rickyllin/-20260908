###############################################################################
# 三項已實測有效的方法（Normal 架構）
#
#   1. lc_wsvd()      加權 SVD，權重用配適期望死亡數 mu_hat 而非觀測 D
#                     （修正 Wilmoth 1993 原式的權重內生性）
#   2. lc_heteropca() HeteroPCA（Zhang, Cai & Wu 2022, AoS）
#                     Gram 矩陣對角線被噪音變異數灌水，以秩一近似反覆填補
#   3. bbp_threshold() 尖峰模型相變門檻，預測 beta 何時不可還原
#
# 實測摘要見檔尾。
###############################################################################

## ===================== 共用：給定 beta 解 kappa 並標準化 ===================
.bk_normalize <- function(Z, b) {
  k <- as.vector(crossprod(Z, b)) / sum(b^2)
  k <- k - mean(k)
  s <- sum(b); if (abs(s) < 1e-12) s <- 1
  list(b = b / s, k = k * s)
}

## ===================== 1. 加權 SVD（權重 = mu_hat）=========================
#' @param wmode "mu" 用配適期望死亡數（建議）；"D" 用觀測死亡數（Wilmoth 原式，
#'              有權重內生性，會使 alpha 產生正向偏誤）
#' @param ridge 加在權重上的小常數，使零格得到小但非零的權重
lc_wsvd <- function(D, E, wmode = c("mu", "D"), ridge = 0.375,
                    outer_it = 12, inner_it = 100, tol = 1e-11) {
  wmode <- match.arg(wmode)
  A <- nrow(D); Tn <- ncol(D)
  lm_ <- log(pmax(D, 0.5) / E)
  W <- pmax(D, 0.5) + ridge
  b <- rep(1 / A, A); k <- NULL

  for (o in seq_len(outer_it)) {
    a <- rowSums(W * lm_) / pmax(rowSums(W), 1e-9)
    Z <- lm_ - a
    for (i in seq_len(inner_it)) {                      # 加權交替最小平方
      k <- colSums(W * Z * b) / pmax(colSums(W * b^2), 1e-12)
      k <- k - mean(k)
      bn <- rowSums(W * Z * matrix(k, A, Tn, byrow = TRUE)) /
            pmax(rowSums(W * matrix(k^2, A, Tn, byrow = TRUE)), 1e-12)
      s <- sum(bn); if (abs(s) < 1e-12) break
      if (max(abs(bn / s - b)) < tol) { b <- bn / s; break }
      b <- bn / s
    }
    k <- colSums(W * Z * b) / pmax(colSums(W * b^2), 1e-12); k <- k - mean(k)
    s <- sum(b); b <- b / s; k <- k * s

    if (wmode == "D") break                             # 不迭代權重
    Wn <- E * exp(outer(a, rep(1, Tn)) + outer(b, k)) + ridge
    if (max(abs(Wn - W)) / max(1, max(W)) < 1e-8) { W <- Wn; break }
    W <- Wn
  }
  names(a) <- names(b) <- rownames(D); names(k) <- colnames(D)
  list(a = a, b = b, k = k, W = W)
}

## ===================== 2. HeteroPCA ========================================
#' 由中心化矩陣 Z 取出主方向（對角線以秩一近似反覆填補）
heteropca_vec <- function(Z, r = 1, iters = 30, tol = 1e-10) {
  G <- tcrossprod(Z) / ncol(Z)
  N <- G; diag(N) <- 0
  for (i in seq_len(iters)) {
    ei  <- eigen(N, symmetric = TRUE)
    idx <- seq_len(r)
    Nd  <- ei$vectors[, idx, drop = FALSE] %*%
           diag(pmax(ei$values[idx], 0), r) %*%
           t(ei$vectors[, idx, drop = FALSE])
    Nn <- G; diag(Nn) <- diag(Nd)
    if (max(abs(Nn - N)) < tol) { N <- Nn; break }
    N <- Nn
  }
  u <- eigen(N, symmetric = TRUE)$vectors[, 1]
  if (sum(u) < 0) u <- -u
  u
}

#' HeteroPCA 版的 LC 配適
#' @param weighted TRUE 時先以 mu_hat 加權（兩者可疊加）
lc_heteropca <- function(D, E, weighted = FALSE, zero_sub = 0.5, ...) {
  lm_ <- log(pmax(D, zero_sub) / E)
  if (!weighted) {
    a <- rowMeans(lm_); Z <- lm_ - a
    bk <- .bk_normalize(Z, heteropca_vec(Z, ...))
  } else {
    f <- lc_wsvd(D, E, "mu")
    W <- f$W
    a <- rowSums(W * lm_) / pmax(rowSums(W), 1e-9); Z <- lm_ - a
    s <- sqrt(W / mean(W))
    u <- heteropca_vec(Z * s, ...)
    b <- u / pmax(rowMeans(s), 1e-12)
    bk <- .bk_normalize(Z, b / max(abs(b)))
  }
  names(a) <- names(bk$b) <- rownames(D); names(bk$k) <- colnames(D)
  list(a = a, b = bk$b, k = bk$k)
}

## ===================== 3. BBP 相變門檻 =====================================
#' beta 是否可由 SVD 還原的事前診斷
#'
#' 尖峰模型（Baik–Ben Arous–Péché；矩形版見 Benaych-Georges & Nadakuditi 2012）：
#' 秩一訊號 theta 加上變異數 sigma^2 的白噪音，當 theta/sigma 低於 (A*T)^(1/4)
#' 時，估計出的奇異向量與真值漸近正交——不是估得差，是與真值無關。
#'
#' @param a,b,k 一組 LC 配適（用來構造訊號強度 theta 與噪音尺度）
#' @return list(theta, sigma, threshold, ratio, recoverable)
bbp_threshold <- function(a, b, k, E) {
  A <- length(a); Tn <- length(k)
  mfit <- exp(outer(a, rep(1, Tn)) + outer(b, k))
  Z    <- log(mfit) - rowMeans(log(mfit))
  theta <- svd(Z)$d[1]
  mu    <- E * mfit
  sigma <- sqrt(mean(1 / pmax(mu, 1e-12)))     # log 尺度的等效噪音 sd
  thr   <- (A * Tn)^0.25
  list(theta = theta, sigma = sigma, snr = theta / sigma,
       threshold = thr, ratio = (theta / sigma) / thr,
       recoverable = (theta / sigma) > thr)
}

#' 反解：給定配適期長度，beta 可還原所需的最小人口
bbp_min_population <- function(a, b, k, w_age, lo = 1e3, hi = 1e9, iters = 60) {
  A <- length(a); Tn <- length(k)
  mfit <- exp(outer(a, rep(1, Tn)) + outer(b, k))
  Z <- log(mfit) - rowMeans(log(mfit))
  theta <- svd(Z)$d[1]; thr <- (A * Tn)^0.25
  for (i in seq_len(iters)) {
    mid <- sqrt(lo * hi)
    sg  <- sqrt(mean(1 / pmax(mid * w_age * mfit, 1e-12)))
    if ((theta / sg) / thr > 1) hi <- mid else lo <- mid
  }
  sqrt(lo * hi)
}

###############################################################################
# 實測摘要（女性 2001-2024 為真值；30 次重複）
#
#  N = 2e5
#    方法                 alpha偏誤  SSE(beta)  beta相關  漂移偏誤
#    SVD                   -0.0200     0.0640     0.533    0.1543
#    加權SVD (w = D)        0.0237     0.0369     0.242    0.1159   <- 權重內生性
#    加權SVD (w = mu_hat)  -0.0232     0.0139     0.634    0.0466   <- 修正後
#    HeteroPCA             -0.0200     0.0244     0.620    0.0539
#    HeteroPCA + 加權       0.0206     0.0194     0.496    0.0087
#    卜瓦松 MLE            -0.0010     0.0150     0.669    0.0133
#
#  N = 1e6（大樣本、模型正確——效率代價檢驗）
#    SVD                   -0.0071     0.0122     0.801    0.0492
#    加權SVD (w = mu_hat)  -0.0091     0.0028     0.895    0.0073
#    HeteroPCA             -0.0071     0.0043     0.849    0.0071
#    卜瓦松 MLE            -0.0024     0.0029     0.893   -0.0000
#    -> 權重設對時，Normal 架構不必付出效率代價
#
#  BBP 門檻（theta = 3.444, (A*T)^0.25 = 4.794）
#    N          theta/sigma   比值    預測        實測 beta 相關
#    5e4            4.09      0.854   不可還原       -0.03
#    2e5            8.19      1.708   可還原          0.53
#
#  各配適期的門檻人口 N*
#    T=10: 1,353,066   T=20: 135,279   T=30: 27,020   T=55: 5,401
#    -> 預測了 0929 報告表 12 的每一格；相關係數由約 0.2 跨到約 0.5 的位置
#       在每一列都落在預測門檻上。
#
#  保留：BBP 理論建立在白噪音上，此處以 sqrt(mean(1/mu)) 作等效 sigma 為近似。
#        落點正確，但正式引用前宜補白化後的對照。
###############################################################################
