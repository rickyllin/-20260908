###############################################################################
# 狀態空間 REML 估計 kappa 的創新變異數，以及截尾類估計量（對照用）
#
#   問題：kappa_hat_t = kappa_t + u_t，故 Var(diff(kappa_hat)) = sigma^2 + Var(diff(u))
#         動差法（報告 §六、八 步驟二）以拔靴估 Var(diff(u*)) 後相減，
#         但 N = 5 萬時兩個大噪音量相減不穩定，僅 10% 可行。
#
#   本檔改以狀態空間表述：
#         觀測  kappa_hat_t = kappa_t + u_t,        u_t ~ N(0, tau_t^2)
#         狀態  kappa_t = kappa_{t-1} + delta + e_t, e_t ~ N(0, sigma^2)
#   其中 tau_t^2 由卜瓦松資訊矩陣直接代入（免估計），
#   sigma^2 在概似中被估計而非相減得到，結構上不可能為負。
#
#   需要：lc_poisson_safe()（或任何回傳 a, b, k 的卜瓦松 LC 配適）
###############################################################################

## ===================== 1. 由卜瓦松配適取出 tau_t^2 =========================
#' kappa_t 的估計變異數（Fisher 資訊的倒數）
#'
#' 給定 alpha, beta，kappa_t 的 Fisher 資訊為 sum_x mu_xt * beta_x^2。
#' 這一項就是卜瓦松 LC 迭代中 kappa 更新式的分母，故無額外計算成本。
tau2_from_info <- function(a, b, k, E) {
  Tn <- ncol(E)
  mu <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
  1 / pmax(as.vector(crossprod(mu, b^2)), 1e-12)
}

## ===================== 2. 狀態空間 ML / REML ===============================
#' 估計 kappa 的創新標準差 sigma
#'
#' 差分後 d_t = diff(kappa_hat)_t = delta + e_t + u_t - u_{t-1}，其共變異為
#'   Var(d_t)      = sigma^2 + tau_t^2 + tau_{t-1}^2
#'   Cov(d_t,d_t+1) = -tau_t^2
#' delta 為固定效果。restricted = TRUE 時對 delta 積分（REML），
#' 即在目標函數中加入 log|X' Sigma^{-1} X| 一項。
#'
#' @param k    估計出的 kappa 序列（長度 T）
#' @param tau2 各年的估計變異數（長度 T），由 tau2_from_info() 取得
#' @param restricted TRUE 為 REML，FALSE 為 ML
#' @return list(sigma, delta, loglik, converged)
sigma_state_space <- function(k, tau2, restricted = TRUE,
                              interval = c(1e-8, 1e4)) {
  d  <- diff(k); n <- length(d)
  t2 <- as.numeric(tau2)
  stopifnot(length(t2) == n + 1)

  negll <- function(log_s2) {
    s2 <- exp(log_s2)
    ## 三對角共變異矩陣
    S <- diag(s2 + t2[-1] + t2[-(n + 1)], n)
    if (n > 1) {
      idx <- cbind(1:(n - 1), 2:n)
      S[idx] <- -t2[2:n]
      S[idx[, c(2, 1)]] <- -t2[2:n]
    }
    ch <- tryCatch(chol(S), error = function(e) NULL)
    if (is.null(ch)) return(1e10)
    Si  <- chol2inv(ch)
    one <- rep(1, n)
    q   <- sum(one %*% Si %*% one)        # X' Sigma^{-1} X
    dh  <- as.numeric((one %*% Si %*% d) / q)
    r   <- d - dh
    val <- 2 * sum(log(diag(ch))) + as.numeric(r %*% Si %*% r)
    if (restricted) val <- val + log(q)
    0.5 * val
  }

  opt <- optimize(negll, interval = log(interval))
  s2  <- exp(opt$minimum)

  ## 取回 delta
  S <- diag(s2 + t2[-1] + t2[-(n + 1)], n)
  if (n > 1) {
    idx <- cbind(1:(n - 1), 2:n)
    S[idx] <- -t2[2:n]; S[idx[, c(2, 1)]] <- -t2[2:n]
  }
  Si <- chol2inv(chol(S)); one <- rep(1, n)
  delta <- as.numeric((one %*% Si %*% d) / sum(one %*% Si %*% one))

  list(sigma = sqrt(s2), delta = delta, loglik = -opt$objective,
       converged = opt$minimum > log(interval[1]) * 0.999 &&
                   opt$minimum < log(interval[2]) * 0.999)
}

#' 一次完成：配適卜瓦松 LC → 取 tau^2 → 估 sigma
#' @param fitter 回傳 list(a, b, k) 的卜瓦松 LC 配適函式
lc_sigma <- function(D, E, fitter = lc_poisson_safe, restricted = TRUE) {
  f  <- fitter(D, E)
  kk <- if (!is.null(f$k)) f$k else f$kappa
  aa <- if (!is.null(f$a)) f$a else f$alpha
  bb <- if (!is.null(f$b)) f$b else f$beta
  t2 <- tau2_from_info(aa, bb, kk, E)
  ss <- sigma_state_space(kk, t2, restricted)
  Tn <- length(kk)
  list(a = aa, b = bb, k = kk,
       sigma = ss$sigma,
       sigma_naive = sd(diff(kk) - (kk[Tn] - kk[1]) / (Tn - 1)),
       drift_ls = (kk[Tn] - kk[1]) / (Tn - 1),   # 端點漂移（報告所用）
       drift_gls = ss$delta,                      # 加權漂移（狀態空間）
       tau = sqrt(t2), converged = ss$converged)
}

## ===================== 3. 截尾類估計量（對照組）===========================
#' 加權卜瓦松 LC；W 為 0/1 權重矩陣
lc_poisson_w <- function(D, E, W, maxit = 3000, tol = 1e-11) {
  A <- nrow(D); Tn <- ncol(D)
  a <- log(pmax(rowSums(W * D), 0.5) / pmax(rowSums(W * E), 1e-9))
  b <- rep(1 / A, A); k <- seq(1, -1, length.out = Tn)
  for (it in seq_len(maxit)) {
    mu <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    a  <- a + log(pmax(rowSums(W * D), 0.5) / pmax(rowSums(W * mu), 1e-12))
    mu <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    k  <- k + colSums(W * (D - mu) * b) / pmax(colSums(W * mu * b^2), 1e-12)
    k  <- k - mean(k)
    mu <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    bn <- b + rowSums(W * (D - mu) * matrix(k, A, Tn, byrow = TRUE)) /
              pmax(rowSums(W * mu * matrix(k^2, A, Tn, byrow = TRUE)), 1e-12)
    s <- sum(bn); if (abs(s) < 1e-10) break
    nb <- bn / s; nk <- k * s
    del <- max(max(abs(nb - b)), max(abs(nk - k)) / max(1, max(abs(k))))
    b <- nb; k <- nk
    if (del < tol) break
  }
  list(a = a, b = b, k = k)
}

#' 最大截尾概似（FAST-TLE 式迭代重加權）
#' @param frac 保留比例
trimmed_lc <- function(D, E, frac = 0.90, iters = 6) {
  W <- matrix(1, nrow(D), ncol(D)); h <- round(frac * length(D))
  for (i in seq_len(iters)) {
    f  <- lc_poisson_w(D, E, W)
    mu <- E * exp(outer(f$a, rep(1, ncol(D))) + outer(f$b, f$k))
    ll <- D * log(pmax(mu, 1e-300)) - mu - lgamma(D + 1)
    thr <- sort(ll, decreasing = TRUE)[h]
    Wn  <- (ll >= thr) * 1
    if (identical(Wn, W)) break
    W <- Wn
  }
  c(lc_poisson_w(D, E, W), list(W = W))
}

#' log 尺度的最小截尾平方
lts_lc <- function(D, E, frac = 0.90, iters = 6, zero_sub = 0.5) {
  lm_ <- log(pmax(D, zero_sub) / E); h <- round(frac * length(D))
  W <- matrix(1, nrow(D), ncol(D))
  for (i in seq_len(iters)) {
    am <- rowSums(W * lm_) / pmax(rowSums(W), 1e-9)
    Z  <- (lm_ - am) * W
    sv <- svd(Z); u1 <- sv$u[, 1]; v1 <- sv$v[, 1]
    if (sum(u1) < 0) { u1 <- -u1; v1 <- -v1 }
    b <- u1 / sum(u1); k <- sv$d[1] * v1 * sum(u1)
    r2 <- (lm_ - am - outer(b, k))^2
    Wn <- (r2 <= sort(r2)[h]) * 1
    if (identical(Wn, W)) break
    W <- Wn
  }
  list(a = am, b = b, k = k, W = W)
}

###############################################################################
# 實測摘要（女性 2001-2024 為真值，sigma 真值 0.7496）
#
#  sigma 的估計：
#    N         方法                中位數    偏誤      可行率
#    5e4       動差法去偏（報告）    2.4360   +1.6865     10%
#              狀態空間 ML          0.6899   -0.0596    100%
#              狀態空間 REML        0.8974   +0.1479    100%
#    2e5       動差法去偏           0.6705   -0.0790     85%
#              狀態空間 REML        0.6058   -0.1438    100%
#    1e6       動差法去偏           0.7081   -0.0415    100%
#              狀態空間 REML        0.6623   -0.0872    100%
#
#  -> N = 5 萬（動差法失效處）狀態空間法偏誤小一個數量級且 100% 可行，
#     把去偏的操作下限由 N >= 20 萬推進到至少 N = 5 萬。
#     N >= 20 萬時動差法的中位數偏誤較小、狀態空間法的 IQR 較小，兩法互補。
#     ML 與 REML 互有勝負（小 N 時 ML 較佳、大 N 時 REML 較佳），宜並列報告。
#
#  預測區間（h = 10、名目 80%）：
#    N = 5e4：不可約下限寬度 1.717；現行作法 7.126（涵蓋 0.40）；
#             動差法 4.081（涵蓋 0.99，靠 2.4 倍寬度換來）；
#             狀態空間 REML 1.640（涵蓋 0.67）
#    -> 寬度已貼合下限，殘餘缺口在位置（kappa_hat_T 這個初始值），非寬度。
#
#  截尾類（對照，一律不建議採用）：
#    N = 5e6（模型正確、無零格）時 LTS 的 SSE(beta) 為卜瓦松 MLE 的 43 倍、
#    beta 相關由 0.981 崩至 0.258。逐年齡檢查顯示 LTS 丟掉 63% 的 1-4 歲組
#    與 58% 的 5-9 歲組——正是 beta_x 最大的兩個年齡；MTL 則偏好丟掉
#    70-94 歲，即死亡數最多的年齡。兩種截尾規則丟掉的都是最具資訊的格。
###############################################################################
