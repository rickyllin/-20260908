###############################################################################
# 對數轉換偏誤的解析校正：由閉式推導而來的 Lee-Carter 估計量
#
#   推導。標準 LC 對 log m_hat = log( max(D, c0) / E ) 作 SVD，而
#   alpha_x 即該列的平均。故 alpha_x 的偏誤為逐格偏誤的平均，
#   而逐格偏誤有閉式：D ~ Poisson(mu) 時
#
#     b(mu; c0) = E[ log max(D, c0) ] - log mu
#               = e^{-mu} log c0 + sum_{d>=1} P(D=d) log d - log mu.      (*)
#
#   此式的三項性質全部可解析得出：
#     (i)   mu -> 0 時 b -> log(c0) - log(mu) -> +infinity（零格替代高估）
#     (ii)  mu 大時 b -> -1/(2 mu)（Jensen 不等式，實測 mu=300 時
#           b/(-1/(2mu)) = 1.003）
#     (iii) 兩者之間有唯一變號點 mu*，c0 = 0.5 時 mu* = 0.9234；
#           且 mu* 隨 c0 移動（c0 = 0.1/0.3/0.5/1.0 -> mu* =
#           0.134/0.533/0.923/1.509），亦即\textbf{變號點由分析者選定的
#           替代值決定，而非資料給定}。
#
#   實測驗證。以 (*) 預測標準 LC 的逐年齡 alpha 偏誤，與模擬實測值的
#   相關係數在三個人口規模下分別為 1.000、1.000、0.988，迴歸斜率
#   1.001、0.990、0.962，殘差約 0.01（蒙地卡羅誤差量級）。
#   亦即標準 LC 的 alpha 偏誤可由單一個死亡數的函數完全解釋。
#
#   由此得到的估計量。既然偏誤可預測，即可扣除。作法是校正\textbf{格}
#   而非校正參數，使校正同時流入 alpha、beta 與 kappa：
#     1. 以標準 LC 起始，得 mu_hat_{x,t} = E_{x,t} exp(alpha_x + beta_x kappa_t)
#     2. 對每一格計算 b(mu_hat_{x,t})
#     3. 以 log m_hat - b(mu_hat) 重配適
#     4. 反覆至收斂
#
#   方法論歸屬。這屬 Cox and Snell (1968) 的「事後修正」而非 Firth (1993)
#   的「預先修正」，但施加的對象不是最大概似估計量的 O(n^{-1}) 偏誤，
#   而是\textbf{對數轉換與零格替代}所產生的偏誤。此一區別是關鍵：
#   Cox-Snell 需要 MLE 先存在，而本問題的 MLE 可為負無窮；
#   而 (*) 完全不依賴 MLE，只依賴 D 的邊際分布。
#   代價是它依賴卜瓦松假設，且 mu_hat 本身有誤差（故需迭代）。
###############################################################################

#' b(mu; c0) 的直接計算（慢，用於建表與驗證）
logbias_exact <- function(mu, c0 = 0.5) {
  vapply(mu, function(m) {
    if (!is.finite(m) || m <= 0) return(NA_real_)
    K <- max(60L, as.integer(ceiling(m + 12 * sqrt(m) + 12)))
    d <- 0:K; p <- dpois(d, m)
    lg <- ifelse(d == 0, log(c0), log(pmax(d, 1e-300)))
    sum(p * lg) + (1 - sum(p)) * log(max(m, 1)) - log(m)
  }, numeric(1))
}

#' 建立 b(mu) 的查表（對數等距格點 + 線性內插），供模擬中大量呼叫
make_logbias <- function(c0 = 0.5, lo = 1e-4, hi = 1e4, n = 600) {
  g  <- exp(seq(log(lo), log(hi), length.out = n))
  bg <- logbias_exact(g, c0)
  lg <- log(g)
  function(mu) {
    out <- numeric(length(mu))
    m   <- as.vector(mu)
    small <- m <= lo
    big   <- m >= hi
    mid   <- !small & !big & is.finite(m) & m > 0
    ## mu 極小：b = log(c0) - log(mu) 為漸近精確（P(D=0) -> 1）
    out[small] <- log(c0) - log(pmax(m[small], 1e-300))
    ## mu 極大：b = -1/(2mu) 為漸近精確
    out[big]   <- -1 / (2 * m[big])
    out[mid]   <- approx(lg, bg, xout = log(m[mid]), rule = 2)$y
    out[!is.finite(out)] <- 0
    dim(out) <- dim(mu)
    out
  }
}

#' 解析偏誤校正的 Lee-Carter
#'
#' @param shrink 校正的收縮係數。1 為完全扣除；<1 為部分扣除，
#'   用以檢驗「mu_hat 的誤差是否使完全扣除過頭」
lc_analytic <- function(D, E, c0 = 0.5, maxit = 20, tol = 1e-9,
                        shrink = 1, bfun = NULL) {
  if (is.null(bfun)) bfun <- make_logbias(c0)
  A <- nrow(D); Tn <- ncol(D)
  lm_ <- log(pmax(D, c0) / E)
  f <- lc_svd_fit(lm_)
  a <- f$a; b <- f$b; k <- f$k
  for (it in seq_len(maxit)) {
    mu_hat <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    corr   <- shrink * bfun(mu_hat)
    f2 <- lc_svd_fit(lm_ - corr)
    dif <- max(abs(f2$a - a), abs(f2$b - b), abs(f2$k - k) / max(1, max(abs(k))))
    a <- f2$a; b <- f2$b; k <- f2$k
    if (dif < tol) break
  }
  names(a) <- names(b) <- rownames(D); names(k) <- colnames(D)
  list(a = a, b = b, k = k, iter = it, c0 = c0, shrink = shrink)
}

#' 只校正 alpha 的版本（配合資訊加權）
#'
#'   動機。由 Proposition（beta 的分解）可知，逐格扣除 b 會同時做兩件事：
#'     (甲) 移除 alpha 的偏誤，其量為列平均 bbar_x；
#'     (乙) 移除 beta 的確定性擾動 Delta_{x,t} = b(mu_xt) - bbar_x。
#'   但實測顯示 (乙) 只佔 beta 總誤差的一成以下，而逐格扣除所注入的
#'   mu_hat 噪音卻使 beta 明顯變差。既然如此，正確的設計是\textbf{只做 (甲)}：
#'   把 alpha_hat 扣掉 bbar_x，而完全不動中心化矩陣，
#'   如此 beta_hat 與 kappa_hat 不受任何噪音注入，改由加權處理。
#'
#'   由於 alpha_hat 是列平均、在奇異值分解之前即已定出，
#'   事後扣除 bbar_x 不影響 beta_hat 與 kappa_hat 的計算，
#'   兩項改善因而可以疊加而不互相干擾。
#'
#' @param weight 是否以配適期望死亡數加權求解秩一結構
lc_alpha_only <- function(D, E, c0 = 0.5, weight = TRUE,
                          maxit = 20, tol = 1e-9, bfun = NULL) {
  if (is.null(bfun)) bfun <- make_logbias(c0)
  A <- nrow(D); Tn <- ncol(D)
  lm_ <- log(pmax(D, c0) / E)
  f <- lc_svd_fit(lm_); a <- f$a; b <- f$b; k <- f$k
  for (it in seq_len(maxit)) {
    mu_hat <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    bbar   <- rowMeans(bfun(mu_hat))          # 只取列平均，不用逐格值
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
    a2 <- am - bbar                            # 只校正 alpha
    dif <- max(abs(a2 - a), abs(b2 - b), abs(k2 - k) / max(1, max(abs(k))))
    a <- a2; b <- b2; k <- k2
    if (dif < tol) break
  }
  names(a) <- names(b) <- rownames(D); names(k) <- colnames(D)
  list(a = a, b = b, k = k, iter = it, c0 = c0, weight = weight)
}
