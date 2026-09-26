###############################################################################
# 優先順序第 4、5、7 項的實作 + BBP 對齊公式（已依原文核對）
#
#   7. fh_shrink()     Fay-Herriot 經驗貝氏收縮 —— 免調參取代 lambda 選擇
#   5. lc_poisson_firth()  Firth (1993) 偏誤減少 —— 消除稀疏資料下的發散
#   4. kalman_jumpoff()    Kalman 濾波的 kappa_{T|T} 與其變異數 P_T
#      bbp_alignment()     BBP 相變：門檻與門檻以上的還原程度
#
#   第 6 項（Gavish-Donoho 最適奇異值收縮）經實測無效，故不提供；
#   原因：它作用在奇異「值」上，而問題在奇異「向量」，beta 估計值完全不變。
###############################################################################

## ===================== 7. Fay-Herriot 經驗貝氏 =============================
#' beta_x 的 Fisher 資訊（即報告式(10)的 W_x，無額外計算成本）
info_beta <- function(a, b, k, E) {
  Tn <- length(k)
  mu <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
  pmax(as.vector(mu %*% k^2), 1e-12)
}

#' Fay-Herriot 收縮
#'   bhat_x = b_x + e_x,  Var(e_x) = v_x（由資訊矩陣得到，已知）
#'   b_x ~ N(b0_x, Atau)，Atau 以 (RE)ML 估計
#'   B_x = Atau/(Atau+v_x) 逐年齡不同：資訊薄的年齡收縮重、資訊厚的幾乎不收縮
#'
#' @param b0 收縮目標，預設 1/A（「各年齡以相同速度改善」的虛無假設）
#' @return list(beta, B, Atau) —— beta 已重新標準化使 sum = 1
fh_shrink <- function(bhat, v, b0 = NULL, restricted = TRUE) {
  A <- length(bhat)
  if (is.null(b0)) b0 <- rep(1 / A, A)
  d <- bhat - b0
  negll <- function(lA) {
    s <- exp(lA) + v
    val <- sum(log(s) + d^2 / s)
    if (restricted) val <- val + log(sum(1 / s))   # 對 b0 的水準積分
    0.5 * val
  }
  opt  <- optimize(negll, c(log(1e-12), log(1e2)))
  Atau <- exp(opt$minimum)
  B    <- Atau / (Atau + v)
  bt   <- b0 + B * d
  s    <- sum(bt)
  list(beta = if (abs(s) > 1e-12) bt / s else bt, B = B, Atau = Atau)
}

#' 一步到位：配適卜瓦松 LC 後做 FH 收縮
#' 注意：不要迭代（收縮 -> 重配適 kappa -> 再收縮）。實測顯示迭代版更差，
#'       且部分重複中 beta 會完全塌陷為常數。
lc_fh <- function(D, E, fitter = lc_poisson_safe, b0 = NULL) {
  f <- fitter(D, E)
  a <- if (!is.null(f$a)) f$a else f$alpha
  b <- if (!is.null(f$b)) f$b else f$beta
  k <- if (!is.null(f$k)) f$k else f$kappa
  v  <- 1 / info_beta(a, b, k, E)
  sh <- fh_shrink(b, v, b0)
  list(a = a, b = sh$beta, b_unshrunk = b, k = k,
       B = sh$B, Atau = sh$Atau)
}

## ===================== 5. Firth 偏誤減少 ===================================
#' 雙線性 LC 的帽子矩陣對角線
#'   Jacobian： d eta_xt/d alpha_a = 1{x=a}
#'              d eta_xt/d beta_a  = kappa_t * 1{x=a}
#'              d eta_xt/d kappa_s = beta_x * 1{t=s}
#'   I = J' W J（W = diag(mu)）；識別條件使 I 秩虧 2，故用廣義逆。
lc_leverage <- function(a, b, k, E) {
  A <- nrow(E); Tn <- ncol(E); n <- A * Tn; p <- 2 * A + Tn
  mu <- as.vector(t(E * exp(outer(a, rep(1, Tn)) + outer(b, k))))  # 列優先
  xi <- rep(seq_len(A), each = Tn); ti <- rep(seq_len(Tn), times = A)
  J <- matrix(0, n, p)
  J[cbind(seq_len(n), xi)]             <- 1
  J[cbind(seq_len(n), A + xi)]         <- k[ti]
  J[cbind(seq_len(n), 2 * A + ti)]     <- b[xi]
  Jw <- J * sqrt(mu)
  Ii <- MASS::ginv(crossprod(Jw))
  h  <- rowSums((Jw %*% Ii) * Jw)
  matrix(h, A, Tn, byrow = TRUE)
}

#' Firth 修正的卜瓦松 LC
#'   把得分方程中的 (D - mu) 換成 (D - mu + h/2)，
#'   等價於對概似乘上 Jeffreys 事前 |I(theta)|^{1/2}。
#'   主要效益：完全消除稀疏資料下的發散（N=1e4 時 35% -> 0%）。
lc_poisson_firth <- function(D, E, maxit = 400, tol = 1e-10, firth = TRUE) {
  A <- nrow(D); Tn <- ncol(D)
  a <- log(pmax(rowSums(D), 0.5) / rowSums(E))
  b <- rep(1 / A, A); k <- seq(1, -1, length.out = Tn)
  for (it in seq_len(maxit)) {
    h  <- if (firth) lc_leverage(a, b, k, E) / 2 else matrix(0, A, Tn)
    mu <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    a  <- a + rowSums(D - mu + h) / pmax(rowSums(mu), 1e-12)
    mu <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    k  <- k + colSums((D - mu + h) * b) / pmax(colSums(mu * b^2), 1e-12)
    k  <- k - mean(k)
    mu <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    bn <- b + rowSums((D - mu + h) * matrix(k, A, Tn, byrow = TRUE)) /
              pmax(rowSums(mu * matrix(k^2, A, Tn, byrow = TRUE)), 1e-12)
    s <- sum(bn); if (abs(s) < 1e-12) break
    nb <- bn / s; nk <- k * s
    del <- max(max(abs(nb - b)), max(abs(nk - k)) / max(1, max(abs(k))))
    b <- nb; k <- nk
    if (del < tol) break
  }
  names(a) <- names(b) <- rownames(D); names(k) <- colnames(D)
  list(a = a, b = b, k = k, iter = it)
}

## ===================== 4. Kalman 初始值與其變異數 ==========================
#' 狀態空間：kappa_hat_t = kappa_t + u_t,  u_t ~ N(0, tau_t^2)
#'           kappa_t = kappa_{t-1} + delta + e_t,  e_t ~ N(0, sigma^2)
#' 回傳最後一點的濾波估計（在 t = T 即等於平滑估計）與其變異數。
#'
#' 實測重點：點估計 kappa_{T|T} 只帶來邊際改善，真正有用的是 P_T
#'           —— 把初始值本身的不確定性納入預測分布後涵蓋率明顯提升。
kalman_jumpoff <- function(khat, tau2, sigma, delta) {
  m <- khat[1]; P <- tau2[1]
  for (t in seq_along(khat)[-1]) {
    mp <- m + delta; Pp <- P + sigma^2
    K  <- Pp / (Pp + tau2[t])
    m  <- mp + K * (khat[t] - mp)
    P  <- (1 - K) * Pp
  }
  list(kappa_T = m, var_T = P)
}

#' e0 的預測分位數，可選是否納入初始值的不確定性
forecast_e0_q <- function(a, b, k_last, drift, sigma, h, life_fun,
                          jumpoff_var = 0, M = 2000, probs = c(.1, .9)) {
  paths <- k_last + t(apply(matrix(rnorm(M * h, drift, sigma), M, h), 1, cumsum))
  if (jumpoff_var > 0) paths <- paths + rnorm(M, 0, sqrt(jumpoff_var))
  e0 <- apply(paths, 2, function(col)
    vapply(col, function(x) life_fun(exp(a + b * x)), numeric(1)))
  apply(e0, 2, quantile, probs)
}

## ===================== BBP 相變（已依原文核對）=============================
#' Benaych-Georges & Nadakuditi (2012), JMVA 111: 120-135, §3.1 式 (9)(10)
#'
#' 設定：n x m 矩陣（n <= m）、噪音各元變異數 1/m、c = n/m。
#'   門檻   theta > c^(1/4)
#'   對齊   |<u_hat, u>|^2 -> 1 - c(1+theta^2)/(theta^2(theta^2+c))，門檻以下為 0
#'
#' 對應本設定：n = A, m = T, 各元噪音變異數 sigma^2 時 theta = s1/(sigma*sqrt(T))。
#'
#' 驗證：theta >= 1.2*c^(1/4) 時公式與模擬誤差 < 0.05；門檻附近與以下有明顯的
#'       有限樣本偏離——理論極限 0，但 A = 22 時隨機單位向量的期望內積約 0.17，
#'       故實務上「門檻以下」應讀成「與隨機方向無法區分」。
bbp_alignment <- function(a, b, k, E) {
  A <- length(a); Tn <- length(k); c_ <- A / Tn
  mfit <- exp(outer(a, rep(1, Tn)) + outer(b, k))
  Z  <- log(mfit) - rowMeans(log(mfit))
  s1 <- svd(Z)$d[1]
  sg <- sqrt(mean(1 / pmax(E * mfit, 1e-12)))     # 等效白噪音 sd（近似）
  th <- s1 / (sg * sqrt(Tn))
  thr <- c_^0.25
  al2 <- if (th >= thr) max(1 - c_ * (1 + th^2) / (th^2 * (th^2 + c_)), 0) else 0
  list(theta = th, threshold = thr, ratio = th / thr,
       alignment = sqrt(al2), recoverable = th > thr,
       noise_floor = sqrt(2 / (pi * A)))          # 隨機方向的期望內積
}

#' 反解：給定配適期長度，beta 可還原所需的最小人口
bbp_min_population <- function(a, b, k, w_age, lo = 1e3, hi = 1e9, iters = 60) {
  A <- length(a); Tn <- length(k); c_ <- A / Tn
  mfit <- exp(outer(a, rep(1, Tn)) + outer(b, k))
  s1 <- svd(log(mfit) - rowMeans(log(mfit)))$d[1]
  for (i in seq_len(iters)) {
    mid <- sqrt(lo * hi)
    sg  <- sqrt(mean(1 / pmax(mid * w_age * mfit, 1e-12)))
    if (s1 / (sg * sqrt(Tn)) > c_^0.25) hi <- mid else lo <- mid
  }
  sqrt(lo * hi)
}

###############################################################################
# 實測摘要（女性 2001-2024 為真值）
#
#  7. Fay-Herriot（免調參）vs 以真值挑出的最佳格點 lambda
#     N = 5e4   未懲罰 SSE 0.07695 / 相關 0.310
#               FH-EB  SSE 0.01024 / 相關 0.418   <- 相關最高
#               LASSO 0.20（格點最佳） SSE 0.01005 / 相關 0.307
#     N = 2e5   FH-EB  SSE 0.00651 / 相關 0.660
#               LASSO 0.03（格點最佳） SSE 0.00622 / 相關 0.680
#     -> 兩個規模都與 oracle 調參的懲罰同水準，而 FH 完全由資料決定。
#        與可實現的方法（交叉驗證，報告表 9 相對損失 3.49）相比優勢更大。
#     逐年齡收縮係數 B_x（N=2e5）：5-9 歲 0.13、80-84 歲 0.94、100+ 0.39
#        —— 資訊薄的重收縮、資訊厚的幾乎不收縮，是導出來的而非調出來的。
#
#  5. Firth
#     N        卜瓦松發散   Firth發散   卜瓦松SSE   FirthSSE
#     1e4      14/40(35%)   0/40        1.6600     0.8943
#     2e4       5/40(13%)   0/40        0.6824     0.3384
#     5e4       0/40        0/40        0.0849     0.0708
#     N = 5e4 時 alpha 偏誤 -0.0064 -> -0.0012，但漂移偏誤 -0.0040 -> -0.0342。
#     等效 pseudo-count h/2 中位數 0.0492（範圍 0.020-0.146），比慣例的 0.5
#     小一個數量級。註：卜瓦松的 SSE 只計收斂樣本，是倖存者樣本，比較對其有利。
#
#  4. Kalman（h=10、名目 80%）
#     N = 5e4   不可約下限 1.727/0.87；用 kappa_hat_T 1.981/0.72；
#               用 kappa_{T|T} 1.974/0.74；再加計 P_T 2.110/0.77
#     N = 2e5   不可約下限 1.715/0.76；用 kappa_hat_T 1.429/0.68；
#               用 kappa_{T|T} 1.360/0.69；再加計 P_T 1.517/0.74
#     -> 點估計只有邊際效益，P_T 才是關鍵（寬度代價僅 6-7%）。
#
#  6. Gavish-Donoho 最適收縮（不提供，實測無效）
#     SSE(beta) 與相關係數一位數字都沒變（0.40456 -> 0.40456），
#     因為收縮作用在 s1 上而 beta = u1/sum(u1) 只取決於 u1。
#     漂移則明顯變差（N=2e5 偏誤 0.0928 -> 0.3139）。
#     作為診斷工具有效（其門檻與 BBP 一致），作為估計方法無效。
###############################################################################
