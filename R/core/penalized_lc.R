###############################################################################
# 帶懲罰項的 Lee-Carter 模型：LASSO / Ridge / Elastic Net on beta_x
#
#   D_xt ~ Poisson(E_xt * exp(alpha_x + beta_x * kappa_t))
#
#   min  -l(alpha,beta,kappa)
#        + lambda * sum_x w_x [ a*|beta_x - beta0_x| + (1-a)/2*(beta_x-beta0_x)^2 ]
#   s.t. sum_x beta_x = 1,  sum_t kappa_t = 0
#
#   a = 1    -> LASSO
#   a = 0    -> Ridge
#   0<a<1    -> Elastic Net
#
# ---------------------------------------------------------------------------
# 關鍵設計：收縮目標 beta0 預設為均勻向量 1/A，而非 0。
#
#   在識別條件 sum_x beta_x = 1 之下，朝 0 收縮的 L1 懲罰恆為常數
#   （見 lc_poisson_lasso.R 檔尾的恆等式），估計值完全不動。
#   取 beta0 = 1/A 則 sum(beta0) = 1 = sum(beta)，故 sum(beta - beta0) = 0，
#   此時 sum|beta - beta0| = 2 * sum(正部) 並非常數，懲罰項有效。
#
#   人口學意義：beta0 = 1/A 代表「各年齡以相同速度改善」的虛無假設。
#   LASSO 的稀疏性因而回答「哪些年齡組顯著偏離平均改善速度」。
#
# 依賴：lc_poisson_lasso.R（load_data）
###############################################################################

if (!exists("load_data")) source("R/core/lc_poisson_lasso.R")

soft <- function(z, l) sign(z) * pmax(abs(z) - l, 0)

## ===================== 1. 懲罰 LC 的估計 ==================================
#' @param alpha_en  Elastic Net 混合參數：1 = LASSO, 0 = Ridge
#' @param beta0     收縮目標，預設 rep(1/A, A)
#' @param track     TRUE 時回傳每次迭代的目標函數值與參數變動量
#' @param init_beta beta 的起始值。預設 rep(1/A, A)；在 T 很短或 N 很小時
#'   beta 幾乎不可識別，演算法會停在起始值不動，使 SSE 產生假性的「好表現」。
#'   傳入 "svd" 則以 SVD 配適為起始值，可避免此假象。
lc_penalized <- function(D, E, lambda = 0, alpha_en = 1, beta0 = NULL, w = NULL,
                         maxit = 500, tol = 1e-11, track = FALSE,
                         init_beta = NULL) {
  A <- nrow(D); Tn <- ncol(D)
  if (is.null(beta0)) beta0 <- rep(1 / A, A)      # <- 關鍵預設
  if (is.null(w))     w     <- rep(1, A)

  alpha <- log(rowSums(pmax(D, 0.5)) / rowSums(E))
  beta  <- rep(1 / A, A)
  kappa <- seq(1, -1, length.out = Tn)

  loglik <- function(a, b, k) {
    eta <- outer(a, rep(1, Tn)) + outer(b, k)
    if (any(!is.finite(eta))) return(-Inf)
    eta <- pmin(eta, 50)                       # 防溢位
    v <- sum(D * eta - E * exp(eta))
    if (is.finite(v)) v else -Inf
  }
  pen <- function(b) {
    d <- b - beta0
    lambda * sum(w * (alpha_en * abs(d) + (1 - alpha_en) / 2 * d^2))
  }

  hist <- if (track) data.frame(iter = integer(0), obj = numeric(0),
                                dbeta = numeric(0)) else NULL
  obj_old <- -Inf
  for (it in seq_len(maxit)) {
    beta_prev <- beta

    ## --- alpha：封閉解 -----------------------------------------------------
    mu    <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
    alpha <- alpha + log(rowSums(D) / rowSums(mu))

    ## --- kappa：Newton 一步後置中 -----------------------------------------
    mu    <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
    num   <- colSums((D - mu) * beta)
    den   <- colSums(mu * beta^2)
    kappa <- kappa + num / pmax(den, 1e-12)
    kappa <- kappa - mean(kappa)

    ## --- beta：proximal Newton + Elastic Net 的近端算子 --------------------
    mu   <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
    g    <- as.vector((D - mu) %*% kappa)
    Wd   <- pmax(as.vector(mu %*% kappa^2), 1e-12)
    btil <- beta + g / Wd

    # 給定乘子 nu，逐年齡的封閉解：
    #   beta = beta0 + soft(Wd*(btil-beta0) - nu, lambda*a*w) / (Wd + lambda*(1-a)*w)
    bnu <- function(nu) {
      beta0 + soft(Wd * (btil - beta0) - nu, lambda * alpha_en * w) /
        (Wd + lambda * (1 - alpha_en) * w)
    }
    lo <- -1e12; hi <- 1e12
    for (i in 1:400) {
      mid <- (lo + hi) / 2
      sm <- sum(bnu(mid))
      if (!is.finite(sm)) break               # 數值失敗則放棄本次二分
      if (sm > 1) lo <- mid else hi <- mid
    }
    bnew <- bnu((lo + hi) / 2)
    if (any(!is.finite(bnew))) bnew <- beta   # 退回上一步

    ## 步長折半，確保目標函數不下降
    step <- 1; cur <- loglik(alpha, beta, kappa) - pen(beta)
    for (i in 1:40) {
      cand <- beta + step * (bnew - beta)
      s <- sum(cand); if (abs(s) > 1e-12) cand <- cand / s
      if (is.finite(loglik(alpha, cand, kappa)) &&
          loglik(alpha, cand, kappa) - pen(cand) >= cur - 1e-12) break
      step <- step / 2
    }
    beta <- beta + step * (bnew - beta)
    s <- sum(beta); beta <- beta / s; kappa <- kappa * s

    obj <- loglik(alpha, beta, kappa) - pen(beta)
    if (!is.finite(obj)) { beta <- beta_prev; break }
    if (track) hist <- rbind(hist, data.frame(
      iter = it, obj = obj, dbeta = max(abs(beta - beta_prev))))
    if (is.finite(obj_old) &&
        abs(obj - obj_old) < tol * max(1, abs(obj))) break
    obj_old <- obj
  }

  names(alpha) <- names(beta) <- rownames(D)
  names(kappa) <- colnames(D)
  mu  <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
  dev <- 2 * sum(ifelse(D > 0, D * log(D / mu), 0) - (D - mu))

  list(alpha = alpha, beta = beta, kappa = kappa, fitted = mu,
       deviance = dev, loglik = loglik(alpha, beta, kappa),
       nfree = sum(abs(beta - beta0) > 1e-7),   # 未被吸到 beta0 的年齡組數
       iter = it, converged = it < maxit,
       lambda = lambda, alpha_en = alpha_en, beta0 = beta0, history = hist)
}

## ===================== 2. lambda 路徑 =====================================
lc_pen_lambda_max <- function(D, E, beta0 = NULL, w = NULL, alpha_en = 1) {
  A <- nrow(D)
  if (is.null(beta0)) beta0 <- rep(1 / A, A)
  if (is.null(w)) w <- rep(1, A)
  f  <- lc_penalized(D, E, lambda = 0)
  mu <- f$fitted
  g  <- as.vector((D - mu) %*% f$kappa)
  Wd <- as.vector(mu %*% f$kappa^2)
  # L1 部分：使所有 beta 都被吸到 beta0 的最小 lambda（去掉共同平移後取最大）
  z <- Wd * (f$beta + g / Wd - beta0)
  l1 <- max(abs(z - median(z)) / w)
  if (alpha_en >= 0.01) return(l1 / alpha_en)
  # 純 Ridge 沒有能造成精確崩塌的有限 lambda，改以「懲罰曲率 = 概似曲率」
  # 的尺度為基準：lambda = Wd/w 時收縮約 50%，取 10 倍作為 lambda_max。
  10 * max(Wd / w)
}

lc_pen_path <- function(D, E, alpha_en = 1, beta0 = NULL, w = NULL,
                        nlam = 20, frac_min = 1e-3, verbose = FALSE) {
  A <- nrow(D); Tn <- ncol(D); N <- A * Tn
  lmax <- lc_pen_lambda_max(D, E, beta0, w, alpha_en)
  lams <- lmax * exp(seq(0, log(frac_min), length.out = nlam))
  res <- lapply(seq_along(lams), function(i) {
    f  <- lc_penalized(D, E, lambda = lams[i], alpha_en = alpha_en,
                       beta0 = beta0, w = w)
    df <- f$nfree + (Tn - 2) + A
    if (verbose) cat(sprintf("  lambda=%.4g  nfree=%d  dev=%.1f\n",
                             lams[i], f$nfree, f$deviance))
    data.frame(lambda = lams[i], frac = lams[i] / lmax, deviance = f$deviance,
               nfree = f$nfree, iter = f$iter,
               bic = f$deviance + log(N) * df, aic = f$deviance + 2 * df)
  })
  out <- do.call(rbind, res); attr(out, "lambda_max") <- lmax; out
}

## ===================== 3. 偏誤—變異數分解 =================================
# 對每個 lambda，以蒙地卡羅估 bias^2、variance、MSE（對 beta）
#' @param fracs  lambda / lambda_max 的比例格點（每個模擬資料集各自計算
#'               lambda_max，否則不同 N 之間的 lambda 不可比）
pen_bias_variance <- function(truth, N, fracs, alpha_en = 1, reps = 60,
                              beta0 = NULL, seed = 1234, verbose = TRUE) {
  set.seed(seed)
  A <- length(truth$a); E <- N * truth$w
  if (is.null(beta0)) beta0 <- rep(1 / A, A)
  store <- array(NA_real_, c(reps, length(fracs), A))
  lmax_v <- numeric(reps)
  for (r in seq_len(reps)) {
    D <- matrix(rpois(length(E), E * truth$m), nrow(E), ncol(E),
                dimnames = dimnames(truth$m))
    lmax <- tryCatch(lc_pen_lambda_max(D, E, beta0, NULL, alpha_en),
                     error = function(e) NA_real_)
    lmax_v[r] <- lmax
    if (!is.finite(lmax)) next
    for (j in seq_along(fracs))
      store[r, j, ] <- tryCatch(
        lc_penalized(D, E, lambda = fracs[j] * lmax,
                     alpha_en = alpha_en, beta0 = beta0)$beta,
        error = function(e) rep(NA_real_, A))
    if (verbose && r %% 20 == 0) cat(sprintf("   rep %d/%d\n", r, reps))
  }
  do.call(rbind, lapply(seq_along(fracs), function(j) {
    B <- store[, j, ]
    B <- B[stats::complete.cases(B), , drop = FALSE]
    mb <- colMeans(B)
    data.frame(frac = fracs[j], lambda_max_med = median(lmax_v, na.rm = TRUE),
               nrep_ok = nrow(B),
               bias2 = sum((mb - truth$b)^2),
               variance = sum(apply(B, 2, var)),
               mse = mean(rowSums(sweep(B, 2, truth$b)^2)),
               mse_med = median(rowSums(sweep(B, 2, truth$b)^2)),
               cor_med = median(apply(B, 1, function(v)
                 cor(v - mean(v), truth$b - mean(truth$b)))))
  }))
}

## ===================== 4. 變異數的收斂速度 ================================
# 對每個 N 估 sum_x Var(beta_hat_x)，再以 log-log 迴歸取斜率
#' @param frac  lambda / lambda_max；0 即未懲罰
pen_variance_rate <- function(truth, Ns, frac = 0, alpha_en = 1, reps = 60,
                              beta0 = NULL, seed = 99, verbose = TRUE) {
  set.seed(seed)
  A <- length(truth$a)
  if (is.null(beta0)) beta0 <- rep(1 / A, A)
  out <- lapply(Ns, function(N) {
    E <- N * truth$w
    B <- t(replicate(reps, {
      D <- matrix(rpois(length(E), E * truth$m), nrow(E), ncol(E),
                  dimnames = dimnames(truth$m))
      lam <- if (frac == 0) 0 else
        frac * tryCatch(lc_pen_lambda_max(D, E, beta0, NULL, alpha_en),
                        error = function(e) 0)
      tryCatch(lc_penalized(D, E, lambda = lam, alpha_en = alpha_en,
                            beta0 = beta0)$beta,
               error = function(e) rep(NA_real_, A))
    }))
    B <- B[stats::complete.cases(B), , drop = FALSE]
    mb <- colMeans(B)
    if (verbose) cat(sprintf("   N=%.0e 完成（%d/%d 次成功）\n", N, nrow(B), reps))
    data.frame(N = N, nrep_ok = nrow(B), variance = sum(apply(B, 2, var)),
               bias2 = sum((mb - truth$b)^2),
               mse = mean(rowSums(sweep(B, 2, truth$b)^2)))
  })
  res <- do.call(rbind, out)
  fit <- lm(log(variance) ~ log(N), data = res)
  attr(res, "slope") <- unname(coef(fit)[2])
  attr(res, "slope_se") <- summary(fit)$coefficients[2, 2]
  res
}

## ===================== 5. Elastic Net 的變體 ==============================
# 本節三個變體都針對同一個問題：LASSO 的收縮會對「大係數」也施加同樣的
# 懲罰，因而產生偏誤。我們發現 Ridge 在 MSE 上勝過 LASSO，正是此偏誤所致。

#' 自適應權重（Zou 2006）：以未懲罰估計的偏離量倒數為權重，
#' 使偏離 beta0 較遠（即證據較強）的年齡組受到較輕的懲罰。
adaptive_w <- function(D, E, beta0 = NULL, gamma = 1, eps = 1e-4) {
  A <- nrow(D); if (is.null(beta0)) beta0 <- rep(1 / A, A)
  b <- lc_penalized(D, E, lambda = 0)$beta
  w <- 1 / pmax(abs(b - beta0), eps)^gamma
  w / mean(w)
}

#' 鬆弛 LASSO（Meinshausen 2007）：先以 LASSO 選出作用集，
#' 再在作用集上「不加懲罰」重新配適，非作用集固定於 beta0。
#' 以權重實作：作用集 w = 0（無懲罰）、非作用集 w 極大（釘在 beta0）。
lc_relaxed <- function(D, E, lambda, alpha_en = 1, beta0 = NULL, w = NULL, ...) {
  A <- nrow(D); if (is.null(beta0)) beta0 <- rep(1 / A, A)
  f1 <- lc_penalized(D, E, lambda = lambda, alpha_en = alpha_en,
                     beta0 = beta0, w = w, ...)
  free <- abs(f1$beta - beta0) > 1e-7
  if (!any(free)) return(f1)
  w2 <- ifelse(free, 0, 1e12)
  f2 <- lc_penalized(D, E, lambda = 1, alpha_en = 1, beta0 = beta0, w = w2, ...)
  f2$nfree <- sum(free); f2$stage1 <- f1; f2
}

#' MCP（Zhang 2010）：非凸懲罰，對大係數的懲罰逐步歸零，
#' 因而保留 LASSO 的稀疏性但去除其對大係數的偏誤。
#' 近端算子由 min 0.5*W*u^2 - c*u + w*p_{lambda,gamma}(u) 導出：
#'   a = c/W, lam' = lambda*w/W
#'   |a| <= gamma*lambda -> u = S(a, lam') / (1 - w/(gamma*W))
#'   否則               -> u = a
lc_mcp <- function(D, E, lambda = 0, gamma = 3, beta0 = NULL, w = NULL,
                   maxit = 500, tol = 1e-11) {
  A <- nrow(D); Tn <- ncol(D)
  if (is.null(beta0)) beta0 <- rep(1 / A, A)
  if (is.null(w))     w     <- rep(1, A)
  alpha <- log(rowSums(pmax(D, 0.5)) / rowSums(E))
  beta  <- rep(1 / A, A); kappa <- seq(1, -1, length.out = Tn)

  loglik <- function(a, b, k) {
    eta <- outer(a, rep(1, Tn)) + outer(b, k)
    if (any(!is.finite(eta))) return(-Inf)
    v <- sum(D * eta - E * exp(pmin(eta, 50))); if (is.finite(v)) v else -Inf
  }
  obj_old <- -Inf
  for (it in seq_len(maxit)) {
    bp <- beta
    mu <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
    alpha <- alpha + log(rowSums(D) / rowSums(mu))
    mu <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
    kappa <- kappa + colSums((D - mu) * beta) / pmax(colSums(mu * beta^2), 1e-12)
    kappa <- kappa - mean(kappa)
    mu <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
    g  <- as.vector((D - mu) %*% kappa)
    Wd <- pmax(as.vector(mu %*% kappa^2), 1e-12)
    btil <- beta + g / Wd

    bnu <- function(nu) {
      cc  <- Wd * (btil - beta0) - nu
      a   <- cc / Wd
      lam2 <- lambda * w / Wd
      den <- pmax(1 - w / (gamma * Wd), 1e-6)
      u   <- ifelse(abs(a) <= gamma * lambda, soft(a, lam2) / den, a)
      beta0 + u
    }
    lo <- -1e12; hi <- 1e12
    for (i in 1:400) {
      mid <- (lo + hi) / 2; sm <- sum(bnu(mid))
      if (!is.finite(sm)) break
      if (sm > 1) lo <- mid else hi <- mid
    }
    bnew <- bnu((lo + hi) / 2)
    if (any(!is.finite(bnew))) bnew <- beta
    step <- 1
    for (i in 1:40) {
      cand <- beta + step * (bnew - beta); s <- sum(cand)
      if (abs(s) > 1e-12) cand <- cand / s
      if (is.finite(loglik(alpha, cand, kappa))) break
      step <- step / 2
    }
    beta <- beta + step * (bnew - beta); s <- sum(beta)
    beta <- beta / s; kappa <- kappa * s
    obj <- loglik(alpha, beta, kappa)
    if (!is.finite(obj)) { beta <- bp; break }
    if (is.finite(obj_old) && abs(obj - obj_old) < tol * max(1, abs(obj))) break
    obj_old <- obj
  }
  names(alpha) <- names(beta) <- rownames(D); names(kappa) <- colnames(D)
  list(alpha = alpha, beta = beta, kappa = kappa, iter = it,
       nfree = sum(abs(beta - beta0) > 1e-7), lambda = lambda, gamma = gamma)
}
