###############################################################################
# lambda 的選擇
#
#   比較四種選法，並與「已知真值」的 oracle 對照：
#     (A) 原文式：對觀察值的 MAE/MSE          -> 必然選到 lambda = 0
#     (B) 留出格交叉驗證（Poisson 離差）      -> 系統性修勻不足
#     (C) 留出年份的下游預測誤差              -> 瞄準預測目標
#     (D) 參數式 bootstrap 的插入式選法       -> 瞄準 beta 的 MSE（推薦）
#
# 需先 source("R/studies/mc_exposure.R")（其中已 source 主檔與修正檔）
###############################################################################

source("R/studies/mc_exposure.R")

## ===================== 0. 可加權的 L1 修勻 =================================
# wt = 0 的格子不進入保真項，即為「留出」
smooth_l1_w <- function(logm, l, wt = NULL) {
  A <- nrow(logm); Tn <- ncol(logm); mid <- make_grid(A)$mid
  Da2 <- Matrix(D2_uneven(mid), sparse = TRUE)
  Dt2 <- Matrix(D2_uneven(seq_len(Tn) - 1), sparse = TRUE)
  Da1 <- Matrix(D1_uneven(mid), sparse = TRUE)
  Dt1 <- Matrix(D1_uneven(seq_len(Tn) - 1), sparse = TRUE)
  
  fid <- Diagonal(A * Tn)
  yv  <- as.vector(t(logm))
  if (!is.null(wt)) { wv <- as.vector(t(wt)); fid <- Diagonal(x = wv); yv <- yv * wv }
  
  R <- rbind(fid,
             l       * kronecker(Da2, Diagonal(Tn)),
             l       * kronecker(Diagonal(A), Dt2),
             (l / 5) * kronecker(Da1, Dt1))
  yext <- c(yv, rep(0, nrow(R) - length(yv)))
  Rc <- methods::as(methods::as(R, "dgCMatrix"), "matrix.csr")
  z  <- rq_sfn_safe(Rc, yext, tau = 0.5)$coef
  matrix(z, A, Tn, byrow = TRUE, dimnames = dimnames(logm))
}

poisson_dev <- function(Dobs, Eobs, z) {
  mu <- pmax(Eobs * exp(z), 1e-10)
  2 * sum(ifelse(Dobs > 0, Dobs * log(Dobs / mu), 0) - (Dobs - mu))
}

## ===================== 1. (A) 原文式準則 ===================================
crit_fit_error <- function(D, E, lambdas, zero_sub = 0.5) {
  logm <- log(pmax(D, zero_sub) / E)
  data.frame(lambda = lambdas,
             MAE = vapply(lambdas, function(l)
               mean(abs(smooth_l1_w(logm, l) - logm)), numeric(1)),
             MSE = vapply(lambdas, function(l)
               mean((smooth_l1_w(logm, l) - logm)^2), numeric(1)))
}

## ===================== 2. (B) 留出格交叉驗證 ===============================
# 注意：lambda = 0 無法被交叉驗證（對留出格產生不了預測），
#       故此準則只能在 lambda > 0 之間排序。
select_lambda_cv <- function(D, E, lambdas, K = 3, zero_sub = 0.5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  stopifnot(all(lambdas > 0))
  logm <- log(pmax(D, zero_sub) / E)
  fold <- matrix(sample(rep_len(seq_len(K), length(D))), nrow(D), ncol(D))
  sc <- vapply(lambdas, function(l) {
    s <- 0
    for (f in seq_len(K)) {
      wt <- (fold != f) * 1
      Z  <- smooth_l1_w(logm, l, wt)
      m_ <- fold == f
      s  <- s + poisson_dev(D[m_], E[m_], Z[m_])
    }
    s
  }, numeric(1))
  list(lambda = lambdas[which.min(sc)],
       table = data.frame(lambda = lambdas, cv_deviance = sc))
}

## ===================== 3. (C) 留出年份的下游預測誤差 =======================
select_lambda_forecast <- function(D, E, lambdas, h = 5, zero_sub = 0.5) {
  Tn <- ncol(D); fT <- Tn - h
  Df <- D[, 1:fT, drop = FALSE]; Ef <- E[, 1:fT, drop = FALSE]
  lmh <- log(pmax(D[, (fT+1):Tn, drop = FALSE], zero_sub) /
               E[, (fT+1):Tn, drop = FALSE])
  logf <- log(pmax(Df, zero_sub) / Ef)
  sc <- vapply(lambdas, function(l) {
    Z <- if (l == 0) logf else smooth_l1_w(logf, l)
    f <- lc_svd_fit(Z)
    dr <- (f$k[fT] - f$k[1]) / (fT - 1)
    pred <- outer(f$a, rep(1, h)) + outer(f$b, f$k[fT] + dr * seq_len(h))
    mean(abs(pred - lmh))
  }, numeric(1))
  list(lambda = lambdas[which.min(sc)],
       table = data.frame(lambda = lambdas, fc_MAE = sc))
}

## ===================== 4. (D) 參數式 bootstrap 插入式選法 ==================
# 直接瞄準估計量（預設 beta）的 MSE：
#   以 pilot lambda 配適得到 m_hat -> 自 m_hat 模擬 B 組 -> 對每個 lambda
#   計算估計量相對於 pilot 值的 SSE -> 取最小。可迭代以降低 pilot 的偏誤。
select_lambda_boot <- function(D, E, lambdas, pilot = 5, B = 20,
                               target = c("beta", "drift", "alpha"),
                               iter = 1, zero_sub = 0.5, seed = NULL) {
  target <- match.arg(target)
  if (!is.null(seed)) set.seed(seed)
  logm <- log(pmax(D, zero_sub) / E); Tn <- ncol(D)
  get_t <- function(f) switch(target,
                              beta = f$b, alpha = f$a, drift = (f$k[Tn] - f$k[1]) / (Tn - 1))
  
  lam <- pilot; tab <- NULL
  for (it in seq_len(iter)) {
    f0 <- lc_svd_fit(smooth_l1_w(logm, lam))
    m0 <- exp(outer(f0$a, rep(1, Tn)) + outer(f0$b, f0$k))
    t0 <- get_t(f0)
    sse <- matrix(NA_real_, B, length(lambdas))
    for (bb in seq_len(B)) {
      Db <- matrix(rpois(length(E), E * m0), nrow(E), ncol(E))
      lb <- log(pmax(Db, zero_sub) / E)
      for (j in seq_along(lambdas))
        sse[bb, j] <- sum((get_t(lc_svd_fit(smooth_l1_w(lb, lambdas[j]))) - t0)^2)
    }
    sc  <- apply(sse, 2, median)
    lam <- lambdas[which.min(sc)]
    tab <- data.frame(lambda = lambdas, boot_sse = sc)
  }
  list(lambda = lam, table = tab)
}

## ===================== 5. 主研究：各選法 vs oracle =========================
lambda_study <- function(truth, Ns = c(5e4, 2e5, 1e6),
                         lambdas = c(1, 2, 5, 10, 20, 40),
                         reps = 20, B = 15, seed = 5, verbose = TRUE) {
  set.seed(seed)
  Tn <- length(truth$k)
  drift_true <- (truth$k[Tn] - truth$k[1]) / (Tn - 1)
  res <- list()
  
  for (N in Ns) {
    E <- N * truth$w
    sse_b <- sse_a <- matrix(NA_real_, reps, length(lambdas))
    pick  <- data.frame(cv = numeric(reps), fc = numeric(reps), bt = numeric(reps))
    sse_pick <- pick
    
    for (r in seq_len(reps)) {
      D <- matrix(rpois(length(E), E * truth$m), nrow(E), ncol(E),
                  dimnames = dimnames(truth$m))
      logm <- log(pmax(D, 0.5) / E)
      for (j in seq_along(lambdas)) {
        f <- lc_svd_fit(smooth_l1_w(logm, lambdas[j]))
        sse_b[r, j] <- sum((f$b - truth$b)^2)
        sse_a[r, j] <- sum((f$a - truth$a)^2)
      }
      pick$cv[r] <- select_lambda_cv(D, E, lambdas)$lambda
      pick$fc[r] <- select_lambda_forecast(D, E, lambdas)$lambda
      pick$bt[r] <- select_lambda_boot(D, E, lambdas, B = B)$lambda
      for (nm in names(pick))
        sse_pick[[nm]][r] <- sse_b[r, match(pick[[nm]][r], lambdas)]
      if (verbose && r %% 5 == 0) cat(sprintf("   N=%.0e rep %d\n", N, r))
    }
    
    med_b <- apply(sse_b, 2, median)
    orc   <- lambdas[which.min(med_b)]
    res[[length(res)+1]] <- data.frame(
      N = N,
      oracle_lambda = orc, oracle_sse = min(med_b),
      oracle_lambda_alpha = lambdas[which.min(apply(sse_a, 2, median))],
      cv_lambda_mode  = as.numeric(names(which.max(table(pick$cv)))),
      cv_sse  = median(sse_pick$cv),  cv_regret  = median(sse_pick$cv)  / min(med_b),
      fc_lambda_mode  = as.numeric(names(which.max(table(pick$fc)))),
      fc_sse  = median(sse_pick$fc),  fc_regret  = median(sse_pick$fc)  / min(med_b),
      bt_lambda_mode  = as.numeric(names(which.max(table(pick$bt)))),
      bt_sse  = median(sse_pick$bt),  bt_regret  = median(sse_pick$bt)  / min(med_b))
    if (verbose) cat(sprintf("N = %.0e 完成\n", N))
  }
  out <- do.call(rbind, res)
  attr(out, "drift_true") <- drift_true
  out
}

## ===================== 6. 執行 =============================================
if (sys.nframe() == 0) {
  truth <- mc_truth("Female", 2001, 2024)
  lams  <- c(1, 2, 5, 10, 20, 40)
  
  ## --- (A) 原文式準則：在真實資料上示範它會選到哪裡 ---
  dat <- load_data(sex = "Female")
  keep <- which(dat$years >= 2001 & dat$years <= 2024)
  Dr <- dat$D[, keep]; Er <- dat$E[, keep]
  cat("原文式準則（對觀察值的 MAE/MSE）：\n")
  print(crit_fit_error(Dr, Er, c(0, lams)), digits = 4, row.names = FALSE)
  cat("→ 兩欄都在 lambda = 0 取最小；此準則無法選出非零的修勻程度。\n\n")
  
  ## --- 主研究（耗時，先用小 reps 試跑）---
  cat("各選法 vs oracle（reps 建議 >= 30，先用 10 試流程）：\n")
  st <- lambda_study(truth, reps = 10, B = 12)
  print(st, digits = 4, row.names = FALSE)
  write.csv(st, "output/tables/lambda_selection.csv", row.names = FALSE)
  
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  matplot(st$N, cbind(st$oracle_lambda, st$cv_lambda_mode,
                      st$fc_lambda_mode, st$bt_lambda_mode),
          type = "b", pch = 16, log = "xy", lty = 1,
          xlab = "Exposure N", ylab = expression(lambda),
          main = "選到的 lambda")
  legend("topright", c("oracle","CV","forecast","bootstrap"),
         col = 1:4, pch = 16, bty = "n", cex = .8)
  matplot(st$N, cbind(st$cv_regret, st$fc_regret, st$bt_regret),
          type = "b", pch = 16, log = "x", lty = 1, col = 2:4,
          xlab = "Exposure N", ylab = "SSE(beta) / oracle",
          main = "相對損失（1 = 追平 oracle）")
  abline(h = 1, lty = 2)
  par(op)
}

###############################################################################
# 先跑過的結果（Female 2001-2024 為真值）
#
#  N        oracle lambda   oracle SSE   CV 選到   CV 的 SSE   相對損失
#  5e4           40           0.0238      0,2,5     0.1701       x7.2
#  2e5           20           0.0057      2         0.0179       x3.1
#  1e6            2           0.0037      2         0.0037       x1.0
#
#  bootstrap 插入式選法：
#  5e4        選 10/20/40    SSE 0.0843   x3.5
#  2e5        選 10/20/40    SSE 0.0046   x0.80   （勝過格點 oracle）
#  1e6        選 2/5/10/20   SSE 0.0032   x0.85
#
# 四個結論：
#  1. 原文式準則（對觀察值的 MAE/MSE）必然選到 lambda = 0，因為修勻一定會
#     增加對觀察值的誤差。這不是方法比較，是準則的定義決定的。
#  2. CV 在大人口完美命中，在小人口修勻不足七倍。原因不是 CV 壞掉，而是
#     CV 衡量「曲面」的預測誤差，而 beta 是曲面的泛函、對雜訊敏感得多，
#     需要的修勻程度遠高於曲面本身。落差在最需要它的地方最大。
#  3. lambda = 0 無法被交叉驗證（對留出格產生不了預測），CV 只能在
#     lambda > 0 之間排序，回答不了「要不要修勻」。
#  4. alpha 與 beta 的最適 lambda 方向相反（alpha 恆在 lambda 約 1-2 最小，
#     lambda = 40 時 SSE 暴增到 7.6-8.7）。曲面上只有一個旋鈕，必然妥協；
#     這是「懲罰應加在參數而非曲面」的最強論據。
###############################################################################