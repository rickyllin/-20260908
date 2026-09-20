###############################################################################
# Rabbi & Mazzuco (2021, EJP 37:97-120) 方法的復刻
#   LC_e0dagger = (1) 先用 2D L1（quantile LASSO）修勻死亡率曲面
#                 (2) 在修勻後的資料上跑標準 LC (SVD)
#                 (3) 調整 kappa_t 使「配適的壽命差異 e0dagger」= 觀察值
#                 (4) ARIMA(0,1,0) with drift 外推
#
# 對照組：LC（k 調整至總死亡數）、LM（Lee-Miller，k 調整至 e0）、LCP（Poisson）
#
# 與原文的差異（務必在報告裡寫明）：
#   * 原文用單齡組 0-100+，本檔是五齡組 22 組。年齡方向的二階差分改用
#     「不等距分割差商」，以組中點為節點。
#   * e0dagger 改用簡略生命表（Chiang）的離散近似，精度低於單齡版本。
#   * 原文用 smoothAPC 套件；本檔直接把式(4)(5) 的堆疊設計矩陣丟給
#     quantreg::rq.fit.sfn，結果相同但可控、且不需要 Nelder-Mead。
#
# 需要：readxl, quantreg, SparseM, Matrix
###############################################################################

source("lc_poisson_lasso.R")     # load_data(), lc_poisson()
library(Matrix)

## ========================= 1. 年齡格點與差分算子 ===========================
make_grid <- function(A) {
  n     <- c(1, 4, rep(5, A - 2))
  start <- c(0, cumsum(n)[-A])
  mid   <- start + c(0.5, 2.0, rep(2.5, A - 2))
  list(n = n, start = start, mid = mid)
}

# 不等距二階差商
D2_uneven <- function(v) {
  n <- length(v); R <- matrix(0, n - 2, n)
  for (i in 2:(n - 1)) {
    h1 <- v[i] - v[i-1]; h2 <- v[i+1] - v[i]
    R[i-1, i-1] <-  2 / (h1 * (h1 + h2))
    R[i-1, i  ] <- -2 / (h1 * h2)
    R[i-1, i+1] <-  2 / (h2 * (h1 + h2))
  }
  R
}
D1_uneven <- function(v) {
  n <- length(v); R <- matrix(0, n - 1, n)
  for (i in 1:(n - 1)) { h <- v[i+1] - v[i]; R[i, i] <- -1/h; R[i, i+1] <- 1/h }
  R
}

## ========================= 2. 2D L1 修勻 ===================================
# min_z ||y - z||_1 + lxx||Dxx z||_1 + lxt||Dxt z||_1 + ltt||Dtt z||_1
# 疊成 ||y_ext - R z||_1，即 tau=0.5 的分位數迴歸（原文式 (4)(5)）
build_R <- function(A, Tn, mid, lxx, lxt, ltt) {
  Ia <- Diagonal(A); It <- Diagonal(Tn)
  Da2 <- Matrix(D2_uneven(mid), sparse = TRUE)
  Dt2 <- Matrix(D2_uneven(seq_len(Tn) - 1), sparse = TRUE)
  Da1 <- Matrix(D1_uneven(mid), sparse = TRUE)
  Dt1 <- Matrix(D1_uneven(seq_len(Tn) - 1), sparse = TRUE)
  rbind(Diagonal(A * Tn),
        lxx * kronecker(Da2, It),
        ltt * kronecker(Ia,  Dt2),
        lxt * kronecker(Da1, Dt1))
}

smooth_l1_2d <- function(logm, lxx = 5, lxt = 1, ltt = 5,
                         method = c("rq", "irls"), maxit = 60, eps = 1e-4) {
  method <- match.arg(method)
  A <- nrow(logm); Tn <- ncol(logm)
  mid <- make_grid(A)$mid
  R <- build_R(A, Tn, mid, lxx, lxt, ltt)
  # z 以「年齡慢、年份快」展開，與 kronecker 的順序一致
  y <- as.vector(t(logm))
  yext <- c(y, rep(0, nrow(R) - length(y)))

  if (method == "rq") {
    if (!requireNamespace("quantreg", quietly = TRUE) ||
        !requireNamespace("SparseM", quietly = TRUE))
      stop("需要 quantreg 與 SparseM；或改用 method = \"irls\"")
    Rc  <- methods::as(methods::as(R, "dgCMatrix"), "matrix.csr")
    fit <- quantreg::rq.fit.sfn(Rc, yext, tau = 0.5)
    z <- fit$coef
  } else {
    ## IRLS 近似 L1：以 1/max(|r|,eps) 為權重反覆解加權最小平方
    z <- y
    for (it in seq_len(maxit)) {
      r <- as.vector(yext - R %*% z)
      w <- 1 / pmax(abs(r), eps)
      W <- Diagonal(x = w)
      znew <- as.vector(solve(crossprod(R, W %*% R) + 1e-10 * Diagonal(ncol(R)),
                              crossprod(R, W %*% yext)))
      if (max(abs(znew - z)) < 1e-8) { z <- znew; break }
      z <- znew
    }
  }
  out <- matrix(z, nrow = A, ncol = Tn, byrow = TRUE, dimnames = dimnames(logm))
  out
}

## ========================= 3. 簡略生命表與壽命差異 =========================
life_table <- function(m, n = NULL) {
  A <- length(m)
  if (is.null(n)) n <- make_grid(A)$n
  a <- n / 2; a[1] <- 0.1; a[2] <- 1.6
  a[A] <- 1 / max(m[A], 1e-8)                       # 開放組
  q <- n * m / (1 + (n - a) * m)
  q[-A] <- pmin(q[-A], 1); q[A] <- 1
  l <- numeric(A + 1); l[1] <- 1
  for (i in 1:A) l[i+1] <- l[i] * (1 - q[i])
  d <- l[-(A+1)] - l[-1]
  L <- n * l[-1] + a * d; L[A] <- l[A] * a[A]
  Tt <- rev(cumsum(rev(L)))
  e  <- Tt / pmax(l[-(A+1)], 1e-300)
  e_end <- c(e[-1], 0)
  # 死於 [x, x+n) 者平均損失 = 區間內剩餘 (1-a/n)*n 年 + 區間末的餘命
  edag <- sum(d * ((1 - a / n) * n + e_end)) / l[1]
  list(e0 = e[1], edag = edag, e = e, lx = l, dx = d, qx = q)
}

## ========================= 4. LC 的各種 kappa 調整 =========================
lc_svd_fit <- function(logm) {
  a  <- rowMeans(logm); Z <- logm - a
  sv <- svd(Z); u1 <- sv$u[,1]; v1 <- sv$v[,1]
  if (sum(u1) < 0) { u1 <- -u1; v1 <- -v1 }
  list(a = a, b = u1 / sum(u1), k = sv$d[1] * v1 * sum(u1),
       varexp = sv$d[1]^2 / sum(sv$d^2))
}

# 通用：逐年解 kappa_t 使某個目標量吻合觀察值
adjust_kappa <- function(a, b, k0, target_fun, target_obs, span = 60) {
  vapply(seq_along(k0), function(j) {
    f <- function(kk) target_fun(exp(a + b * kk)) - target_obs[j]
    tryCatch(uniroot(f, c(k0[j] - span, k0[j] + span),
                     extendInt = "yes", tol = 1e-10)$root,
             error = function(e) NA_real_)
  }, numeric(1))
}

fit_variant <- function(D, E, variant = c("LC", "LM", "LCP", "LCedag"),
                        lxx = 5, lxt = 1, ltt = 5, smooth_method = "rq") {
  variant <- match.arg(variant)
  mobs <- D / E
  logm_obs <- log(mobs)

  if (variant %in% c("LC", "LM", "LCP")) {
    logm_use <- logm_obs
  } else {
    logm_use <- smooth_l1_2d(logm_obs, lxx, lxt, ltt, method = smooth_method)
  }

  if (variant == "LCP") {
    f <- lc_poisson(D, E, lambda = 0)
    return(list(a = f$alpha, b = f$beta, k = f$kappa, logm_use = logm_use))
  }

  f <- lc_svd_fit(logm_use)
  a <- f$a; b <- f$b; k <- f$k

  if (variant == "LC") {                       # 原始 LC：對齊總死亡數
    k <- vapply(seq_along(k), function(j) {
      Dt <- sum(D[, j])
      uniroot(function(kk) sum(E[, j] * exp(a + b * kk)) - Dt,
              c(k[j] - 60, k[j] + 60), extendInt = "yes", tol = 1e-10)$root
    }, numeric(1))
  } else if (variant == "LM") {                # Lee-Miller：對齊 e0
    e0_obs <- apply(mobs, 2, function(m) life_table(m)$e0)
    k <- adjust_kappa(a, b, k, function(m) life_table(m)$e0, e0_obs)
  } else if (variant == "LCedag") {            # 本文復刻：對齊 e0dagger
    ed_obs <- apply(mobs, 2, function(m) life_table(m)$edag)
    k <- adjust_kappa(a, b, k, function(m) life_table(m)$edag, ed_obs)
  }
  names(k) <- colnames(D)
  list(a = a, b = b, k = k, logm_use = logm_use)
}

## ========================= 5. 外推與樣本外評估 =============================
# ARIMA(0,1,0) with drift
forecast_k <- function(k, h) {
  n <- length(k); drift <- (k[n] - k[1]) / (n - 1)
  k[n] + drift * seq_len(h)
}

# jump-off 修正：LM 與 LCedag 用實際資料當起點（原文 3.2 節）
oos_evaluate <- function(D, E, h = 10, variants = c("LC","LM","LCP","LCedag"),
                         lxx = 5, lxt = 1, ltt = 5, smooth_method = "rq",
                         jumpoff = c(LC = FALSE, LM = TRUE,
                                     LCP = FALSE, LCedag = TRUE)) {
  Tn  <- ncol(D); fitT <- Tn - h
  Df  <- D[, 1:fitT, drop = FALSE]; Ef <- E[, 1:fitT, drop = FALSE]
  Dh  <- D[, (fitT+1):Tn, drop = FALSE]; Eh <- E[, (fitT+1):Tn, drop = FALSE]
  logm_hold <- log(Dh / Eh)
  e0_hold   <- apply(Dh / Eh, 2, function(m) life_table(m)$e0)
  ed_hold   <- apply(Dh / Eh, 2, function(m) life_table(m)$edag)

  res <- lapply(variants, function(v) {
    f <- fit_variant(Df, Ef, v, lxx, lxt, ltt, smooth_method)
    a <- f$a
    if (isTRUE(jumpoff[[v]]))                 # 以最後一年實際值為起點
      a <- log(Df[, fitT] / Ef[, fitT]) - f$b * f$k[fitT]
    kf   <- forecast_k(f$k, h)
    pred <- outer(a, rep(1, h)) + outer(f$b, kf)
    e0p  <- apply(exp(pred), 2, function(m) life_table(m)$e0)
    edp  <- apply(exp(pred), 2, function(m) life_table(m)$edag)
    err  <- pred - logm_hold
    data.frame(method = v,
               MAE  = mean(abs(err)), MSE = mean(err^2),
               MAE_h1  = mean(abs(err[, 1])),  MAE_h10 = mean(abs(err[, h])),
               ME_e0   = mean(e0_hold - e0p),        # 正 = 低估存活改善
               MAE_e0  = mean(abs(e0_hold - e0p)),
               MAE_edag = mean(abs(ed_hold - edp)),
               drift = (f$k[fitT] - f$k[1]) / (fitT - 1),
               k_lin_resid = sd(residuals(lm(f$k ~ seq_len(fitT)))))
  })
  do.call(rbind, res)
}

## ========================= 6. 執行 =========================================
if (sys.nframe() == 0) {

  dat <- load_data(sex = "Total")      # 亦可 "Female" / "Male"
  D <- dat$D; E <- dat$E
  cat(sprintf("資料：%d 個五齡組 x %d 年 (%d-%d)\n",
              nrow(D), ncol(D), min(dat$years), max(dat$years)))

  ## --- (a) 觀察到的 e0 與 e0dagger ---
  mobs <- D / E
  e0   <- apply(mobs, 2, function(m) life_table(m)$e0)
  ed   <- apply(mobs, 2, function(m) life_table(m)$edag)
  cat(sprintf("e0:  %.2f (%d) -> %.2f (%d)\n", e0[1], dat$years[1],
              e0[length(e0)], max(dat$years)))
  cat(sprintf("e0†: %.2f -> %.2f   (corr with e0 = %.3f)\n",
              ed[1], ed[length(ed)], cor(e0, ed)))

  ## --- (b) 修勻前後 ---
  logm <- log(mobs)
  sm   <- smooth_l1_2d(logm, lxx = 5, lxt = 1, ltt = 5, method = "rq")
  mid  <- make_grid(nrow(D))$mid
  ra   <- function(M) mean(abs(D2_uneven(mid) %*% M))
  rt   <- function(M) mean(abs(M %*% t(D2_uneven(seq_len(ncol(M)) - 1))))
  cat(sprintf("\n修勻：MAE(vs 觀察)=%.5f；年齡方向粗糙度 %.4f->%.4f；時間方向 %.4f->%.4f\n",
              mean(abs(sm - logm)), ra(logm), ra(sm), rt(logm), rt(sm)))

  ## --- (c) 四種變體的 kappa ---
  fits <- lapply(c("LC","LM","LCP","LCedag"),
                 function(v) fit_variant(D, E, v))
  names(fits) <- c("LC","LM","LCP","LCedag")
  cat("\n各變體的 kappa 性質：\n")
  for (v in names(fits)) {
    k <- fits[[v]]$k; n <- length(k)
    cat(sprintf("  %-7s drift=%.4f  離線性殘差 sd=%.3f  sd(diff)=%.3f\n",
                v, (k[n]-k[1])/(n-1),
                sd(residuals(lm(k ~ seq_len(n)))), sd(diff(k))))
  }

  ## --- (d) 樣本外評估：留最後 10 年 ---
  cat("\n樣本外評估（留最後 10 年）：\n")
  oos <- oos_evaluate(D, E, h = 10)
  print(oos, digits = 4, row.names = FALSE)
  write.csv(oos, "rabbi_mazzuco_oos.csv", row.names = FALSE)

  ## --- 圖 ---
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  matplot(dat$years, t(logm), type = "l", col = rgb(0,0,0,.15), lty = 1,
          xlab = "Year", ylab = "log m", main = "觀察值")
  matplot(dat$years, t(sm), type = "l", col = rgb(0,0,1,.3), lty = 1,
          xlab = "Year", ylab = "log m", main = "L1 修勻後")
  plot(dat$years, fits$LC$k, type = "l", lwd = 2,
       ylim = range(sapply(fits, `[[`, "k")),
       xlab = "Year", ylab = expression(kappa[t]), main = "kappa_t 各調整方式")
  cols <- c(LC="black", LM="blue", LCP="grey50", LCedag="red")
  for (v in names(fits)) lines(dat$years, fits[[v]]$k, col = cols[v], lwd = 2)
  legend("topright", names(cols), col = cols, lwd = 2, bty = "n", cex = .8)
  plot(e0, ed, type = "b", pch = 16, xlab = expression(e[0]),
       ylab = expression(e[0]^"†"),
       main = sprintf("e0 vs e0† (r = %.3f)", cor(e0, ed)))
  par(op)
}

###############################################################################
# 我在你的資料上先跑過的結果（Total, 1970-2024, 五齡組）
#
#   e0   68.73 -> 80.55 ; e0†  14.70 -> 11.83（反向關係成立，與原文一致）
#
#   kappa 的離線性殘差 sd：  原始 SVD 0.849 / Lee-Miller(e0) 0.772 / e0† 1.420
#   kappa 的 sd(diff)     ：  0.582 / 0.575 / 0.890
#
# 亦即：**e0† 調整讓 kappa 變得比原始 LC 更不線性、更抖**，與原文在瑞典
# 資料上「k 更規則、標準誤更小」的結論相反。兩個可能的原因，值得在報告裡
# 拆開檢驗：
#   (1) 五齡組 + 100+ 開放組使 e0† 的離散近似誤差偏大，而 e0† 本來就比 e0
#       更依賴死亡年齡分布的細節；
#   (2) 台灣 1970-2024 的死亡轉型比瑞典 1950-2016 快且不規則。
# 用單齡組資料重跑一次即可分離這兩者——這正是原文沒做、而你可以補的檢驗。
###############################################################################
