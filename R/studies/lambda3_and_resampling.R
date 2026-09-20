###############################################################################
# (1) 三個 lambda 方向分開研究
# (2) 重抽樣模組：參數不確定性與預測不確定性
#
# 需先 source("R/studies/mc_exposure.R")
###############################################################################

source("R/studies/mc_exposure.R")

###############################################################################
# 第一部分：三個方向分開的 lambda
#
#   min ||y - z||_1 + lxx||Dxx z||_1 + lxt||Dxt z||_1 + ltt||Dtt z||_1
#
#   Dxx 罰年齡方向的曲率、Dtt 罰時間方向的曲率、
#   Dxt 罰「年齡型態隨時間的變化」。在 LC 底下
#       d^2 log m / dx dt = beta'_x * kappa'_t
#   所以 Dxt 本質上就是對 beta_x 粗糙度的懲罰。
###############################################################################

smooth3 <- function(logm, lxx = 0, lxt = 0, ltt = 0) {
  A <- nrow(logm); Tn <- ncol(logm); mid <- make_grid(A)$mid
  blocks <- list(Diagonal(A * Tn))
  if (lxx > 0) blocks <- c(blocks, lxx * kronecker(
    Matrix(D2_uneven(mid), sparse = TRUE), Diagonal(Tn)))
  if (ltt > 0) blocks <- c(blocks, ltt * kronecker(
    Diagonal(A), Matrix(D2_uneven(seq_len(Tn) - 1), sparse = TRUE)))
  if (lxt > 0) blocks <- c(blocks, lxt * kronecker(
    Matrix(D1_uneven(mid), sparse = TRUE),
    Matrix(D1_uneven(seq_len(Tn) - 1), sparse = TRUE)))
  if (length(blocks) == 1L) return(logm)          # 三個都是 0 = 不修勻
  
  R <- do.call(rbind, blocks)
  y <- as.vector(t(logm))
  yext <- c(y, rep(0, nrow(R) - length(y)))
  Rc <- methods::as(methods::as(R, "dgCMatrix"), "matrix.csr")
  z  <- rq_sfn_safe(Rc, yext, tau = 0.5)$coef
  matrix(z, A, Tn, byrow = TRUE, dimnames = dimnames(logm))
}

# 三方向格點搜尋（蒙地卡羅，需要真值）
lambda_grid3 <- function(truth, N,
                         lxx = c(0, 1, 5, 20), lxt = c(0, 2, 5),
                         ltt = c(0, 1, 5, 20),
                         reps = 10, seed = 31, verbose = TRUE) {
  set.seed(seed)
  E <- N * truth$w; Tn <- length(truth$k)
  drift_true <- (truth$k[Tn] - truth$k[1]) / (Tn - 1)
  grid <- expand.grid(lxx = lxx, lxt = lxt, ltt = ltt)
  acc <- lapply(seq_len(nrow(grid)), function(i)
    list(b = numeric(reps), a = numeric(reps), d = numeric(reps)))
  
  for (r in seq_len(reps)) {
    D <- matrix(rpois(length(E), E * truth$m), nrow(E), ncol(E),
                dimnames = dimnames(truth$m))
    logm <- log(pmax(D, 0.5) / E)
    for (i in seq_len(nrow(grid))) {
      Z <- smooth3(logm, grid$lxx[i], grid$lxt[i], grid$ltt[i])
      f <- lc_svd_fit(Z)
      acc[[i]]$b[r] <- sum((f$b - truth$b)^2)
      acc[[i]]$a[r] <- sum((f$a - truth$a)^2)
      acc[[i]]$d[r] <- ((f$k[Tn] - f$k[1]) / (Tn - 1) - drift_true)^2
    }
    if (verbose) cat(sprintf("  rep %d/%d\n", r, reps))
  }
  cbind(grid,
        sse_b = vapply(acc, function(x) median(x$b), numeric(1)),
        sse_a = vapply(acc, function(x) median(x$a), numeric(1)),
        mse_drift = vapply(acc, function(x) mean(x$d), numeric(1)))
}

###############################################################################
# 第二部分：重抽樣
#
#  (A) 參數不確定性：以重抽樣得到 alpha / beta / kappa / drift 的標準誤與信賴區間
#      - "poisson" 參數式：D* ~ Poisson(E * m_hat)
#      - "pearson" 半參數式：重抽 Pearson 殘差，D* = m_hat*E + r* * sqrt(m_hat*E)
#        （能反映過度離散；若資料確實是 Poisson，兩者結果應接近）
#
#  (B) 預測不確定性：重抽 kappa 的一階差（創新項）
#      - "iid"   獨立重抽
#      - "block" 移動區塊重抽，保留時間相依性
###############################################################################

## ---- (A) 參數不確定性 -----------------------------------------------------
resample_params <- function(D, E, B = 500,
                            scheme = c("poisson", "pearson"),
                            smoother = NULL, seed = NULL) {
  scheme <- match.arg(scheme)
  if (!is.null(seed)) set.seed(seed)
  A <- nrow(D); Tn <- ncol(D)
  fit_once <- function(Dm) {
    Z <- log(pmax(Dm, 0.5) / E)
    if (!is.null(smoother)) Z <- smoother(Z)
    lc_svd_fit(Z)
  }
  f0 <- fit_once(D)
  mu <- E * exp(outer(f0$a, rep(1, Tn)) + outer(f0$b, f0$k))
  res_p <- as.vector((D - mu) / sqrt(pmax(mu, 1e-8)))   # Pearson 殘差
  
  A_ <- matrix(NA_real_, B, A); B_ <- matrix(NA_real_, B, A)
  K_ <- matrix(NA_real_, B, Tn); dr <- numeric(B)
  for (b in seq_len(B)) {
    Db <- if (scheme == "poisson")
      matrix(rpois(length(mu), mu), A, Tn)
    else
      pmax(round(mu + sample(res_p, length(mu), TRUE) * sqrt(pmax(mu, 1e-8))), 0)
    f <- fit_once(Db)
    A_[b, ] <- f$a; B_[b, ] <- f$b; K_[b, ] <- f$k
    dr[b] <- (f$k[Tn] - f$k[1]) / (Tn - 1)
  }
  ci <- function(M) t(apply(M, 2, quantile, c(.025, .5, .975), na.rm = TRUE))
  list(fit = f0, scheme = scheme,
       alpha = data.frame(age = rownames(D), est = f0$a,
                          se = apply(A_, 2, sd), ci(A_)),
       beta  = data.frame(age = rownames(D), est = f0$b,
                          se = apply(B_, 2, sd), ci(B_)),
       kappa = data.frame(year = colnames(D), est = f0$k,
                          se = apply(K_, 2, sd), ci(K_)),
       drift = c(est = (f0$k[Tn] - f0$k[1]) / (Tn - 1), se = sd(dr),
                 quantile(dr, c(.025, .5, .975))),
       draws = list(a = A_, b = B_, k = K_, drift = dr))
}

## ---- (B) 預測不確定性：重抽 kappa 的創新項 -------------------------------
resample_forecast <- function(fit, h = 26, B = 2000,
                              scheme = c("iid", "block"), block = 4,
                              seed = NULL) {
  scheme <- match.arg(scheme)
  if (!is.null(seed)) set.seed(seed)
  k <- fit$k[!is.na(fit$k)]; n <- length(k)
  dr <- (k[n] - k[1]) / (n - 1)
  eps <- diff(k) - dr                       # 去除 drift 後的創新項
  
  draw <- function() {
    if (scheme == "iid") sample(eps, h, TRUE)
    else {                                   # 移動區塊
      nb <- ceiling(h / block)
      st <- sample(seq_len(length(eps) - block + 1), nb, TRUE)
      as.vector(unlist(lapply(st, function(s) eps[s:(s + block - 1)])))[1:h]
    }
  }
  paths <- t(replicate(B, k[n] + cumsum(dr + draw())))
  e0 <- apply(paths, 2, function(col)
    vapply(col, function(x) life_table(exp(fit$a + fit$b * x))$e0, numeric(1)))
  list(kappa = paths, e0 = e0, drift = dr, scheme = scheme,
       kappa_q = apply(paths, 2, quantile, c(.1, .5, .9)),
       e0_q    = apply(e0,    2, quantile, c(.1, .5, .9)))
}

## ---- 重抽樣標準誤 vs 蒙地卡羅真實標準差 ----------------------------------
# 回答「在真實資料上算出來的 SE 可不可信」
resample_calibration <- function(truth, N, B = 60, reps = 40,
                                 smoother = NULL, seed = 77) {
  set.seed(seed)
  E <- N * truth$w
  gen <- function() matrix(rpois(length(E), E * truth$m), nrow(E), ncol(E),
                           dimnames = dimnames(truth$m))
  fit_once <- function(Dm) {
    Z <- log(pmax(Dm, 0.5) / E); if (!is.null(smoother)) Z <- smoother(Z)
    lc_svd_fit(Z)
  }
  MC <- t(replicate(reps, fit_once(gen())$b))          # 真實抽樣分布
  bs <- resample_params(gen(), E, B = B, smoother = smoother)$beta$se
  data.frame(age = truth$ages, sd_true = apply(MC, 2, sd), se_boot = bs,
             ratio = bs / apply(MC, 2, sd))
}

###############################################################################
# 執行
###############################################################################
if (sys.nframe() == 0) {
  
  truth <- mc_truth("Female", 2001, 2024)
  
  ## ---- 三方向 lambda（耗時；先用小 reps）----
  cat("三方向 lambda 格點搜尋，N = 2e5\n")
  g <- lambda_grid3(truth, 2e5, reps = 8)
  g <- g[order(g$sse_b), ]
  cat("\n依 SSE(beta) 排序，前 8 名：\n")
  print(head(g, 8), digits = 4, row.names = FALSE)
  cat("\n最差 3 名：\n"); print(tail(g, 3), digits = 4, row.names = FALSE)
  cat("\n依 SSE(alpha) 排序，前 3 名：\n")
  print(head(g[order(g$sse_a), ], 3), digits = 4, row.names = FALSE)
  write.csv(g, "output/tables/lambda_grid3.csv", row.names = FALSE)
  
  ## ---- 重抽樣：在真實資料上 ----
  dat  <- load_data(sex = "Female")
  keep <- which(dat$years >= 2001 & dat$years <= 2024)
  D <- dat$D[, keep]; E <- dat$E[, keep]
  
  sm <- function(Z) smooth3(Z, lxx = 1, lxt = 2, ltt = 1)   # 上面選出的組合
  
  cat("\n參數不確定性（無修勻 vs 修勻），B = 300\n")
  rp_raw <- resample_params(D, E, B = 300, scheme = "poisson", seed = 1)
  rp_sm  <- resample_params(D, E, B = 300, scheme = "poisson",
                            smoother = sm, seed = 1)
  cmp <- data.frame(age = dat$ages,
                    beta_raw = rp_raw$beta$est, se_raw = rp_raw$beta$se,
                    beta_sm  = rp_sm$beta$est,  se_sm  = rp_sm$beta$se,
                    se_ratio = rp_sm$beta$se / rp_raw$beta$se)
  print(cmp, digits = 3, row.names = FALSE)
  cat(sprintf("\ndrift：無修勻 %.4f (SE %.4f) / 修勻 %.4f (SE %.4f)\n",
              rp_raw$drift["est"], rp_raw$drift["se"],
              rp_sm$drift["est"],  rp_sm$drift["se"]))
  
  ## Pearson 半參數式，檢查是否有過度離散
  rp_pe <- resample_params(D, E, B = 300, scheme = "pearson", seed = 1)
  cat(sprintf("Pearson 方案的 drift SE = %.4f（若遠大於 Poisson 的 %.4f，表示過度離散）\n",
              rp_pe$drift["se"], rp_raw$drift["se"]))
  
  ## ---- 預測不確定性 ----
  cat("\n預測區間（kappa 創新項重抽），到 2050\n")
  fc_i <- resample_forecast(rp_sm$fit, h = 26, scheme = "iid",   seed = 2)
  fc_b <- resample_forecast(rp_sm$fit, h = 26, scheme = "block", seed = 2)
  cat(sprintf("  iid   : e0(2050) 80%%PI = [%.2f, %.2f]\n",
              fc_i$e0_q[1, 26], fc_i$e0_q[3, 26]))
  cat(sprintf("  block : e0(2050) 80%%PI = [%.2f, %.2f]\n",
              fc_b$e0_q[1, 26], fc_b$e0_q[3, 26]))
  
  ## ---- 重抽樣 SE 的校準 ----
  cat("\n重抽樣 SE vs 蒙地卡羅真實 SD（N = 2e5，修勻版）\n")
  cal <- resample_calibration(truth, 2e5, B = 40, reps = 30, smoother = sm)
  print(cal, digits = 3, row.names = FALSE)
  cat(sprintf("比值中位數 = %.3f（< 1 表示重抽樣低估了不確定性）\n",
              median(cal$ratio, na.rm = TRUE)))
  
  ## ---- 圖 ----
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  A <- nrow(D); mid <- make_grid(A)$mid
  plot(mid, rp_raw$beta$est, type = "l", lwd = 2, xlab = "Age",
       ylab = expression(beta[x]), main = "beta 與 95% 重抽樣區間",
       ylim = safe_range(rp_raw$beta[, 4:6], rp_sm$beta[, 4:6]))
  polygon(c(mid, rev(mid)), c(rp_raw$beta[,4], rev(rp_raw$beta[,6])),
          col = rgb(0,0,0,.12), border = NA)
  lines(mid, rp_sm$beta$est, col = 2, lwd = 2)
  polygon(c(mid, rev(mid)), c(rp_sm$beta[,4], rev(rp_sm$beta[,6])),
          col = rgb(1,0,0,.15), border = NA)
  legend("topright", c("無修勻","修勻"), col = 1:2, lwd = 2, bty = "n", cex = .8)
  
  plot(mid, cmp$se_ratio, type = "b", pch = 16, xlab = "Age",
       ylab = "SE(修勻) / SE(無修勻)", main = "修勻對 SE 的效果")
  abline(h = 1, lty = 2)
  
  yf <- 2024 + 1:26
  plot(yf, fc_i$e0_q[2, ], type = "n", ylim = safe_range(fc_i$e0_q, fc_b$e0_q),
       xlab = "Year", ylab = expression(e[0]), main = "e0 預測區間")
  polygon(c(yf, rev(yf)), c(fc_i$e0_q[1,], rev(fc_i$e0_q[3,])),
          col = rgb(0,0,1,.15), border = NA)
  polygon(c(yf, rev(yf)), c(fc_b$e0_q[1,], rev(fc_b$e0_q[3,])),
          col = rgb(1,0,0,.12), border = NA)
  lines(yf, fc_i$e0_q[2,], lwd = 2, col = 4)
  legend("topleft", c("iid","block"), fill = c(rgb(0,0,1,.3), rgb(1,0,0,.3)),
         bty = "n", cex = .8)
  
  hist(rp_sm$draws$drift, breaks = 30, col = "grey80", border = "white",
       xlab = "drift", main = "drift 的重抽樣分布")
  abline(v = rp_sm$drift["est"], col = 2, lwd = 2)
  par(op)
}

###############################################################################
# 三方向格點的結果（先跑過，reps = 8）
#
#  N = 2e5，依 SSE(beta) 排序
#    lxx  lxt  ltt   SSE(beta)   SSE(alpha)
#     20    2    1     0.00654      5.131
#      1    2   20     0.00694      0.219
#      1    2    1     0.00752      0.203     <- 兼顧兩者的選擇
#      0    0    0     0.04831      0.296     （不修勻）
#     20    0    0     0.12353      5.976     （比不修勻更糟）
#
#  N = 1e6 的排序型態相同。
#
#  三個結論：
#   1. lxt 是主力。兩個暴露數的前八名全部 lxt > 0，最差三名全部 lxt = 0。
#      純年齡方向的修勻 (20,0,0) 比完全不修勻還糟。
#   2. lxx 大會毀掉 alpha（SSE 從 0.2 跳到 5-7）。五齡組資料的年齡方向已被
#      聚合平滑過，再罰就是純粹的過度修勻。
#   3. 先前「alpha 與 beta 要相反的 lambda」是三方向耦合造成的假象。
#      拆開後 (1,2,1) 的 SSE(beta) 只比最佳差一點，SSE(alpha) 好二十倍。
#
#  為什麼 lxt 最有效：在 LC 底下
#      d^2 log m / dx dt = beta'_x * kappa'_t
#  所以交叉項本質上就是對 beta_x 粗糙度的懲罰，是三個方向裡唯一直接
#  瞄準 beta 的。這在數學上是「beta 的差分懲罰」，與老師的兩項指示
#  （懲罰加在 beta 上、放棄差分懲罰）落在交界處，報告時宜如實呈現證據，
#  由老師決定方向。
###############################################################################