###############################################################################
# 蒙地卡羅實驗：暴露數 x 估計方法
#
#   真值   ：以 Female 2001-2024 的 LC 配適為 (alpha, beta, kappa)
#   生成   ：E_xt = N * w_xt（w 為各年的年齡結構），
#            D_xt ~ Poisson(E_xt * exp(alpha_x + beta_x * kappa_t))
#   比較   ：純 LC（SVD）vs LASSO 修勻後的 LC，多組 lambda
#   評估   ：alpha/beta/kappa 的偏誤、變異數、MSE；drift；崩潰率
#
# 需先 source("rabbi_mazzuco_replication.R") 與 source("patch_adjust_kappa.R")
###############################################################################

source("rabbi_mazzuco_replication.R")
source("Patch adjust kappa.R")

## ===================== 1. 真值 =============================================
mc_truth <- function(sex = "Female", y0 = 2001, y1 = 2024) {
  dat  <- load_data(sex = sex)
  keep <- which(dat$years >= y0 & dat$years <= y1)
  D <- dat$D[, keep, drop = FALSE]; E <- dat$E[, keep, drop = FALSE]
  f <- lc_svd_fit(log(D / E))
  list(a = f$a, b = f$b, k = f$k,
       m = exp(outer(f$a, rep(1, ncol(D))) + outer(f$b, f$k)),
       w = sweep(E, 2, colSums(E), "/"),      # 各年的年齡結構（和為 1）
       ages = dat$ages, years = dat$years[keep],
       varexp = f$varexp)
}

## ===================== 2. 單次模擬與估計 ===================================
# lambda 的三個方向以 (l, l/5, l) 連動：年齡、交叉、時間
lam_vec <- function(l) c(lxx = l, lxt = l / 5, ltt = l)

mc_one <- function(truth, N, lambdas, zero_sub = 0.5) {
  E <- N * truth$w
  D <- matrix(rpois(length(E), E * truth$m), nrow(E), ncol(E),
              dimnames = dimnames(truth$m))
  logm <- log(pmax(D, zero_sub) / E)
  
  est <- list()
  est[["LC"]] <- lc_svd_fit(logm)
  for (l in lambdas) {
    lv <- lam_vec(l)
    S  <- smooth_l1_2d(logm, lv["lxx"], lv["lxt"], lv["ltt"], method = "rq")
    est[[sprintf("LASSO%g", l)]] <- lc_svd_fit(S)
  }
  list(est = est, nzero = sum(D == 0))
}

## ===================== 3. 主迴圈 ===========================================
cor_centered <- function(v, t) cor(v - mean(v), t - mean(t))

mc_run <- function(truth,
                   Ns = c(1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6),
                   lambdas = c(2, 5, 10, 20),
                   reps = 50, seed = 2026, verbose = TRUE) {
  set.seed(seed)
  A <- length(truth$a); Tn <- length(truth$k)
  drift_true <- (truth$k[Tn] - truth$k[1]) / (Tn - 1)
  methods <- c("LC", sprintf("LASSO%g", lambdas))
  out <- list()
  
  for (N in Ns) {
    store <- setNames(lapply(methods, function(m)
      list(a = matrix(NA_real_, reps, A), b = matrix(NA_real_, reps, A),
           k = matrix(NA_real_, reps, Tn))), methods)
    nzero <- numeric(reps)
    
    for (r in seq_len(reps)) {
      one <- mc_one(truth, N, lambdas)
      nzero[r] <- one$nzero
      for (m in methods) {
        store[[m]]$a[r, ] <- one$est[[m]]$a
        store[[m]]$b[r, ] <- one$est[[m]]$b
        store[[m]]$k[r, ] <- one$est[[m]]$k
      }
    }
    
    for (m in methods) {
      S <- store[[m]]
      dr  <- (S$k[, Tn] - S$k[, 1]) / (Tn - 1)
      cb  <- apply(S$b, 1, cor_centered, t = truth$b)
      row <- data.frame(
        N = N, method = m, zero_cells = mean(nzero),
        ## alpha
        a_bias = mean(colMeans(S$a) - truth$a),
        a_bias2 = sum((colMeans(S$a) - truth$a)^2),
        a_var  = sum(apply(S$a, 2, var)),
        a_mse_med = median(rowSums(sweep(S$a, 2, truth$a)^2)),
        ## beta
        b_bias2 = sum((colMeans(S$b) - truth$b)^2),
        b_var  = sum(apply(S$b, 2, var)),
        b_mse_med = median(rowSums(sweep(S$b, 2, truth$b)^2)),
        b_cor_med = median(cb),
        breakdown = mean(cb < 0.8),          # 崩潰率
        ## kappa / drift
        k_bias2 = sum((colMeans(S$k) - truth$k)^2),
        k_var  = sum(apply(S$k, 2, var)),
        drift_med = median(dr),
        drift_bias = median(dr) - drift_true,
        drift_iqr = IQR(dr))
      out[[length(out) + 1]] <- row
    }
    if (verbose) cat(sprintf("N = %9.0f 完成（平均零格 %.1f / %d）\n",
                             N, mean(nzero), A * Tn))
  }
  res <- do.call(rbind, out)
  attr(res, "drift_true") <- drift_true
  res
}

## ===================== 4. 變異數估計的準確度（選用）=======================
# 用參數式 bootstrap 在「單一模擬資料集」上估 beta 的標準誤，
# 再跟蒙地卡羅的真實標準差比較，回答「變異數估得準不準」
mc_se_check <- function(truth, N, lambda = 5, B = 100, reps = 50, seed = 99) {
  set.seed(seed)
  A <- length(truth$a)
  ## (i) 蒙地卡羅真實標準差
  MC <- t(replicate(reps, mc_one(truth, N, lambda)$est[[sprintf("LASSO%g", lambda)]]$b))
  sd_mc <- apply(MC, 2, sd)
  ## (ii) 單一資料集上的參數式 bootstrap
  E <- N * truth$w
  D <- matrix(rpois(length(E), E * truth$m), nrow(E), ncol(E))
  lv <- lam_vec(lambda)
  f0 <- lc_svd_fit(smooth_l1_2d(log(pmax(D, .5) / E), lv[1], lv[2], lv[3]))
  m0 <- exp(outer(f0$a, rep(1, ncol(E))) + outer(f0$b, f0$k))
  BS <- t(replicate(B, {
    Db <- matrix(rpois(length(E), E * m0), nrow(E), ncol(E))
    lc_svd_fit(smooth_l1_2d(log(pmax(Db, .5) / E), lv[1], lv[2], lv[3]))$b
  }))
  sd_bs <- apply(BS, 2, sd)
  data.frame(age = truth$ages, sd_mc = sd_mc, sd_boot = sd_bs,
             ratio = sd_bs / sd_mc)
}

## ===================== 5. 繪圖 =============================================
mc_plots <- function(res, file = NULL) {
  ms <- unique(res$method); cols <- setNames(seq_along(ms), ms)
  pick <- function(m, col) res[[col]][res$method == m]
  Ns <- sort(unique(res$N))
  if (!is.null(file)) png(file, 1200, 950, res = 120)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  plot(Ns, pick(ms[1], "b_mse_med"), type = "n", log = "xy",
       ylim = safe_range(res$b_mse_med), xlab = "Exposure N",
       ylab = "median SSE(beta)", main = "beta 的 MSE（中位數）")
  for (m in ms) lines(Ns, pick(m, "b_mse_med"), type = "b", pch = 16, col = cols[m])
  legend("bottomleft", ms, col = cols, lwd = 1, pch = 16, bty = "n", cex = .7)
  
  plot(Ns, pick(ms[1], "a_bias"), type = "n", log = "x",
       ylim = safe_range(res$a_bias), xlab = "Exposure N",
       ylab = "mean bias(alpha)", main = "alpha 的平均偏誤")
  for (m in ms) lines(Ns, pick(m, "a_bias"), type = "b", pch = 16, col = cols[m])
  abline(h = 0, lty = 2)
  
  dt <- attr(res, "drift_true")
  plot(Ns, pick(ms[1], "drift_med"), type = "n", log = "x",
       ylim = safe_range(res$drift_med, dt), xlab = "Exposure N",
       ylab = "median drift", main = sprintf("drift（真值 %.4f）", dt))
  for (m in ms) lines(Ns, pick(m, "drift_med"), type = "b", pch = 16, col = cols[m])
  abline(h = dt, lty = 2, col = 2)
  
  plot(Ns, pick(ms[1], "breakdown"), type = "n", log = "x", ylim = c(0, 1),
       xlab = "Exposure N", ylab = "P(cor(beta) < 0.8)", main = "崩潰率")
  for (m in ms) lines(Ns, pick(m, "breakdown"), type = "b", pch = 16, col = cols[m])
  par(op); if (!is.null(file)) dev.off()
}

## ===================== 6. 執行 =============================================
if (sys.nframe() == 0) {
  truth <- mc_truth("Female", 2001, 2024)
  cat(sprintf("真值：%d 齡組 x %d 年，第一成分解釋 %.3f，drift = %.4f\n",
              length(truth$a), length(truth$k), truth$varexp,
              (truth$k[length(truth$k)] - truth$k[1]) / (length(truth$k) - 1)))
  
  ## 先小規模試跑確認流程，再放大 reps
  res <- mc_run(truth, reps = 50)
  print(res[, c("N","method","zero_cells","a_bias","b_mse_med","b_cor_med",
                "breakdown","drift_med","drift_iqr")],
        digits = 4, row.names = FALSE)
  write.csv(res, "mc_results.csv", row.names = FALSE)
  mc_plots(res, file = "mc_plots.png")
  
  ## 變異數估計的準確度（挑一個中等暴露數）
  cat("\n參數式 bootstrap 的 SE 與蒙地卡羅 SD 比較（N = 2e5, lambda = 5）：\n")
  print(mc_se_check(truth, 2e5, lambda = 5, B = 60, reps = 40), digits = 3,
        row.names = FALSE)
}

###############################################################################
# 我先跑過的結果（Female 2001-2024 為真值，reps = 120）
#
#  N        零格/528   alpha 偏誤   beta 相關中位數(LC / LASSO5)   drift 中位數
#  1e4       213.8      +0.442        -0.245  /  -0.429            -0.12 / -0.11
#  5e4        89.6      +0.031        -0.072  /  +0.036            -0.07 / -0.25
#  2e5        21.5      -0.050        +0.531  /  +0.602            -0.22 / -0.42
#  1e6         0.4      -0.021        +0.851  /  +0.866            -0.35 / -0.41
#                                                        （drift 真值 -0.3678）
#
# 四個要點：
#  1. N <= 2e4 時零格佔四成，主導一切的是零格處理而非估計方法。alpha 的
#     正向偏誤 (+0.44) 來自 0.5 取代，lambda 從 0 加到 20 都救不了。
#  2. N <= 2e5 時兩種方法的 beta 都已崩潰（相關中位數 < 0.62，崩潰率 > 0.8）。
#     管線本身沒問題：N = 1e8 時相關為 0.999。
#  3. drift 被系統性壓平——小人口會低估死亡改善速度。這是有方向的偏誤，
#     比變異數變大嚴重。LASSO 修勻在 N >= 5e4 時把 drift 拉回真值附近，
#     且四分位距減半，這是它最明確的貢獻。
#  4. alpha 與 beta 的最適 lambda 方向相反：N = 1e6 時 beta 的 MSE 隨 lambda
#     遞減，alpha 的 MSE 卻從 0.064 暴增到 5.15（過度修勻）。單一 lambda
#     必然是妥協。
#
# 另一個負面結果（已試過，不要再走）：把零格排除在保真項之外的加權版本
# 更差（N = 2e5 時 beta 相關掉到 0.11），因為丟掉零格等於丟掉「此處死亡率
# 很低」的訊息。
###############################################################################