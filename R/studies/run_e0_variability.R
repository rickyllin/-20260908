###############################################################################
# 平均餘命的偏誤與離散（回應 0929 會議：圖 5 除偏誤外，各方法的變異散佈
# 差異似乎也很大，可補 summary 表；並想看平均餘命是否與死亡率一樣，
# 在小人口下有劇烈震盪）
#
#   本檔產生兩項結果。
#   其一是各估計量零歲平均餘命的分布摘要：中位偏誤、平均偏誤、標準差、
#   四分位距、全距，以及落在真值正負半歲內的重複比例。偏誤與離散是兩件
#   不同的事，前者由估計方程的中心決定、後者由其變異數決定，
#   本文的校正只動前者，故這張表正是檢驗「校正是否以變異數為代價」的地方。
#
#   其二是平均餘命的逐年序列。死亡率的震盪眾所周知，但平均餘命是對
#   整條死亡率曲線的加權平均，震盪未必同樣劇烈。本檔以單次實現畫出
#   各年的平均餘命，並以逐年變動的平均絕對值量化震盪幅度。
#
# 輸出：output/tables/tableM3_e0_variability.csv
#       output/tables/tableM3_e0_oscillation.csv
#       output/figures/figT8_e0_oscillation.png
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/fig_axis_utils.R")

SEED <- 20260926; REPS <- 100; NS <- c(1e4, 5e4, 2e5); YEARS <- 2001:2024
dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0); yrs <- as.numeric(colnames(D0))
truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage <- E0[, Tn] / sum(E0[, Tn]); BF <- make_logbias(0.5)
e0_true_t <- apply(mtrue, 2, function(m) life_table(m)$e0)
e0_true   <- e0_true_t[Tn]

fit_cellw <- function(D, E) {
  g <- lc_analytic(D, E, bfun = BF, shrink = 1)
  mu <- E * exp(outer(g$a, rep(1, ncol(D))) + outer(g$b, g$k))
  lmc <- log(pmax(D, 0.5) / E) - BF(mu); W <- mu / mean(mu)
  am <- rowSums(W * lmc) / rowSums(W); Z <- (lmc - am) * sqrt(W)
  sv <- svd(Z); u1 <- sv$u[, 1]; v1 <- sv$v[, 1]
  if (sum(u1) < 0) { u1 <- -u1; v1 <- -v1 }
  bb <- u1 / sum(u1); kk <- sv$d[1] * v1 * sum(u1) / sqrt(pmax(colMeans(W), 1e-12))
  kk <- kk - mean(kk); sb <- sum(bb); list(a = am, b = bb / sb, k = kk * sb)
}
FIT <- list("標準 LC"            = function(D, E) lc_svd_fit(log(pmax(D, 0.5) / E)),
            "Poisson MLE"        = function(D, E) lc_poisson_firth(D, E, firth = FALSE),
            "Firth"              = function(D, E) lc_poisson_firth(D, E, firth = TRUE),
            "中心化"             = function(D, E) lc_analytic(D, E, bfun = BF, shrink = 1),
            "中心化＋加權"       = fit_cellw,
            "僅校正 alpha＋加權" = function(D, E) lc_alpha_only(D, E, bfun = BF, weight = TRUE))
K <- length(FIT)

## ---- 一、分布摘要，並同時記錄平均餘命的逐年序列 --------------------------
SEL <- c("標準 LC", "Firth", "中心化＋加權")
res <- list(); osc <- list(); SER1 <- vector("list", length(NS))
for (ni in seq_along(NS)) {
  N <- NS[ni]; E <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu <- E * mtrue
  eS <- matrix(NA_real_, REPS, K)
  SERIES <- array(NA_real_, c(REPS, Tn, length(SEL) + 1),
                  dimnames = list(NULL, NULL, c("觀測（零格以 0.5 替代）", SEL)))
  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000 * ni + r)              # 與第伍節各表同一組種子
    D <- matrix(rpois(length(E), mu), A, dimnames = dimnames(E))
    SERIES[r, , 1] <- apply(pmax(D, 0.5) / E, 2, function(m) life_table(m)$e0)
    for (j in seq_len(K)) {
      f <- tryCatch(FIT[[j]](D, E), error = function(e) NULL)
      if (is.null(f) || !all(is.finite(f$a)) || !all(is.finite(f$b))) next
      mh <- exp(outer(f$a, rep(1, Tn)) + outer(f$b, f$k))
      eS[r, j] <- tryCatch(life_table(mh[, Tn])$e0, error = function(e) NA_real_)
      nm <- names(FIT)[j]
      if (nm %in% SEL)
        SERIES[r, , nm] <- tryCatch(apply(mh, 2, function(m) life_table(m)$e0),
                                    error = function(e) rep(NA_real_, Tn))
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }
  bias <- eS - e0_true
  res[[ni]] <- data.frame(
    N = N, 估計量 = names(FIT),
    中位偏誤 = round(apply(bias, 2, median, na.rm = TRUE), 3),
    平均偏誤 = round(colMeans(bias, na.rm = TRUE), 3),
    標準差   = round(apply(bias, 2, sd, na.rm = TRUE), 3),
    四分位距 = round(apply(bias, 2, function(v) IQR(v, na.rm = TRUE)), 3),
    全距     = round(apply(bias, 2, function(v) diff(range(v, na.rm = TRUE))), 3),
    命中率   = round(apply(bias, 2, function(v) mean(abs(v) <= 0.5, na.rm = TRUE)), 3),
    均方根誤差 = round(apply(bias, 2, function(v) sqrt(mean(v^2, na.rm = TRUE))), 3),
    stringsAsFactors = FALSE)
  cat(sprintf("\n===== N = %.0e =====\n", N)); print(res[[ni]][, -1], row.names = FALSE)

  ## 震盪幅度：每次重複算一條序列的逐年變動，再取 100 次的中位數
  flux <- function(M, f) median(apply(M, 1, function(v)
            if (all(is.finite(v))) f(abs(diff(v))) else NA_real_), na.rm = TRUE)
  rng  <- function(M) median(apply(M, 1, function(v)
            if (all(is.finite(v))) diff(range(v)) else NA_real_), na.rm = TRUE)
  nmz <- dimnames(SERIES)[[3]]
  osc[[ni]] <- data.frame(N = N, 序列 = c("真值", nmz),
    逐年變動均值 = round(c(mean(abs(diff(e0_true_t))),
                     sapply(nmz, function(z) flux(SERIES[, , z], mean))), 3),
    逐年變動最大 = round(c(max(abs(diff(e0_true_t))),
                     sapply(nmz, function(z) flux(SERIES[, , z], max))), 3),
    粗糙度 = round(c(mean(abs(diff(e0_true_t, differences = 2))),
                 sapply(nmz, function(z) median(apply(SERIES[, , z], 1, function(v)
                   if (all(is.finite(v))) mean(abs(diff(v, differences = 2))) else NA_real_),
                   na.rm = TRUE))), 3),
    全距 = round(c(diff(range(e0_true_t)),
                   sapply(nmz, function(z) rng(SERIES[, , z]))), 3),
    stringsAsFactors = FALSE)
  SER1[[ni]] <- c(list("真值" = e0_true_t),
                  lapply(nmz, function(z) SERIES[1, , z]))
  names(SER1[[ni]]) <- c("真值", nmz)
}
tab <- do.call(rbind, res)
write.csv(tab, "output/tables/tableM3_e0_variability.csv", row.names = FALSE)
oscT <- do.call(rbind, osc)
write.csv(oscT, "output/tables/tableM3_e0_oscillation.csv", row.names = FALSE)
print(oscT, row.names = FALSE)

## ---- 二、單次實現的逐年序列圖 ----------------------------------------------
COL <- c("真值" = "black", "觀測（零格以 0.5 替代）" = "#B03A2E",
         "標準 LC" = "#D98880", "Firth" = "#7D3C98", "中心化＋加權" = "#1E8449")
LTY <- c(1, 1, 2, 4, 1); LWD <- c(3.0, 1.6, 2.2, 2.2, 2.6)
png_cjk("output/figures/figT8_e0_oscillation.png", width = 2250, height = 860, res = 160)
par(mfrow = c(1, 3), mar = c(4.6, 4.6, 3.2, 0.8), mgp = c(2.8, 0.8, 0))
for (ni in seq_along(NS)) {
  ser <- SER1[[ni]]; yl <- range(unlist(ser), na.rm = TRUE)
  yl <- yl + c(-1, 1) * 0.06 * diff(yl)
  plot(yrs, ser[[1]], type = "n", ylim = yl, xlab = "年份",
       ylab = "零歲平均餘命（歲）",
       main = sprintf("N = %s", formatC(NS[ni], format = "d", big.mark = ",")))
  grid(col = "grey88", lty = 1)
  for (j in seq_along(ser))
    lines(yrs, ser[[j]], col = COL[names(ser)[j]], lwd = LWD[j], lty = LTY[j],
          type = if (j == 2) "b" else "l", pch = 16, cex = 0.6)
  if (ni == 1) legend("bottomright", bty = "n", cex = 1.0, legend = names(ser),
                      col = COL[names(ser)], lwd = LWD, lty = LTY)
}
dev.off()
cat("\n已輸出 tableM3 與 figT8\n")
