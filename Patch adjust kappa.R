###############################################################################
# 修正 adjust_kappa：e0dagger 對 kappa 是單峰而非單調，
# 寬區間 + extendInt 會兩端同號而失敗。改成網格掃描 + 取最近的根。
#
# 用法：source("rabbi_mazzuco_replication.R") 之後再 source 本檔，
#       或直接把下面的函式覆蓋原檔中的同名函式。
###############################################################################

adjust_kappa <- function(a, b, k0, target_fun, target_obs,
                         span = 30, ngrid = 400, verbose = TRUE) {
  n <- length(k0)
  out <- rep(NA_real_, n)
  nroots <- integer(n)
  
  for (j in seq_len(n)) {
    g <- seq(k0[j] - span, k0[j] + span, length.out = ngrid)
    v <- vapply(g, function(x) target_fun(exp(a + b * x)), numeric(1)) - target_obs[j]
    ok <- is.finite(v)
    if (sum(ok) < 2) next
    g <- g[ok]; v <- v[ok]
    
    idx <- which(sign(v[-length(v)]) * sign(v[-1]) < 0)
    nroots[j] <- length(idx)
    if (!length(idx)) next
    
    roots <- vapply(idx, function(i)
      uniroot(function(x) target_fun(exp(a + b * x)) - target_obs[j],
              c(g[i], g[i + 1]), tol = 1e-10)$root, numeric(1))
    out[j] <- roots[which.min(abs(roots - k0[j]))]   # 留在與資料同一側的分支
  }
  
  if (verbose) {
    if (any(is.na(out)))
      warning(sprintf("adjust_kappa: %d 年在 span=%g 內找不到根（索引 %s）",
                      sum(is.na(out)), span,
                      paste(which(is.na(out)), collapse = ",")))
    if (any(nroots > 1))
      message(sprintf("adjust_kappa: %d 年有多重根，已取最接近原始 kappa 者",
                      sum(nroots > 1)))
  }
  out
}

## ---- 診斷工具：畫出目標量對 kappa 的曲線，確認單調性 ----------------------
plot_target_curve <- function(a, b, target_fun, klim = c(-80, 80),
                              k_obs = NULL, main = "target vs kappa") {
  g <- seq(klim[1], klim[2], length.out = 400)
  v <- vapply(g, function(x) target_fun(exp(a + b * x)), numeric(1))
  plot(g, v, type = "l", lwd = 2, xlab = expression(kappa), ylab = "target",
       main = main)
  if (!is.null(k_obs)) {
    rug(k_obs, col = 2)
    abline(v = range(k_obs), col = 2, lty = 2)
  }
  i <- which.max(v)
  points(g[i], v[i], pch = 19, col = 4)
  text(g[i], v[i], sprintf(" 峰值 %.2f @ k=%.1f", v[i], g[i]), pos = 4, cex = .8)
  invisible(data.frame(k = g, value = v))
}

## ---- 繪圖修正：ylim 濾掉非有限值，標題改 ASCII 避免中文亂碼 ---------------
safe_range <- function(...) range(c(...), na.rm = TRUE, finite = TRUE)