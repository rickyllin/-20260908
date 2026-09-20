###############################################################################
# 繪圖用的座標軸工具：把 R 預設的電腦式科學記號（1e-03、5e+04）
# 改為學術書寫的 10^{-3}、5 x 10^{4}
###############################################################################

#' 依數值產生 plotmath 標籤：1e4 -> 10^4；2e5 -> 2 x 10^5
pow_lab <- function(v, digits = 1) {
  lapply(v, function(x) {
    if (!is.finite(x) || x == 0) return(bquote(0))
    e <- floor(log10(abs(x))); m <- x / 10^e
    if (abs(m - 1) < 1e-8) bquote(10^.(e))
    else bquote(.(signif(m, digits)) %*% 10^.(e))
  })
}

#' 在 side 邊畫出以 10 的次方標示的座標軸
axis_pow <- function(side, at, digits = 1, ...) {
  axis(side, at = at, labels = as.expression(pow_lab(at, digits)), ...)
}

#' 自動選取涵蓋 range(v) 的 10 的整數次方刻度
pow_ticks <- function(v, n = 5) {
  v <- v[is.finite(v) & v > 0]
  lo <- floor(log10(min(v))); hi <- ceiling(log10(max(v)))
  e <- lo:hi
  if (length(e) > n) e <- e[seq(1, length(e), length.out = n)]
  10^round(e)
}
