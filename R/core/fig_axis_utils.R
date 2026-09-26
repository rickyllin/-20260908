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

## ===================== 中文圖檔輸出 ========================================
#' 開啟可正確顯示中文的 png 裝置
#'
#' macOS 上 R 的 quartz 裝置未必可用（capabilities("quartz") 為 FALSE），
#' 預設的 png 裝置找不到 CJK 字型時會把中文畫成方框（tofu）。
#' 本函式改用 cairo 後端並指定字型家族，找不到時退回預設。
#'
#' @param prefer 字型家族的優先序；第一個能由 fontconfig 解析者勝出
png_cjk <- function(file, width = 1900, height = 700, res = 150,
                    prefer = c("Heiti TC", "PingFang TC", "Hiragino Sans GB",
                               "Arial Unicode MS")) {
  fam <- NULL
  if (capabilities("cairo")) {
    avail <- tryCatch(system("fc-list :lang=zh-tw family", intern = TRUE),
                      error = function(e) character(0), warning = function(w) character(0))
    avail <- unique(trimws(unlist(strsplit(avail, ","))))
    hit <- prefer[prefer %in% avail]
    fam <- if (length(hit)) hit[1] else NULL
  }
  if (!is.null(fam)) {
    png(file, width = width, height = height, res = res, type = "cairo", family = fam)
  } else {
    warning("找不到中文字型，圖中的中文可能顯示為方框")
    png(file, width = width, height = height, res = res)
  }
  invisible(fam)
}
