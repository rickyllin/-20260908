###############################################################################
# 四張 EDA 診斷圖：從資料分布決定該用哪個估計量
#
#   ①②③ 為純粹的 EDA，拿到資料即可計算
#   ④   需要一次便宜的預備配適（取 mu_hat 以估噪音尺度）
#
#   用法：
#     dg <- lc_eda(D, E)        # D, E 為 年齡 x 年份 的矩陣
#     print(dg$summary)         # 數值摘要與建議
#     plot_lc_eda(dg)           # 四格圖
###############################################################################

## ===================== 診斷計算 ============================================
#' @param D,E  年齡 x 年份 的死亡數與暴露數矩陣
#' @param fitter 預備配適函式，需回傳 list(a,b,k)。建議用 Firth 版
#'               （見 fh_firth_kalman.R 的 lc_poisson_firth），零格多時較穩。
lc_eda <- function(D, E, fitter = NULL, age_mid = NULL) {
  A <- nrow(D); Tn <- ncol(D); cc <- A / Tn
  if (is.null(age_mid)) {
    n <- c(1, 4, rep(5, A - 2))
    age_mid <- c(0, cumsum(n)[-A]) + c(0.5, 2.0, rep(2.5, A - 2))
  }
  logm <- log(pmax(D, 0.5) / E)

  ## ① 各年齡的死亡數量級
  d_med <- apply(D, 1, median)

  ## ② 零格落在哪些年齡
  zero_by_age <- rowMeans(D == 0)

  ## ③ 由資料直接估的噪音變異：Var(二階差)/6
  ##    二階差消去線性趨勢；係數 6 來自 Var(Delta^2 X) = 6 sigma^2
  d2 <- t(apply(logm, 1, diff, differences = 2))
  noise_var <- apply(d2, 1, var)  / 6

  ## ④ 奇異值譜與噪音上緣
  Z  <- logm - rowMeans(logm)
  sv <- svd(Z)$d
  if (is.null(fitter)) {
    ## 無預備配適：改由觀測 D 估 sigma（偏保守，門檻請放寬到 1.2）
    sigma <- sqrt(mean(1 / pmax(D, 0.5)))
    sigma_src <- "觀測 D（保守）"
  } else {
    f  <- fitter(D, E)
    mu <- E * exp(outer(f$a, rep(1, Tn)) + outer(f$b, f$k))
    sigma <- sqrt(mean(1 / pmax(mu, 1e-9)))
    sigma_src <- "預備配適的 mu_hat（建議）"
  }
  bulk  <- sigma * sqrt(Tn) * (1 + sqrt(cc))
  ratio <- sv[1] / bulk

  ## ---- 摘要與建議 ----
  verdict <- if (ratio > 2)        "訊號充分：beta 用 SVD 亦可"
        else if (ratio > 1.2)      "可還原：卜瓦松或加權 SVD(w=mu_hat)，加 FH 收縮"
        else if (ratio > 0.6)      "邊緣：方法選擇效益最大，用卜瓦松/Firth + FH"
        else                       "不可還原：不要估 beta，只報 alpha（Firth）與區間"
  need_firth <- min(d_med) < 5
  het_span   <- diff(range(log10(pmax(noise_var, 1e-12))))

  summary <- data.frame(
    項目 = c("最小單年死亡中位數", "死亡數跨年齡落差（倍）", "零死亡格比例（%）",
             "零格最集中的年齡", "噪音變異跨度（數量級）",
             "s1", "噪音上緣", "s1 / 上緣", "sigma 來源", "判讀"),
    數值 = c(sprintf("%.0f", min(d_med)),
             sprintf("%.0f", max(d_med) / pmax(min(d_med), 1)),
             sprintf("%.1f", 100 * mean(D == 0)),
             if (any(zero_by_age > 0)) rownames(D)[which.max(zero_by_age)] else "無",
             sprintf("%.1f", het_span),
             sprintf("%.2f", sv[1]), sprintf("%.2f", bulk),
             sprintf("%.2f", ratio), sigma_src, verdict),
    stringsAsFactors = FALSE)

  advice <- c(
    if (need_firth) "① 最小死亡數 < 5 → Firth 為必要（防發散與幼年 alpha 偏誤）"
    else            "① 最小死亡數 >= 5 → 卜瓦松 MLE 即可",
    if (max(d_med) / pmax(min(d_med), 1) > 10)
      "① 跨年齡落差 > 10 倍 → 不可用等權重 SVD，改卜瓦松或 w = mu_hat 的加權 SVD" else NULL,
    if (any(zero_by_age > 0.2))
      sprintf("② 下列年齡的 alpha 在未修正方法下不可用：%s",
              paste(rownames(D)[zero_by_age > 0.2], collapse = ", ")) else NULL,
    if (het_span > 2) "③ 噪音變異跨度 > 2 個數量級 → 等權重方法被誤設" else NULL,
    sprintf("④ %s", verdict),
    "區間：卜瓦松 + 狀態空間 REML（無條件建議）；s1/上緣 < 1 時再加 Kalman 初始值變異數")

  list(A = A, Tn = Tn, c = cc, age_mid = age_mid, D = D, logm = logm,
       d_med = d_med, zero_by_age = zero_by_age, noise_var = noise_var,
       mu_med = if (is.null(fitter)) pmax(d_med, 0.5) else apply(mu, 1, median),
       sv = sv, sigma = sigma, bulk = bulk, ratio = ratio,
       summary = summary, advice = advice)
}

## ===================== 繪圖 ================================================
#' 四格診斷圖。可傳入多個 lc_eda() 結果做對照（最多三組）。
plot_lc_eda <- function(..., labels = NULL) {
  dgs <- list(...)
  if (is.null(labels)) labels <- paste0("資料", seq_along(dgs))
  stopifnot(length(dgs) <= 3)
  ## 有序量（資料量大小）用單一色相的 ordinal ramp，深 = 資料多
  pal <- c("#0d366b", "#2a78d6", "#86b6ef")[seq_along(dgs)]

  op <- par(mfrow = c(1, 4), mar = c(4.2, 4.4, 3.2, 1.2), bg = "#fcfcfb",
            col.axis = "#52514e", col.lab = "#52514e", fg = "#d6d5d0",
            cex.main = 1.05, font.main = 1)
  on.exit(par(op))
  grid_ <- function() grid(NA, NULL, col = "#e8e7e2", lty = 1, lwd = 0.7)

  ## ① 死亡數量級
  yl <- range(unlist(lapply(dgs, function(d) pmax(d$d_med, 0.4))))
  plot(NA, xlim = range(dgs[[1]]$age_mid), ylim = yl, log = "y",
       xlab = "年齡", ylab = "單年死亡數（中位數）",
       main = "① 各年齡的死亡數量級", col.main = "#0b0b0b")
  grid_()
  rect(par("usr")[1], par("usr")[3], par("usr")[2], log10(5),
       col = "#eb68341a", border = NA)
  abline(h = 5, col = "#b8b7b2", lty = 3)
  for (i in seq_along(dgs))
    lines(dgs[[i]]$age_mid, pmax(dgs[[i]]$d_med, 0.4), col = pal[i], lwd = 2, type = "b", pch = 16, cex = .5)
  legend("topleft", labels, col = pal, lwd = 2, bty = "n", cex = .85)

  ## ② 零格落在哪些年齡
  plot(NA, xlim = range(dgs[[1]]$age_mid), ylim = c(0, 100),
       xlab = "年齡", ylab = "零死亡年份的比例（%）",
       main = "② 零格落在哪些年齡", col.main = "#0b0b0b")
  grid_()
  for (i in seq_along(dgs))
    lines(dgs[[i]]$age_mid, 100 * dgs[[i]]$zero_by_age, col = pal[i], lwd = 2, type = "b", pch = 16, cex = .5)
  legend("topright", labels, col = pal, lwd = 2, bty = "n", cex = .85)

  ## ③ 噪音變異 vs 死亡數
  xs <- unlist(lapply(dgs, function(d) pmax(d$mu_med, 0.4)))
  ys <- unlist(lapply(dgs, function(d) pmax(d$noise_var, 1e-5)))
  plot(NA, xlim = range(xs), ylim = range(ys), log = "xy",
       xlab = "該年齡的期望死亡數", ylab = "由資料估得的噪音變異",
       main = "③ 異質變異：噪音隨死亡數變動", col.main = "#0b0b0b")
  grid_()
  xx <- 10^seq(log10(min(xs)), log10(max(xs)), length.out = 50)
  lines(xx, 1 / xx, col = "#7d7c78", lty = 2, lwd = 1.4)
  for (i in seq_along(dgs))
    points(pmax(dgs[[i]]$mu_med, 0.4), pmax(dgs[[i]]$noise_var, 1e-5),
           col = pal[i], pch = 16, cex = .85)
  legend("topright", c("理論 1/D", labels), col = c("#7d7c78", pal),
         lty = c(2, rep(NA, length(dgs))), pch = c(NA, rep(16, length(dgs))),
         lwd = c(1.4, rep(NA, length(dgs))), bty = "n", cex = .85)

  ## ④ 訊號是否穿出噪音
  K <- 6
  yl4 <- range(unlist(lapply(dgs, function(d) c(d$sv[1:K], d$bulk))))
  plot(NA, xlim = c(0.4, K + 2.6), ylim = yl4, log = "y",
       xlab = "第 k 個奇異值", ylab = "奇異值",
       main = "④ 訊號是否穿出噪音", col.main = "#0b0b0b")
  grid_()
  for (i in seq_along(dgs)) {
    lines(1:K, dgs[[i]]$sv[1:K], col = pal[i], lwd = 2, type = "b", pch = 16, cex = .6)
    abline(h = dgs[[i]]$bulk, col = pal[i], lwd = 1.3, lty = 2)
  }
  ## 判讀靠的是 s1 與上緣的「比值」，故把比值直接寫進圖例，不在圖面標註
  legend("bottomleft",
         sprintf("%s：s1/上緣 = %.2f", labels, vapply(dgs, `[[`, 0, "ratio")),
         col = pal, lwd = 2, pch = 16, bty = "n", cex = .78)
  mtext("實線 = 觀測奇異值；虛線 = 同尺度純噪音能產生的上緣", side = 1, line = 2.9,
        adj = 0, cex = .62, col = "#52514e")
  invisible(NULL)
}

###############################################################################
# 判讀門檻速查
#
#  ① 最小單年死亡中位數 < 5  -> Firth 為必要；< 20 -> 建議
#     跨年齡落差 > 10 倍       -> 不可用等權重 SVD
#  ② 零格比例 > 20% 的年齡    -> 該年齡的 alpha 在未修正方法下不可用
#  ③ 噪音變異跨度 > 2 個數量級 -> 等權重方法被誤設
#     點系統性高於 1/D 線       -> 過度離散，考慮負二項
#  ④ s1 / 噪音上緣
#       > 2      訊號充分，SVD 可用
#       1.2-2    可還原：卜瓦松或加權 SVD(w=mu_hat) + FH
#       0.6-1.2  邊緣：方法選擇效益最大
#       < 0.6    不可還原：不要估 beta
#
#  注意：s1 的絕對值不可判讀。小人口的 s1 反而較大（噪音灌水），
#        必須與上緣相比。本檔測試資料：
#          全國 1,187 萬  s1=3.44  上緣=0.53  比值=6.53
#          模擬 20 萬     s1=4.64  上緣=4.04  比值=1.15
#          模擬 5 萬      s1=4.70  上緣=8.36  比值=0.56
#
#  上緣公式 sigma*sqrt(T)*(1+sqrt(A/T)) 來自 Benaych-Georges & Nadakuditi
#  (2012), JMVA 111: 120-135 式(9)：在 BBP 相變門檻上，觀測到的最大奇異值
#  恰好等於純噪音矩陣的上緣。
###############################################################################
