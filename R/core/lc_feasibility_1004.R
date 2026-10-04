###############################################################################
# 小人口 LC 的事前可行性檢定
#
#   用途：在做任何估計、任何模擬之前，只憑「曝露數 + 一組粗略死亡率」
#         就先算出該地區的 LC 參數會偏多少，並給出分級建議。
#
#   依據：進度報告（1005）Proposition 1 / Corollary 1
#         E[alpha_hat_x] - alpha_x = (1/T) sum_t b(mu_{x,t}; c0)      恆等式
#         b(mu;c0) = e^{-mu} log c0 + sum_{d>=1} e^{-mu} mu^d/d! log d - log mu
#
#   這是閉式，不需要模擬、不需要參考母體、不需要真值。
#   報告第肆節之三已驗證：三種差異顯著的年齡結構下，閉式對逐年齡偏誤的
#   相關係數為 0.998 / 0.999 / 1.000。
#
#   用法：
#     fs <- lc_feasibility(E_matrix, m_hat_matrix)   # 年齡 x 年份
#     print(fs)
###############################################################################

## ========================= 偏誤函數 =========================================

#' b(mu; c0)：標準 LC 估計方程在真值處的期望（式 13）
#' @param mu 期望死亡數（可為向量）
#' @param c0 零格替代值。慣例 0.5，但見 mu_star() —— 它會改變偏誤的方向
b_bias <- function(mu, c0 = 0.5) {
  sapply(mu, function(m) {
    if (m <= 0) return(Inf)
    dmax <- ceiling(m + 12*sqrt(m) + 12)          # 截斷誤差 < 機器精度
    d <- 1:dmax
    exp(-m)*log(c0) + sum(dpois(d, m)*log(d)) - log(m)
  })
}

#' 偏誤的變號點：mu < mu* 時 alpha 偏高，mu > mu* 時偏低
#' 注意這個門檻由「分析者選的 c0」決定，不是資料的性質
mu_star <- function(c0 = 0.5)
  uniroot(function(m) b_bias(m, c0), c(1e-8, 50), tol = 1e-10)$root

#' 中位數／分位數類估計量的崩潰門檻：pi = e^{-mu} > 1/2
#' 與 c0 無關，是資料的性質（Donoho & Huber 1983 的 50% 上限代入混合表示）
MU_BREAKDOWN <- log(2)   # 0.6931

## ========================= 可行性檢定 =======================================

#' @param E      年齡 x 年份 的曝露數矩陣
#' @param m_hat  同維度的粗略死亡率估計。可用全國率、鄰近地區率、
#'               或該地區自己的 D/E（零格先以非零者內插）——閉式對這個起點不敏感，
#'               因為它只透過 mu = E*m 進入，而分級門檻以數量級為單位。
#' @param c0     零格替代值，須與後續實際估計時所用者一致
lc_feasibility <- function(E, m_hat, c0 = 0.5, age_labels = rownames(E)) {
  stopifnot(all(dim(E) == dim(m_hat)))
  A <- nrow(E); Tn <- ncol(E)
  if (is.null(age_labels)) age_labels <- paste0("age", seq_len(A))

  mu <- E * m_hat                                   # 各格期望死亡數
  B  <- matrix(b_bias(as.vector(mu), c0), A, Tn)
  bx <- rowMeans(B)                                 # Corollary 1：alpha 的偏誤
  ms <- mu_star(c0)

  per_age <- data.frame(
    age          = age_labels,
    mu_mean      = rowMeans(mu),
    mu_min       = apply(mu, 1, min),
    pi_zero      = rowMeans(exp(-mu)),              # 預期零格比例
    alpha_bias   = bx,
    rate_factor  = exp(bx),                         # 死亡率被放大幾倍
    ## ASCII 代碼：print(data.frame) 在非 UTF-8 locale 下會把中文變義序列
    status       = ifelse(rowMeans(mu) < MU_BREAKDOWN, "BRK",   # 中位數已崩潰
                   ifelse(rowMeans(mu) < ms,           "POS",   # 偏誤為正
                                                       "neg")), # 偏誤為負
    stringsAsFactors = FALSE, row.names = NULL)

  f_break <- mean(mu < MU_BREAKDOWN)
  f_pos   <- mean(mu < ms)
  tier <- if (f_pos < 0.02) "A 可直接編表"      else
          if (f_pos < 0.15) "B 建議校正"        else
          if (f_pos < 0.35) "C 必須校正"        else
                            "D 不建議單獨編表"

  advice <- switch(substr(tier,1,1),
    A = c("標準 LC 的 alpha 偏誤已可忽略。",
          "仍建議加資訊加權 w = mu_hat：它處理的是 beta 的異質變異，與人口規模無關。"),
    B = c("建議採用「僅校正 alpha + 加權」：兩者正交可疊加。",
          "校正只施加在列平均（alpha），加權只施加在中心化矩陣（beta, kappa）。"),
    C = c("必須校正。優先序：先加權（對 e0 影響最大），再校正 alpha（對參數本身影響最大）。",
          "若已有型態相近的參考母體，可先以 SMR*m_ref 填補零格，再施以中心化。",
          "請同時檢查秩一結構是否成立——校正依賴 mu_hat，而 mu_hat 由秩一配適算出。"),
    D = c("不建議單獨編製生命表。閉式顯示偏誤已大到任何校正都只是在偏誤與變異之間重分配。",
          "可行的替代：與鄰近地區合併、採多母體共同因子模型、或只報告經年齡標準化的彙總指標。"))

  if (f_break > 0.10)
    advice <- c(advice,
      sprintf("有 %.0f%% 的格子已超過中位數崩潰門檻（期望死亡數 < %.3f）：",
              100*f_break, MU_BREAKDOWN),
      "  分位數迴歸、中位數型穩健估計量在這些年齡上不適用，不是調參數可以救的。")

  structure(list(
    A = A, Tn = Tn, c0 = c0, mu_star = ms,
    frac_breakdown = f_break, frac_positive_bias = f_pos,
    alpha_bias_median = median(bx), alpha_bias_max = max(abs(bx)),
    worst_age = age_labels[which.max(abs(bx))],
    per_age = per_age, tier = tier, advice = advice),
    class = "lc_feas")
}

print.lc_feas <- function(x, ...) {
  cat(sprintf("小人口 LC 可行性檢定   A=%d  T=%d  c0=%.2f\n", x$A, x$Tn, x$c0))
  cat(sprintf("  偏誤變號點 mu* = %.4f（由 c0 決定）；崩潰門檻 ln2 = %.4f（資料性質）\n",
              x$mu_star, MU_BREAKDOWN))
  cat(sprintf("  期望死亡數 < mu* 的格子：%.1f%%；< ln2 的格子：%.1f%%\n",
              100*x$frac_positive_bias, 100*x$frac_breakdown))
  cat(sprintf("  預測 alpha 偏誤：中位 %+.4f，最大 %+.4f（%s）\n",
              x$alpha_bias_median, x$alpha_bias_max, x$worst_age))
  cat(sprintf("\n  分級：%s\n\n", x$tier))
  for (a in x$advice) cat("  -", a, "\n")
  cat("\n逐年齡（rate_factor = exp(alpha 偏誤) = 死亡率被放大的倍數）\n")
  cat("  status: BRK = 中位數型估計量已崩潰 | POS = 偏誤為正 | neg = 偏誤為負\n")
  p <- x$per_age
  p$mu_mean <- round(p$mu_mean,2); p$mu_min <- round(p$mu_min,2)
  p$pi_zero <- round(p$pi_zero,3); p$alpha_bias <- round(p$alpha_bias,4)
  p$rate_factor <- round(p$rate_factor,3)
  print(p, row.names = FALSE)
  invisible(x)
}

###############################################################################
# 分級門檻的來源（進度報告 1005，臺灣女性 2001-2024，重複 100 次）
#
#   N        <mu* 比例   標準 LC 的 alpha 最大偏誤   e0 偏誤     分級
#   10,000     42%              2.33              -3.62 歲     D
#   50,000     22%              0.87              -1.48 歲     C
#   200,000     4%              0.18              -0.90 歲     B
#   1,000,000   0%              ~0                 ~0          A
#
# 重要：e0 的改善主要來自「加權」而非「校正 alpha」。報告表 17 顯示
#   僅校正 alpha：      alpha 最大偏誤 0.8658 -> 0.0439，但 e0 僅 -1.480 -> -1.466
#   僅校正 alpha + 加權：e0 -1.466 -> +0.007
# 故若只能實作一項，實務上應先做加權（w = mu_hat），再做校正。
# 校正的價值在於 alpha 參數本身（跨區比較、作為共變數時）。
###############################################################################
