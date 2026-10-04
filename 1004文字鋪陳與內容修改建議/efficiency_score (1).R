###############################################################################
# 校正後最小平方相對於 Poisson MLE 的效率，與 Neyman 正交性診斷
#
#   對應《五步驟架構與你的論文》一文的第肆節
#
#   主要結果（本檔已數值驗證至 1e-10）：
#
#     令 g(D) = log max(D, c0)，D ~ Poisson(mu)，連結為 mu = E exp(eta)。
#     校正後估計方程  psi = g(D) - log E - eta - b(mu(eta))  滿足 E[psi] = 0。
#
#     靈敏度   A = -E[d psi / d eta] = 1 + b'(mu) mu = mu * m1'(mu) = Cov(g(D), D)
#     變異數   B = Var(g(D))
#     三明治   Var(eta_hat) = B / A^2
#
#     相對效率 = (1/mu) / (B/A^2) = Cov(g,D)^2 / (Var(g) Var(D)) = corr(g(D), D)^2
#
#   這條恆等式用到 Poisson 的 d/dmu E[h(D)] = Cov(h(D), D)/mu。
#
#   數值性質：效率在 mu -> 0 與 mu -> inf 兩端皆趨近 1，最低點在 mu = 4.47
#   處為 0.8776，亦即效率損失全程不超過約 12%，且**最差處在中段**而非資料最稀處。
#
#   另：對 mu_hat 相對誤差的一階敏感度 = Cov(g,D) - 1，在 mu -> 0 時趨近 -1。
#   這是「校正在最需要它的地方自動失效」的解析式。
###############################################################################

#' g(D) = log max(D, c0) 的前兩階動差與其與 D 的共變數
#' @return c(Eg, Vg, Cov_gD)
g_moments <- function(mu, c0 = 0.5) {
  if (mu <= 0) stop("mu must be positive")
  dmax <- ceiling(mu + 16*sqrt(mu) + 20)        # 截斷誤差遠小於機器精度
  d <- 0:dmax
  p <- dpois(d, mu)
  g <- ifelse(d == 0, log(c0), log(pmax(d, 1)))
  Eg <- sum(p*g)
  c(Eg = Eg, Vg = sum(p*g^2) - Eg^2, Cov = sum(p*g*d) - Eg*mu)
}

#' 偏誤函數 b(mu; c0) = E[g(D)] - log mu    （報告式 13）
b_bias <- function(mu, c0 = 0.5) unname(g_moments(mu, c0)["Eg"]) - log(mu)

#' 靈敏度 A = Cov(g(D), D)。亦等於 1 + b'(mu)*mu，本檔的 verify_identity() 檢查之
sensitivity_A <- function(mu, c0 = 0.5) unname(g_moments(mu, c0)["Cov"])

#' 校正後最小平方相對於 Poisson MLE 的效率 = corr(g(D), D)^2
#' @note 這是 oracle 效率（b 中的 mu 視為已知）。插入 mu_hat 的額外代價見 plugin_sensitivity()
rel_efficiency <- function(mu, c0 = 0.5) {
  m <- g_moments(mu, c0)
  unname(m["Cov"]^2 / (m["Vg"] * mu))
}

#' 對 mu_hat 相對誤差的一階敏感度：mu_hat = mu(1+eps) 使 alpha_hat 偏移約 此值 * eps
#' 不為零即表示估計方程**不是 Neyman 正交**的
plugin_sensitivity <- function(mu, c0 = 0.5) sensitivity_A(mu, c0) - 1

## ========================= 驗證 ============================================

verify_identity <- function(c0 = 0.5, h = 1e-6) {
  mus <- c(0.05,0.1,0.3,0.5,log(2),0.9234,1,2,3,5,10,20,50,100)
  A1 <- sapply(mus, function(m) 1 + (b_bias(m+h,c0)-b_bias(m-h,c0))/(2*h)*m)
  A2 <- sapply(mus, sensitivity_A, c0 = c0)
  cat(sprintf("恆等式 A = 1 + b'(mu)mu = Cov(g,D)：最大絕對差 %.2e\n",
              max(abs(A1-A2))))
  invisible(data.frame(mu = mus, A_numeric = A1, A_closed = A2))
}

#' 效率曲線與最低點
efficiency_profile <- function(c0 = 0.5) {
  lo <- optimize(function(m) rel_efficiency(m, c0), c(0.3, 40))
  cat(sprintf("效率最低點：mu = %.4f，相對效率 = %.4f（損失 %.1f%%）\n",
              lo$minimum, lo$objective, 100*(1-lo$objective)))
  f <- function(m) rel_efficiency(m, c0) - 0.9
  a <- tryCatch(uniroot(f, c(0.01, lo$minimum))$root, error = function(e) NA)
  b <- tryCatch(uniroot(f, c(lo$minimum, 400))$root, error = function(e) NA)
  cat(sprintf("效率 < 0.9 的區間：mu in (%.3f, %.3f)\n", a, b))
  invisible(list(min_at = lo$minimum, min_eff = lo$objective, below90 = c(a,b)))
}

#' 在實際的 年齡x年份 期望死亡數矩陣上評估
#' @param mu  年齡 x 年份 的期望死亡數矩陣（= E * m_hat）
lc_efficiency_report <- function(mu, c0 = 0.5, age_labels = rownames(mu)) {
  if (is.null(age_labels)) age_labels <- paste0("age", seq_len(nrow(mu)))
  re <- matrix(sapply(as.vector(mu), rel_efficiency, c0 = c0), nrow(mu))
  se <- matrix(sapply(as.vector(mu), plugin_sensitivity, c0 = c0), nrow(mu))
  cat(sprintf("oracle 相對效率：格平均 %.4f，資訊加權 %.4f，最低 %.4f\n",
              mean(re), weighted.mean(re, mu), min(re)))
  cat(sprintf("插入敏感度 |Cov(g,D)-1|：最大 %.4f（%s）\n\n",
              max(abs(se)), age_labels[which.max(abs(rowMeans(se)))]))
  data.frame(
    age         = age_labels,
    mu_mean     = rowMeans(mu),
    b           = sapply(rowMeans(mu), b_bias, c0 = c0),
    rel_eff     = rowMeans(re),
    plugin_sens = rowMeans(se),
    row.names = NULL)
}

###############################################################################
# 已驗證的參考值（c0 = 0.5）
#
#   mu       b(mu)      Cov(g,D)   相對效率   插入敏感度
#   0.05    +2.3372      0.0346     0.9999     -0.9654
#   0.30    +0.7176      0.2046     0.9973     -0.7954
#   0.6931  +0.1416      0.4468     0.9863     -0.5532   <- 中位數崩潰門檻
#   0.9234   0.0000      0.5691     0.9773     -0.4309   <- 偏誤變號點
#   2.00    -0.1874      0.9445     0.9283     -0.0555
#   4.47    -0.1413      1.1240     0.8776     +0.1240   <- 效率最低點
#   10.0    -0.0552      1.0618     0.9347     +0.0618
#   100     -0.0050      1.0051     0.9949     +0.0051
#
# 讀法：
#   b 在 mu 小處最大（偏誤問題）
#   相對效率在 mu 中段最低（效率問題）
#   插入敏感度在 mu 小處最大（正交性問題）
#   三者的最差處互不重疊，這是論文可以強調的結構。
###############################################################################
