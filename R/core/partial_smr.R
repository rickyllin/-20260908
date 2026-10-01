###############################################################################
# Partial SMR 與 Whittaker ratio 修勻（Lee 2003；Yue, Wang and Wang 2019）
#
#   本檔實作老師先前文章所採用的參考母體修勻法，作為本文中心化校正的對照組。
#   兩法的共同前提是存在一個人數較多、死亡率型態與目標小區域相近的參考母體；
#   本文的校正則不需要參考母體。兩者因此不是同一類方法，比較的意義在於
#   「以外部資訊換取穩定」與「以分配假設換取穩定」何者在臺灣資料上較有利。
#
#   Partial SMR（Lee 2003, eq. 2-3 of Yue et al. 2019）：
#
#     SMR   = sum_x d_x / sum_x e_x ,          e_x = P_x * mR_x
#     h2    = max( ( sum_x (d_x - e_x SMR)^2 - sum_x d_x ) /
#                  ( SMR^2 * sum_x e_x^2 ), 0 )
#     v_x   = mR_x * exp( ( d_x h2 log(d_x/e_x) + (1 - d_x/sum d) log SMR ) /
#                         ( d_x h2 + (1 - d_x/sum d) ) )
#
#   d_x = 0 時第一項以 0 計（0 * log 0 := 0），此時 v_x = SMR * mR_x，
#   亦即完全借用參考母體的年齡型態。零格因此不需要 c0 替代值，
#   這正是老師在 0929 會議提出的問題：若以 PSMR 處理零格，
#   是否也會得到與 c0 = 0.5 插補相同的實證結果。
#
#   Whittaker ratio（Yue et al. 2019, eq. 4）：對比值 s_x = m_x / mR_x 作
#   Whittaker 修勻，min_r sum_x w_x (r_x - s_x)^2 + h sum_x (Delta^z r_x)^2，
#   再以 r_x * mR_x 為修勻後死亡率。
###############################################################################

#' Partial SMR 修勻（單一年度）
#'
#' @param d  長度 A 的死亡數觀察值
#' @param P  長度 A 的曝露數
#' @param mR 長度 A 的參考母體中央死亡率
#' @return 長度 A 的修勻後死亡率
psmr_year <- function(d, P, mR, floor_m = 1e-12) {
  ## 注意：修勻後的死亡率 SMR * mR_x 恆為正，故不需要 c0 之類的替代值。
  ## 這正是參考母體修勻相對於零格替代的優勢，floor_m 只作數值保護。
  e  <- P * mR
  Sd <- sum(d); Se <- sum(e)
  if (Sd <= 0 || Se <= 0) return(pmax(mR, floor_m))
  SMR <- Sd / Se
  h2  <- max((sum((d - e * SMR)^2) - Sd) / (SMR^2 * sum(e^2)), 0)

  w_obs <- d * h2                       # 觀測的權重：死亡數愈多、異質性愈大則愈重
  w_ref <- 1 - d / Sd                   # 參考的權重：死亡數愈多則愈輕
  lo    <- ifelse(d > 0, log(pmax(d, .Machine$double.xmin) / pmax(e, 1e-300)), 0)
  num   <- w_obs * lo + w_ref * log(SMR)
  den   <- w_obs + w_ref
  r     <- ifelse(den > 0, num / den, log(SMR))
  pmax(mR * exp(r), floor_m)
}

#' Partial SMR 修勻（整個 A x T 矩陣，逐年施作）
psmr_matrix <- function(D, E, MR, floor_m = 1e-12) {
  out <- matrix(NA_real_, nrow(D), ncol(D), dimnames = dimnames(D))
  for (j in seq_len(ncol(D)))
    out[, j] <- psmr_year(D[, j], E[, j], MR[, j], floor_m)
  out
}

#' Whittaker 修勻（對任一序列，z 階差分懲罰）
whittaker_1d <- function(y, w, h, z = 2) {
  A <- length(y)
  Dm <- diag(A)
  for (i in seq_len(z)) Dm <- diff(Dm)
  solve(diag(w) + h * crossprod(Dm), w * y)
}

#' Whittaker ratio 修勻（單一年度）
whit_ratio_year <- function(d, P, mR, h = NULL, z = 2, c0 = 0.5) {
  ## Whittaker 修勻可能給出非正的比值。此處採與標準 LC 相同的零格慣例，
  ## 即修勻後的死亡率不得低於 c0 / P_x；如此三個方法在零格上的處理一致，
  ## 而 Partial SMR 則因 SMR * mR_x 恆為正而完全不需要這道下限。
  s <- (d / P) / pmax(mR, 1e-300)
  if (is.null(h)) h <- mean(P)
  r <- whittaker_1d(s, w = P, h = h, z = z)
  pmax(r * mR, c0 / pmax(P, 1e-300))
}

whit_ratio_matrix <- function(D, E, MR, h = NULL, z = 2, c0 = 0.5) {
  out <- matrix(NA_real_, nrow(D), ncol(D), dimnames = dimnames(D))
  for (j in seq_len(ncol(D)))
    out[, j] <- whit_ratio_year(D[, j], E[, j], MR[, j], h, z, c0)
  out
}

#' 修勻後再配適 Lee-Carter（Yue et al. 2019 的「先修勻、後配模型」路線）
lc_psmr <- function(D, E, MR, floor_m = 1e-12) {
  lc_svd_fit(log(psmr_matrix(D, E, MR, floor_m)))
}

lc_whit_ratio <- function(D, E, MR, h = NULL, z = 2, c0 = 0.5) {
  lc_svd_fit(log(whit_ratio_matrix(D, E, MR, h, z, c0)))
}

#' 只以 SMR 填補零格，其餘各格保留觀測死亡率
#'
#'   這是 0929 會議中老師提出的問題：若把零格的 c0 = 0.5 替代改成以
#'   標準化死亡比填補，是否會得到與 c0 插補相同的實證結果。
#'   零格填入 SMR_t * mR_{x,t}，即「若該年齡的死亡率與參考母體成同一比例，
#'   應有的死亡率」，其餘各格不動。
zero_fill_smr <- function(D, E, MR) {
  out <- D / E
  for (j in seq_len(ncol(D))) {
    e   <- E[, j] * MR[, j]
    SMR <- sum(D[, j]) / sum(e)
    z   <- D[, j] == 0
    if (any(z)) out[z, j] <- SMR * MR[z, j]
  }
  out
}

lc_zero_smr <- function(D, E, MR) lc_svd_fit(log(pmax(zero_fill_smr(D, E, MR), 1e-12)))

###############################################################################
# 死亡率比值情境（Yue et al. 2019, Figure 2 的七種情境）
#   s_x = m^small_x / m^ref_x；參考母體的真實死亡率即 m^small_x / s_x。
###############################################################################
ratio_scenarios <- function(A) {
  half <- ceiling(A / 2)
  v    <- c(seq(1.5, 0.5, length.out = half),
            seq(0.5, 1.5, length.out = A - half + 1)[-1])
  rv   <- c(seq(0.5, 1.5, length.out = half),
            seq(1.5, 0.5, length.out = A - half + 1)[-1])
  list("s=0.8"    = rep(0.8, A),
       "s=1.0"    = rep(1.0, A),
       "s=1.2"    = rep(1.2, A),
       "遞增"     = seq(0.5, 1.5, length.out = A),
       "遞減"     = seq(1.5, 0.5, length.out = A),
       "V 型"     = v,
       "倒 V 型"  = rv)
}

###############################################################################
# 標準 LC 的第二階段 kappa 重估
#   Lee and Carter (1992) 原文在奇異值分解之後重解 kappa_t，使配適的總死亡數
#   等於觀測總死亡數；Lee and Miller (2001) 改為使配適的零歲平均餘命等於
#   觀測值。零格存在時觀測死亡率須以 c0 替代，故一併以 c0 參數化。
###############################################################################
lc_stage2_deaths <- function(a, b, k, D, E, span = 60) {
  vapply(seq_along(k), function(j) {
    Dt <- sum(D[, j])
    if (Dt <= 0) return(k[j])
    tryCatch(uniroot(function(kk) sum(E[, j] * exp(a + b * kk)) - Dt,
                     c(k[j] - span, k[j] + span), extendInt = "yes",
                     tol = 1e-10)$root, error = function(e) k[j])
  }, numeric(1))
}

lc_stage2_e0 <- function(a, b, k, D, E, c0 = 0.5, span = 60) {
  mobs <- pmax(D, c0) / E
  e0o  <- apply(mobs, 2, function(m) life_table(m)$e0)
  out  <- adjust_kappa(a, b, k, function(m) life_table(m)$e0, e0o, span = span)
  ifelse(is.finite(out), out, k)
}

#' 標準 LC（含第二階段調整）
lc_two_stage <- function(D, E, c0 = 0.5, stage2 = c("none", "deaths", "e0")) {
  stage2 <- match.arg(stage2)
  f <- lc_svd_fit(log(pmax(D, c0) / E))
  k <- switch(stage2,
              none   = f$k,
              deaths = lc_stage2_deaths(f$a, f$b, f$k, D, E),
              e0     = lc_stage2_e0(f$a, f$b, f$k, D, E, c0))
  k <- k - mean(k)
  list(a = f$a, b = f$b, k = k, stage2 = stage2)
}
