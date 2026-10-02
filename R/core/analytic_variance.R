###############################################################################
# 對數轉換的二階解析量：變異數、偏誤的導數，與校正所付的變異數代價
#
# 動機。`analytic_bias_lc.R` 只處理一階量 b(mu; c0) = E[log max(D,c0)] - log mu，
# 亦即估計方程期望的偏離。但 1005 版的實證另留下兩件未解釋的事：
#   (i)  標準 LC 的 alpha 帶在幼年組「異常地窄」（人數五萬時帶寬中位 0.218，
#        而該處偏誤達 0.866）；
#   (ii) 中心化之後帶反而變寬（0.218 -> 0.277），且幼年組變寬得更多。
# 這兩件都是二階（變異數）現象，故須另立閉式。本檔提供三個量。
#
# 一、對數尺度的變異數
#
#   v(mu; c0) = Var[ log max(D, c0) ],   D ~ Poisson(mu).
#
#   小 mu 時 log max(D, c0) 以機率 e^{-mu} 取 log c0、以機率 mu e^{-mu} 取 0，
#   故 v ~ mu (log c0)^2 -> 0。這正是 (i) 的解析解釋：零格主導的年齡，
#   其對數死亡率幾乎必然等於常數 log(c0/E)，重複之間幾無變動，
#   因而帶窄；但該常數離真值很遠，因而偏誤大。
#   「帶窄」與「準」在此完全脫鉤。
#
# 二、偏誤函數的導數（精確式，非數值微分）
#
#   對任意 g，d/dmu E[g(D)] = E[g(D+1)] - E[g(D)]。取 g = log max(., c0)，
#   並注意 D+1 >= 1 使 max 失效，得
#
#     b'(mu; c0) = E[log(D+1)] - E[log max(D, c0)] - 1/mu.                (**)
#
#   兩端的極限皆可直接讀出：
#     mu -> 0 時 E[log(D+1)] -> 0、E[log max(D,c0)] -> log c0，
#       故 mu * b'(mu) -> -1；
#     mu -> inf 時 E[log(D+1)] - E[log D] ~ 1/mu，故 mu * b'(mu) -> 0。
#
# 三、校正的變異數膨脹因子
#
#   校正須代入 mu_hat，而 mu_hat = E exp(alpha_hat + beta_hat kappa_hat) 依賴
#   alpha_hat 本身。令 A = alpha_hat_x - alpha_x，則一階展開給出
#
#     b_bar(mu_hat) = b_bar(mu) + c_x A + (來自 beta_hat kappa_hat 的項),
#     c_x = mean_t [ mu_{x,t} b'(mu_{x,t}; c0) ].                         (***)
#
#   若暫時忽略後一項（即把 beta kappa 當已知），校正後的偏差為
#   (1 - c_x) A - b_bar_x，故變異數膨脹因子為 (1 - c_x)^2。
#   由 (**) 的極限，零格主導的年齡 c_x -> -1，膨脹因子 -> 4；
#   期望死亡數大的年齡 c_x -> 0，膨脹因子 -> 1。
#   亦即\textbf{校正所付的變異數代價，恰好集中在它收益最大的年齡}。
#
# 四、以臺灣女性資料核對（REPS = 100，種子同 run_analytic_bias.R）
#
#   v(mu; c0) 對標準 LC 逐年齡 alpha 的抽樣標準差：閉式與實測的比值在
#   22 個年齡上落在 0.82 至 1.15 之間，兩個規模的中位數分別為 0.97 與 1.02。
#   亦即\textbf{標準 LC 的 alpha 變異數亦可由單一個死亡數的函數預測}，
#   與偏誤的情形相同。
#
#   膨脹因子 (1 - c_x)^2 的核對結果分為兩段，須分別陳述。
#     mu >~ 1 的年齡（N = 5e4 時為第 7 至 22 組）：閉式與實測幾乎重合，
#       例如第 13 至 21 組的預測為 0.92/0.95/0.97/0.98/0.98/0.98/0.99/0.98/0.95，
#       實測為 0.92/0.95/0.97/0.98/0.98/0.98/0.99/0.98/0.95。
#       值得注意的是該值\textbf{小於 1}：mu > 2 時 b'(mu) > 0 使 c_x > 0，
#       故校正在期望死亡數大的年齡反而\textbf{降低}變異數。
#     mu < 1 的年齡：閉式\textbf{高估}膨脹，例如 N = 5e4 的第 3 組預測 3.33
#       而實測 2.54。原因是一階展開在該處失效，b 於 mu 小時高度凸，
#       且 mu_hat 受零格替代的牽制而不如展開式所假設的那樣自由變動。
#   因此 (***) 不是上界也不是下界，而是在 mu >~ 1 時精確、在 mu < 1 時偏高。
#   二階展開或直接以 mu_hat 的模擬分布計算，列為本分支的待辦。
#
#   另有一項與偏誤無關的發現：v(mu; c0) 對 mu 非單調，於 mu ~ 2 處取極大
#   （c0 = 0.5 時 v(2) = 0.480），兩端皆趨於零。
#   亦即對數尺度上最吵的年齡不是零格主導的年齡，而是每年約兩人死亡的年齡。
###############################################################################

#' E[ g(D) ] 的截斷精確計算，g 取 log max(., c0) 與其平方
#' @return list(m1 = E[log max(D,c0)], m2 = E[(log max(D,c0))^2],
#'              lp1 = E[log(D+1)])
pois_log_moments <- function(mu, c0 = 0.5) {
  vapply(mu, function(m) {
    if (!is.finite(m) || m <= 0) return(c(NA_real_, NA_real_, NA_real_))
    K  <- max(80L, as.integer(ceiling(m + 12 * sqrt(m) + 12)))
    d  <- 0:K
    p  <- dpois(d, m)
    p  <- p / sum(p)                      # 截尾後重新正規化
    lg <- ifelse(d == 0, log(c0), log(pmax(d, 1e-300)))
    c(sum(p * lg), sum(p * lg^2), sum(p * log(d + 1)))
  }, numeric(3))
}

#' v(mu; c0) = Var[ log max(D, c0) ]
logvar_exact <- function(mu, c0 = 0.5) {
  M <- pois_log_moments(mu, c0)
  pmax(M[2, ] - M[1, ]^2, 0)
}

#' b'(mu; c0)，依式 (**) 精確計算
logbias_deriv <- function(mu, c0 = 0.5) {
  M <- pois_log_moments(mu, c0)
  M[3, ] - M[1, ] - 1 / mu
}

#' mu * b'(mu)，即式 (***) 中逐格的貢獻；極限為 -1（mu->0）與 0（mu->inf）
logbias_elast <- function(mu, c0 = 0.5) mu * logbias_deriv(mu, c0)

#' 逐年齡的解析預測：alpha_hat 的變異數、校正的膨脹因子與校正後變異數下界
#' @param MU A x T 的期望死亡數矩陣
#' @return data.frame(age, sd_raw, c_x, inflate, sd_corr_lower)
alpha_var_predict <- function(MU, c0 = 0.5) {
  A <- nrow(MU); T <- ncol(MU)
  v  <- matrix(logvar_exact(as.vector(MU), c0), A, T)
  el <- matrix(logbias_elast(as.vector(MU), c0), A, T)
  sd_raw <- sqrt(rowSums(v)) / T          # Var(alpha_hat) = T^-2 sum_t v
  c_x    <- rowMeans(el)
  data.frame(age = seq_len(A), sd_raw = sd_raw, c_x = c_x,
             inflate = (1 - c_x)^2,
             sd_corr_lower = abs(1 - c_x) * sd_raw)
}
