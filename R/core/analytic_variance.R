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
# 三、估計方程的靈敏度，與校正的變異數代價
#
#   校正後的估計方程為 psi = g_{c0}(D) - log E - eta - b(mu(eta))，其中
#   mu = E exp(eta)。該方程的期望為零（這正是校正的定義），其靈敏度為
#
#     A(mu) = -E[d psi / d eta] = 1 + mu b'(mu; c0).                     (***)
#
#   由 Poisson 的共變數恆等式 d/dmu E[h(D)] = Cov(h(D), D)/mu 另可得
#   A(mu) = Cov(g_{c0}(D), D)，故 A 同時是靈敏度與共變數，兩種讀法等價。
#
#   依 M-估計的三明治公式，校正後估計量的逐格變異數為 v(mu)/A(mu)^2，
#   而未校正者為 v(mu)。故\textbf{變異數的代價為 1/A^2}，其中
#   A 的兩端極限由 (**) 直接讀出：mu -> 0 時 A -> 0 故代價發散；
#   mu -> inf 時 A -> 1 故代價趨於 1。A 的極大在 mu = 4.4288 處為 1.1259，
#   對應代價的極小 0.7888，亦即\textbf{校正在期望死亡數大的年齡反而降低變異數}。
#
#   須注意這是\textbf{迭代（不動點）}版本的代價，亦即 R/core/analytic_bias_lc.R
#   的 lc_analytic() 與 lc_alpha_only() 所實作者：該處 mu_hat 由校正後的參數
#   重新算出，反覆至收斂。若只做一步校正（mu_hat 取自未校正的配適），
#   代價改為 (1 - mu b')^2，量級小得多。兩者不可混用，見第四節的實測。
#
# 四、相對於 Poisson 最大概似的效率
#
#   Poisson MLE 的逐格變異數為 1/mu，故校正後最小平方相對於它的效率為
#
#     RE(mu) = (1/mu) / (v(mu)/A^2) = A^2 / (mu v(mu))
#            = Cov(g,D)^2 / (Var(g) Var(D)) = corr(g_{c0}(D), D)^2.       (****)
#
#   由最後一個等式立得 0 <= RE <= 1。數值上兩端皆趨近 1，
#   極小值為 0.8776（mu = 4.4700），亦即\textbf{取對數與零格替代合起來，
#   其效率損失全程不超過約 12%}，而最差處在中段而非資料最稀處。
#   後者的理由是 mu 小時 D 幾乎只取 {0,1}，而 g_{c0} 在該兩點上是雙射，
#   相關係數因而趨近 1。
#
#   一項順帶的結論：g_{c0}(d) = log max(d, c0) 在 0 < c0 < 1 時是
#   {0,1,2,...} 上的雙射，故 sigma(g(D)) = sigma(D)，
#   \textbf{零格替代與取對數合起來不損失任何資訊}。標準 LC 的全部損失
#   因而來自估計方程的期望不為零，而非來自轉換丟掉了訊息。
#
# 五、以臺灣女性資料核對（REPS = 100，種子同 run_analytic_bias.R）
#
#   v(mu; c0) 對標準 LC 逐年齡 alpha 的抽樣標準差：閉式與實測的比值在
#   22 個年齡上落在 0.82 至 1.15 之間，兩個規模的中位數分別為 0.97 與 1.02。
#   亦即\textbf{標準 LC 的 alpha 變異數亦可由單一個死亡數的函數預測}。
#
#   1/A^2 對迭代校正的變異數代價：mu >~ 0.2 的年齡與實測相當接近
#   （N = 5e4 時第 3 至 9 組的預測為 32.7/31.0/8.45/3.52/1.98/1.28/0.92，
#   實測為 31.5/26.7/8.48/3.77/2.05/1.54/1.18）；mu < 0.1 的年齡則閉式
#   大幅高估（N = 1e4 第 3 組預測 797 而實測 97），因為一階展開在
#   A 趨近零處失效。
#
#   另有一項與偏誤無關的發現：v(mu; c0) 對 mu 非單調，於 mu ~ 2 處取極大
#   （c0 = 0.5 時 v(2.0365) = 0.4805），兩端皆趨於零。
#   亦即對數尺度上最吵的年齡不是零格主導的年齡，而是每年約兩人死亡的年齡。
#
#   四個量的最差處互不重疊：b 在 mu -> 0 最大、v 在 mu = 2.04 最大、
#   1/A^2 在 mu -> 0 發散而於 mu = 4.43 最小、RE 在 mu = 4.47 最小。

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

#' 估計方程的靈敏度 A(mu) = 1 + mu b'(mu) = Cov(g(D), D)
sens_A <- function(mu, c0 = 0.5) 1 + logbias_elast(mu, c0)

#' 迭代校正的變異數代價 1/A^2
var_cost <- function(mu, c0 = 0.5) 1 / sens_A(mu, c0)^2

#' 相對於 Poisson MLE 的效率 = corr(log max(D,c0), D)^2
rel_efficiency <- function(mu, c0 = 0.5)
  sens_A(mu, c0)^2 / (mu * logvar_exact(mu, c0))

#' 逐年齡的解析預測：alpha_hat 的變異數、校正的變異數代價與校正後的變異數
#' @param MU A x T 的期望死亡數矩陣
#' @return data.frame(age, sd_raw, c_x, inflate, sd_corr_lower)
alpha_var_predict <- function(MU, c0 = 0.5) {
  A <- nrow(MU); T <- ncol(MU)
  v  <- matrix(logvar_exact(as.vector(MU), c0), A, T)
  el <- matrix(logbias_elast(as.vector(MU), c0), A, T)
  sd_raw <- sqrt(rowSums(v)) / T          # Var(alpha_hat) = T^-2 sum_t v
  c_x    <- rowMeans(el)
  Ax     <- 1 + c_x                       # 靈敏度
  re     <- rowMeans(matrix(rel_efficiency(as.vector(MU), c0), A, T))
  data.frame(age = seq_len(A), sd_raw = sd_raw, c_x = c_x, A = Ax,
             cost = 1 / Ax^2,             # 迭代校正的變異數代價
             cost_onestep = (1 - c_x)^2,  # 一步校正的對照
             sd_corr = sd_raw / abs(Ax),
             rel_eff = re)
}
