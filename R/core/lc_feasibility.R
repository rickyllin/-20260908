###############################################################################
# 事前可行性檢定：在配適之前算出標準 LC 會偏多少
#
#   本檔為兩版的合併。分級判準與建議文字取自 R/core/lc_feasibility_1004.R
#   （1001 精算建議的附件），其判準為「期望死亡數低於變號點 mu* 的格子比例」；
#   另補上偏誤量級、兩條分支的分離，以及 e0 的 alpha 通道。
#
#   動機。閉式 b(mu; c0) = E[log max(D,c0)] - log mu 只需要兩樣東西：
#   曝露數 E_{x,t}，以及一組粗略的死亡率 m^R_{x,t}（全國表或縣市表即可）。
#   兩者在配適之前都已知，故標準 LC 的 alpha 偏誤可以\textbf{事前}算出，
#   不需要模擬、不需要真值、也不需要把資料配適過一次。
#
#   與信賴度理論的對照。古典的 full credibility 標準（期望索賠數 1082 件，
#   對應 k = 5%、p = 90%）控制的是估計的\textbf{變異}，其推導假設估計量無偏。
#   本檔提供的是缺的那一半：在同一個「期望次數」的尺度上給出\textbf{偏誤}。
#
#   分級。依期望死亡數低於 mu* 的格子比例分為四級。
#   之所以不以 max|b_bar| 分級，是因為 b 有兩條性質完全不同的分支：
#     正偏誤（mu < mu*）源自零格替代，隨 mu -> 0 無上界，是小人口特有的病灶；
#     負偏誤（mu > mu*）源自 Jensen 不等式，絕對值的上界為
#       |b(2.2886)| = 0.19104，與人口規模無關，任何取對數的模型都有。
#   年齡剖面上總有某一組的期望死亡數落在 b 取極小處，故 max|b_bar| 在大人口
#   時卡在 0.19 附近不再下降，以它分級會使 A 級幾乎無法達到。
#   以零格分支的比例分級則直接對應病灶本身。四級的門檻與對應的 N 為
#     A（< 2%）   可直接編表           N > 215,300
#     B（2-15%）  建議校正             95,800 - 215,300
#     C（15-35%） 必須校正             24,400 - 95,800
#     D（>= 35%） 不建議單獨編表       N < 24,400
#   （對應的 N 以臺灣女性 2024 年齡結構算出，見 R/studies/run_feasibility.R）
#
#   須明白寫出的限制。本檔預測的是 alpha 這一條通道。零歲平均餘命的偏誤另有
#   一條來自 kappa 壓縮的通道，不在 b(mu; c0) 的涵蓋範圍內，
#   故 `e0_alpha_channel` 一欄只是 e0 偏誤的一部分，不可當成 e0 偏誤的預測。
#   另須注意該通道的量隨 N 變化甚大（N = 5 萬時僅 +0.07 歲，
#   N = 1 萬時達 -1.64 歲），不可由單一規模外推。
###############################################################################

MU_BREAKDOWN <- log(2)    # 中位數型估計量的崩潰門檻，與 c0 無關

#' 偏誤變號點 mu*：b(mu; c0) = 0 的唯一正根
mustar <- function(c0 = 0.5)
  uniroot(function(m) logbias_exact(m, c0), c(1e-6, 50), tol = 1e-10)$root

#' 事前可行性檢定
#' @param E    A x T 曝露數矩陣，或長度 A 的向量（視為各年相同）
#' @param mref A x T 參考死亡率矩陣，或長度 A 的向量
#' @param c0   零格替代值，須與後續實際估計時所用者一致
#' @param ages 年齡組名稱
lc_feasibility <- function(E, mref, c0 = 0.5, ages = NULL) {
  if (is.null(dim(E)))    E    <- matrix(E, ncol = 1)
  if (is.null(dim(mref))) mref <- matrix(mref, nrow = nrow(E), ncol = ncol(E))
  stopifnot(all(dim(E) == dim(mref)))
  A <- nrow(E); Tn <- ncol(E)
  if (is.null(ages)) ages <- rownames(E)
  if (is.null(ages)) ages <- as.character(seq_len(A))

  MU <- E * mref
  ms <- mustar(c0)
  B  <- matrix(logbias_exact(as.vector(MU), c0), A, Tn)
  b_bar <- rowMeans(B)

  byage <- data.frame(
    age        = ages,
    mu_median  = apply(MU, 1, median),
    pi_zero    = rowMeans(exp(-MU)),     # 預期零格比例
    share_low  = rowMeans(MU < ms),      # 該年齡落在變號點之下的年份比例
    b_bar      = b_bar,
    m_ratio    = exp(b_bar),             # 死亡率被放大的倍數
    status     = ifelse(apply(MU, 1, median) < MU_BREAKDOWN, "崩潰",
                 ifelse(apply(MU, 1, median) < ms,           "偏高", "偏低")),
    stringsAsFactors = FALSE, row.names = NULL)

  f_pos   <- mean(MU < ms)              # 分級的判準
  f_break <- mean(MU < MU_BREAKDOWN)
  worst <- which.max(abs(b_bar))
  amax  <- abs(b_bar[worst])
  grade <- if (f_pos < 0.02) "A" else if (f_pos < 0.15) "B" else
           if (f_pos < 0.35) "C" else "D"
  advice <- switch(grade,
    A = "可直接編表；仍建議加資訊加權，其處理的是 beta 的異質變異",
    B = "建議採用僅校正 alpha 加資訊加權，兩者正交可疊加",
    C = "必須校正並揭露不確定性；若有型態相近的參考母體，可先以 SMR 填補零格",
    D = "不建議單獨編表（改採合併鄰近地區、多母體共同因子，或只報彙總指標）")
  driver <- if (b_bar[worst] > 0) "零格高估" else "Jensen 低估"

  # alpha 通道對 e0 的影響：把預測偏誤加到參考死亡率上重編生命表。
  # 這\textbf{不是} e0 偏誤的預測，kappa 壓縮的通道不在其中。
  mlast <- mref[, Tn]
  e0_bias <- life_table(mlast * exp(b_bar))$e0 - life_table(mlast)$e0

  list(byage = byage,
       summary = data.frame(
         N_total        = sum(E[, Tn]),
         share_below_mustar = f_pos,
         share_below_ln2    = f_break,
         alpha_max_pred = amax,
         worst_age      = ages[worst],
         b_pos_max      = max(b_bar),
         b_neg_min      = min(b_bar),
         driver         = driver,
         grade          = grade,
         advice         = advice,
         e0_alpha_channel = e0_bias,
         mustar         = ms,
         stringsAsFactors = FALSE))
}
