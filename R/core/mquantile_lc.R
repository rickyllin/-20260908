###############################################################################
# M-分位數 Lee-Carter：把「穩健程度」由開關改成連續刻度
#
#   動機。變體 F（中位數 LC）與變體 A（標準 LC）是兩個極端：前者的崩潰點
#   為 50%、效率低，後者崩潰點為 0、效率高。第伍節之六的結果顯示，中位數
#   在零格比例低於 50% 時有效、高於 50% 時失效，因此真正該問的不是「要不要
#   穩健」，而是「要多穩健」。M-分位數（M-quantile）恰好提供這個刻度。
#
#   損失函數（Breckling and Chambers 1988；Bianchi et al. 2018）：
#
#       rho_{tau,c}(u) = 2 * w_tau(u) * psi_c(u)
#       w_tau(u)       = tau          若 u >  0
#                      = 1 - tau      若 u <= 0
#       psi_c(u)       = u^2 / 2                若 |u| <= c      （Huber）
#                      = c|u| - c^2/2           若 |u| >  c
#
#   兩個極限把既有的變體全部納入同一族：
#       c -> 0    rho 正比於 w_tau(u)|u| = 檢查函數  -> 分位數迴歸
#                 tau = 0.5 時即變體 F（中位數 LC）
#       c -> Inf  rho = 2 w_tau(u) u^2 / 2          -> 非對稱最小平方
#                 （expectile，Newey and Powell 1987）
#                 tau = 0.5 時即變體 A（標準 LC，對數尺度最小平方）
#       tau = 0.5, c 有限                            -> Huber M-估計（Huber 1964）
#
#   因此 (tau, c) 的掃描等於在變體 A 與變體 F 之間畫出一條連續路徑，
#   可直接檢驗「最適穩健程度隨人口規模改變」這項預測。
#
#   非對稱截點（Xu and Chen 2018）。零格替代把觀測推「高」而非推低，
#   故被扭曲的是右尾；對稱的 Huber 截點會以相同力道壓抑兩側。允許
#   c_neg != c_pos 可只壓抑被扭曲的一側，代價是多一個調參。
#
#   演算法。psi_c 的一階導數為 min(|u|,c)sign(u)，故 rho' (u) = W(u) u，
#   其中 W(u) = 2 w_tau(u) min(1, c/|u|)。目標函數因而可寫成迭代加權最小
#   平方（IRLS），這正是 Hunter and Lange (2000) 為分位數迴歸提出的 MM
#   演算法在 c > 0 時的可微版本：以二次函數優化（majorize）原損失，
#   每一步的加權最小平方有封閉解，故目標函數單調不增。
#   相對於變體 F 所用的線性規劃（Koenker and Park 1996 的內點法亦屬此類），
#   IRLS 的優點是雙線性結構下兩個區塊都有封閉解，無須呼叫線性規劃求解器。
###############################################################################

#' Huber 尺度的一致性常數 beta_c = E[psi_c(Z)^2]，Z ~ N(0,1)
#'
#'   psi_c(z) = z 若 |z| <= c，否則 c*sign(z)，故
#'     beta_c = E[Z^2 1{|Z|<=c}] + c^2 P(|Z|>c)
#'            = 2(Phi(c) - 1/2) - 2 c phi(c) + 2 c^2 (1 - Phi(c))
huber_beta <- function(c) {
  if (is.infinite(c)) return(1)
  2 * (pnorm(c) - 0.5) - 2 * c * dnorm(c) + 2 * c^2 * (1 - pnorm(c))
}

#' 由一致性方程解尺度（Huber 的 Proposal 2）
#'
#'   Zoubir et al. (2018, sec. 3.5) 依 Ollila (2016) 的 M-Lasso 估計方程，
#'   把迴歸與尺度定義為同一組零次梯度方程的解，其尺度方程為
#'     (1 / (N * 2 alpha)) sum_i chi( r_i / sigma ) = 1 / gamma,
#'   對 Huber 損失而言即等價於熟知的
#'     (1 / N) sum_i psi_c( r_i / sigma )^2 = beta_c.
#'   此式對 sigma 單調，故以不動點迭代求解：
#'     sigma^2 <- sigma^2 * (1 / (N beta_c)) sum_i psi_c(r_i/sigma)^2
#'
#'   與原先每次迭代以 MAD 重估尺度的差別：MAD 在常態下的效率僅約 37%，
#'   且與損失函數的截點無關；一致性方程所解出的尺度則與 c 相容，
#'   使 (beta, sigma) 成為同一個目標函數的聯合 M-估計量。
#'
#' @param W 各格的資料權重（資訊加權）；NULL 表示等權重
huber_scale <- function(r, c, W = NULL, s0 = NULL, maxit = 50, tol = 1e-8) {
  r <- as.vector(r); ok <- is.finite(r)
  if (!is.null(W)) { w <- as.vector(W)[ok] } else { w <- rep(1, sum(ok)) }
  r <- r[ok]
  if (!length(r)) return(1e-8)
  if (is.infinite(c)) return(max(sqrt(sum(w * r^2) / sum(w)), 1e-8))
  bc <- huber_beta(c)
  s  <- if (is.null(s0)) max(mad(r), 1e-8) else max(s0, 1e-8)
  for (i in seq_len(maxit)) {
    u   <- r / s
    psi <- pmin(pmax(u, -c), c)
    s2  <- s^2 * sum(w * psi^2) / (sum(w) * bc)
    sn  <- max(sqrt(s2), 1e-10)
    if (abs(sn - s) < tol * max(1, s)) { s <- sn; break }
    s <- sn
  }
  s
}

#' M-分位數的 IRLS 權重
#'
#' @param u   殘差矩陣
#' @param tau 非對稱參數；0.5 為對稱
#' @param cn  負殘差側的 Huber 截點（可為 Inf）
#' @param cp  正殘差側的 Huber 截點（可為 Inf）
#' @param eps Hunter and Lange (2000) 的擾動項。c -> 0 時 1/|u| 在 u = 0 附近
#'            發散，MM 演算法以 |u| + eps 取代 |u| 使優化函數處處可微，
#'            eps 隨迭代縮小即收斂至原問題的解。
mq_weight <- function(u, tau = 0.5, cn = Inf, cp = Inf, eps = 0) {
  pos <- u > 0
  wt  <- ifelse(pos, tau, 1 - tau)
  cc  <- ifelse(pos, cp, cn)
  au  <- abs(u) + eps
  hub <- ifelse(is.infinite(cc), 1, pmin(1, cc / pmax(au, 1e-300)))
  2 * wt * hub
}

#' M-分位數損失的目標函數值
mq_obj <- function(u, W, tau = 0.5, cn = Inf, cp = Inf) {
  pos <- u > 0
  wt  <- ifelse(pos, tau, 1 - tau)
  cc  <- ifelse(pos, cp, cn)
  au  <- abs(u)
  psi <- ifelse(au <= cc, au^2 / 2, cc * au - cc^2 / 2)
  sum(2 * W * wt * psi)
}

#' M-分位數 Lee-Carter
#'
#' @param D,E     死亡數與曝露數矩陣（年齡 x 年份）
#' @param tau     非對稱參數。0.5 為對稱（Huber M-估計）
#' @param k_c     Huber 截點，以殘差的穩健尺度為單位：c = k_c * s，
#'                s 取 MAD。k_c = 0 為分位數迴歸（變體 F 的極限），
#'                k_c = Inf 為非對稱最小平方（tau=0.5 時即變體 A）
#' @param k_cn,k_cp 若給定則覆寫 k_c，允許兩側不同截點（Xu and Chen 2018）
#' @param wmode   "none" 等權重；"mu" 以配適期望死亡數加權（資訊加權）；
#'                "pow" 以 mu^wpow 加權
#' @param wpow    wmode = "pow" 時的冪次。
#'
#'   為何需要冪次。最小平方的最適權重是觀測變異數的倒數，卜瓦松下
#'   Var(log m_hat) ~ 1/mu，故 w = mu（即 wpow = 1）。但 Koenker (2005)
#'   第 5 章定理 5.1 指出，檢查函數的最適權重不是 1/sigma 而是
#'   「該分位處的局部密度」f_i(xi_i)：
#'
#'     Rather than weighting by the reciprocals of the standard deviations
#'     of the observations, quantile regression weights should be
#'     proportional to the local density evaluated at the quantile of
#'     interest.  (Koenker 2005, sec. 5.3)
#'
#'   log(D/E) 在 mu 不太小時近似 N(log m, 1/mu)，其中位數處的密度為
#'   sqrt(mu / 2pi)，故檢查函數的最適權重應為 w = sqrt(mu)，即 wpow = 0.5。
#'
#'   可檢驗的預測：最適的 wpow 應隨 c 由 0.5（c -> 0，檢查函數）移動到
#'   1.0（c -> Inf，最小平方），中間的 c 則有中間的冪次。
#' @param zero_sub 零格替代值
lc_mquantile <- function(D, E, tau = 0.5, k_c = 1.345,
                         k_cn = NULL, k_cp = NULL,
                         wmode = c("none", "mu", "pow"), wpow = 1,
                         scale_mode = c("mad", "joint"), zero_sub = 0.5,
                         maxit = 60, tol = 1e-9, init = NULL) {
  wmode <- match.arg(wmode); scale_mode <- match.arg(scale_mode)
  A <- nrow(D); Tn <- ncol(D)
  lm_ <- log(pmax(D, zero_sub) / E)
  if (is.null(k_cn)) k_cn <- k_c
  if (is.null(k_cp)) k_cp <- k_c
  ## k_c = 0 為分位數迴歸的極限。加權最小平方對權重的整體尺度不變，
  ## 故以一個極小的截點代入即得檢查函數的 IRLS 形式。
  KTINY <- 1e-6
  if (k_cn == 0) k_cn <- KTINY
  if (k_cp == 0) k_cp <- KTINY

  if (is.null(init)) {
    f0 <- lc_svd_fit(lm_); a <- f0$a; b <- f0$b; k <- f0$k
  } else { a <- init$a; b <- init$b; k <- init$k }

  Wd <- matrix(1, A, Tn)                      # 資料權重（資訊加權）
  u  <- lm_ - outer(a, rep(1, Tn)) - outer(b, k)
  ## 尺度：MAD（原作法）或由一致性方程聯合求解（Zoubir et al. 2018, sec. 3.5）
  sc_fun <- function(u, Wd, s0 = NULL) {
    if (scale_mode == "mad") max(mad(as.vector(u)), 1e-8)
    else huber_scale(u, 0.5 * (k_cn + k_cp), W = Wd, s0 = s0)
  }
  s  <- sc_fun(u, Wd)
  obj <- mq_obj(u, Wd, tau, k_cn * s, k_cp * s)

  for (it in seq_len(maxit)) {
    if (wmode != "none") {
      mu_hat <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
      pw <- if (wmode == "mu") 1 else wpow
      Wd <- if (pw == 1) mu_hat else mu_hat^pw
      Wd <- Wd / mean(Wd)
    }
    u  <- lm_ - outer(a, rep(1, Tn)) - outer(b, k)
    s  <- sc_fun(u, Wd, s)                    # 尺度
    cn <- k_cn * s; cp <- k_cp * s
    pert <- 1e-4 * s / it                     # Hunter-Lange 擾動，隨迭代縮小
      Wr <- Wd * mq_weight(u, tau, cn, cp, pert)  # 總權重 = 資料權重 x 穩健權重

    ## --- 步驟一：固定 kappa，逐年齡加權最小平方（截距與斜率）---
    for (x in seq_len(A)) {
      w <- Wr[x, ]; y <- lm_[x, ]
      sw <- sum(w); if (!is.finite(sw) || sw <= 0) next
      mk <- sum(w * k) / sw; my <- sum(w * y) / sw
      vk <- sum(w * (k - mk)^2)
      if (!is.finite(vk) || vk < 1e-12) next
      bx <- sum(w * (k - mk) * (y - my)) / vk
      b[x] <- bx; a[x] <- my - bx * mk
    }
    sb <- sum(b); if (!is.finite(sb) || abs(sb) < 1e-12) break
    b <- b / sb; k <- k * sb

    ## --- 步驟二：固定 alpha、beta，逐年份加權最小平方（無截距）---
    u  <- lm_ - outer(a, rep(1, Tn)) - outer(b, k)
    s  <- sc_fun(u, Wd, s)
    Wr <- Wd * mq_weight(u, tau, k_cn * s, k_cp * s, 1e-4 * s / it)
    Zc <- lm_ - a
    for (t in seq_len(Tn)) {
      w <- Wr[, t]
      den <- sum(w * b^2)
      if (!is.finite(den) || den < 1e-14) next
      k[t] <- sum(w * b * Zc[, t]) / den
    }
    k <- k - mean(k)
    sb <- sum(b); b <- b / sb; k <- k * sb

    u <- lm_ - outer(a, rep(1, Tn)) - outer(b, k)
    s <- sc_fun(u, Wd, s)
    obj_new <- mq_obj(u, Wd, tau, k_cn * s, k_cp * s)
    if (is.finite(obj) && abs(obj - obj_new) < tol * max(1, abs(obj))) {
      obj <- obj_new; break
    }
    obj <- obj_new
  }
  names(a) <- names(b) <- rownames(D); names(k) <- colnames(D)
  list(a = a, b = b, k = k, obj = obj, iter = it, scale = s,
       tau = tau, k_cn = k_cn, k_cp = k_cp, wmode = wmode,
       scale_mode = scale_mode,
       wpow = if (wmode == "pow") wpow else if (wmode == "mu") 1 else 0)
}
