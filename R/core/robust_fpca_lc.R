###############################################################################
# 對「離群年份」穩健的 Lee-Carter：Hyndman and Ullah (2007) 取徑的核心步驟
#
#   Hyndman and Ullah (2007) 在函數型資料架構下提出穩健版的 LC，是死亡率
#   文獻中最直接的穩健化嘗試。其穩健化有兩個步驟：
#     (i)  先以懲罰式樣條把各年的 log m_{x,t} 在年齡方向修勻；
#     (ii) 再以穩健主成分分析取代一般 SVD，作法是依「每一年整條曲線」的
#          積分平方誤差判斷該年是否為離群年，並降低其權重。
#   其設想的離群值是戰爭、流感大流行、地震等使「整個年度」偏離趨勢的事件，
#   Hyndman and Ullah (2007) 所舉的例子即法國 1914-1918 與 1943-1945。
#
#   本檔實作步驟 (ii)。步驟 (i) 屬於年齡方向的修勻，本文第肆節之一已另行
#   檢驗，且本研究採五齡組（22 組）而非單齡組，修勻空間本就有限。
#
#   為何值得單獨檢驗。第伍節之五已顯示逐格截尾（LTS/MTL）失效，但那是
#   「格」為單位；Hyndman and Ullah (2007) 是以「年」為單位，兩者並不相同：
#   若小人口的噪音在某些年份剛好同向累積，逐年穩健化仍可能有效。
#   反之，若噪音在每一格獨立且每一年同樣嘈雜，則不存在離群年，
#   逐年穩健化應毫無作用——這是一項可事先寫下、方向明確的預測。
#
#   逐年權重的計算（對應 Hyndman-Ullah 的積分平方誤差判準）：
#       v_t = sum_x ( log m_{x,t} - alpha_x - beta_x kappa_t )^2
#   以 v_t 的中位數與 MAD 作穩健定位與尺度，
#       rule = "huber"：w_t = min(1, lambda * MAD / (v_t - med)_+)
#       rule = "trim" ：把 v_t 最大的 n_trim 年權重設為 0
#   再以逐年權重重解秩一結構，反覆至收斂。
###############################################################################

#' 對離群年份穩健的 LC（Hyndman-Ullah 取徑的穩健 PCA 步驟）
#'
#' @param rule   "huber" 平滑降權；"trim" 硬性剔除
#' @param lambda rule = "huber" 時的截點倍數
#' @param n_trim rule = "trim" 時剔除的年數
lc_robust_year <- function(D, E, rule = c("huber", "trim"),
                           lambda = 3, n_trim = 2, zero_sub = 0.5,
                           maxit = 50, tol = 1e-10) {
  rule <- match.arg(rule)
  A <- nrow(D); Tn <- ncol(D)
  lm_ <- log(pmax(D, zero_sub) / E)

  f0 <- lc_svd_fit(lm_); a <- f0$a; b <- f0$b; k <- f0$k
  wt <- rep(1, Tn); obj <- Inf

  for (it in seq_len(maxit)) {
    ## --- 逐年權重：以整年曲線的積分平方誤差判斷離群年 ---
    R  <- lm_ - outer(a, rep(1, Tn)) - outer(b, k)
    v  <- colSums(R^2)
    md <- median(v); sc <- max(mad(v), 1e-12)
    if (rule == "huber") {
      ex <- pmax(v - md, 0)
      wt <- pmin(1, lambda * sc / pmax(ex, 1e-12))
    } else {
      wt <- rep(1, Tn)
      if (n_trim > 0) wt[order(v, decreasing = TRUE)[seq_len(n_trim)]] <- 0
    }

    ## --- 以逐年權重重解秩一結構 ---
    sw <- sum(wt); if (sw <= 0) break
    a  <- as.vector(lm_ %*% wt) / sw                 # 加權的逐年齡平均
    Z  <- (lm_ - a) * rep(sqrt(wt), each = A)        # 權重吸收進行（年）方向
    sv <- svd(Z); u1 <- sv$u[, 1]; v1 <- sv$v[, 1]
    if (sum(u1) < 0) { u1 <- -u1; v1 <- -v1 }
    b  <- u1 / sum(u1)
    kk <- sv$d[1] * v1 * sum(u1)
    k  <- ifelse(wt > 0, kk / sqrt(pmax(wt, 1e-12)), NA_real_)
    ## 被剔除的年份仍需一個 kappa 值（否則無法算漂移項）：
    ## 以該年的加權最小平方投影回補，這是 Hyndman-Ullah 的作法
    for (t in which(!is.finite(k))) k[t] <- sum(b * (lm_[, t] - a)) / sum(b^2)
    k  <- k - mean(k)
    sb <- sum(b); b <- b / sb; k <- k * sb

    R  <- lm_ - outer(a, rep(1, Tn)) - outer(b, k)
    obj_new <- sum(rep(wt, each = A) * R^2)
    if (is.finite(obj) && abs(obj - obj_new) < tol * max(1, abs(obj))) {
      obj <- obj_new; break
    }
    obj <- obj_new
  }
  names(a) <- names(b) <- rownames(D); names(k) <- colnames(D)
  list(a = a, b = b, k = k, obj = obj, iter = it,
       year_weight = wt, n_down = sum(wt < 0.999), rule = rule)
}
