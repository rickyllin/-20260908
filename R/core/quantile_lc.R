###############################################################################
# 分位數 Lee-Carter：以檢查函數取代平方損失
#
#   min_{alpha,beta,kappa}  sum_{x,t} w_xt * rho_tau( log m_xt - alpha_x - beta_x kappa_t )
#
#   其中 rho_tau(u) = u (tau - 1{u<0})。tau = 0.5 時即（加權）中位數迴歸。
#
#   為何值得做：檢查函數是「等權重」的——一格 16,184 死與一格 37 死在目標
#   函數裡貢獻相同的權重，與 SVD 的等權重最小平方同構。因此可預測：
#     * 中位數對「小計數取對數造成的左尾扭曲」穩健 -> 應改善 alpha
#     * 但它沒有處理加權失衡                        -> 不應改善漂移項
#   加權版（w = mu_hat）則同時具備兩者，用以檢驗「加權是否為必要條件」。
#
#   演算法：雙線性結構使目標函數對 (alpha_x, beta_x) 與對 kappa_t 分別為
#   凸的分段線性，故採交替最小化。
#     步驟一（固定 kappa）：各年齡各自是一條簡單分位數迴歸，以 quantreg::rq 解。
#     步驟二（固定 alpha, beta）：各年份各自是一維加權分位數問題，
#            tau=0.5 時解為 {z_xt / beta_x} 以 w_xt|beta_x| 為權重的加權中位數。
#   兩步皆為各自區塊的全域最小，故目標函數單調不增。
###############################################################################

invisible(lapply(c("SparseM", "quantreg"), requireNamespace, quietly = TRUE))

#' 加權分位數（tau = 0.5 時即加權中位數）
wquantile <- function(x, w, tau = 0.5) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]; w <- w[ok]
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= tau)[1]]
}

#' 目標函數值
qlc_obj <- function(lm_, a, b, k, W, tau = 0.5) {
  u <- lm_ - outer(a, rep(1, length(k))) - outer(b, k)
  sum(W * u * (tau - (u < 0)))
}

#' 分位數 Lee-Carter
#' @param wmode "none" 等權重（變體 F）；"mu" 以配適期望死亡數加權（變體 G）
#' @param init  起始值，需含 a, b, k；預設以 SVD 配適起始
lc_quantile <- function(D, E, tau = 0.5, wmode = c("none", "mu"),
                        zero_sub = 0.5, maxit = 40, tol = 1e-8, init = NULL) {
  wmode <- match.arg(wmode)
  A <- nrow(D); Tn <- ncol(D)
  lm_ <- log(pmax(D, zero_sub) / E)

  if (is.null(init)) {
    am <- rowMeans(lm_); Z <- lm_ - am
    sv <- svd(Z); u1 <- sv$u[, 1]; v1 <- sv$v[, 1]
    if (sum(u1) < 0) { u1 <- -u1; v1 <- -v1 }
    a <- am; b <- u1 / sum(u1); k <- sv$d[1] * v1 * sum(u1)
  } else { a <- init$a; b <- init$b; k <- init$k }

  W <- matrix(1, A, Tn)
  obj <- qlc_obj(lm_, a, b, k, W, tau)

  for (it in seq_len(maxit)) {
    if (wmode == "mu") {
      W <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
      W <- W / mean(W)                      # 正規化，避免目標函數尺度漂移
    }
    ## --- 步驟一：固定 kappa，各年齡一條分位數迴歸 ---
    for (x in seq_len(A)) {
      fit <- tryCatch(
        quantreg::rq(lm_[x, ] ~ k, tau = tau, weights = W[x, ],
                     method = "br"),
        error = function(e) NULL)
      if (!is.null(fit)) { a[x] <- coef(fit)[1]; b[x] <- coef(fit)[2] }
    }
    s <- sum(b); if (!is.finite(s) || abs(s) < 1e-12) break
    b <- b / s; k <- k * s

    ## --- 步驟二：固定 alpha、beta，各年份一維加權分位數 ---
    Zc <- lm_ - a
    for (t in seq_len(Tn)) {
      nz <- abs(b) > 1e-10
      k[t] <- wquantile(Zc[nz, t] / b[nz], W[nz, t] * abs(b[nz]), tau)
    }
    k <- k - mean(k)
    s <- sum(b); b <- b / s; k <- k * s

    obj_new <- qlc_obj(lm_, a, b, k, W, tau)
    if (is.finite(obj) && abs(obj - obj_new) < tol * max(1, abs(obj))) {
      obj <- obj_new; break
    }
    obj <- obj_new
  }
  names(a) <- names(b) <- rownames(D); names(k) <- colnames(D)
  list(a = a, b = b, k = k, obj = obj, iter = it, wmode = wmode)
}
