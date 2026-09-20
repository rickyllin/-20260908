###############################################################################
# Penalized Poisson Lee-Carter — L1 penalty on beta_x  (NO difference operator)
#
#   D_xt ~ Poisson(E_xt * exp(alpha_x + beta_x * kappa_t))
#
#   max  l(alpha, beta, kappa) - lambda * sum_x w_x |beta_x - beta0_x|
#   s.t. sum_x beta_x = 1,  sum_t kappa_t = 0
#
#   beta0 = 0        -> 老師指定的原始形式（見檔案末尾的說明：此時懲罰項退化）
#   beta0 = 參考曲線 -> 朝參考母體收縮，L1 真正產生選擇效果
#
# 依賴：readxl。資料檔同 lee_carter.R。
###############################################################################

library(readxl)

## ======================== 1. 讀檔（Age 欄修復同前）=========================
age_labels <- c("0","1-4","5-9","10-14",
                paste0(seq(15,105,by=5), "-", seq(19,109,by=5)), "110+")

read_mort <- function(path, labs) {
  df <- as.data.frame(read_excel(path, col_types = c("numeric","text",
                                                     "numeric","numeric","numeric")))
  names(df)[1:2] <- c("Year","Age")
  df <- df[order(df$Year), ]
  stopifnot(all(table(df$Year) == length(labs)))
  df$Age <- rep(labs, times = length(unique(df$Year)))
  df
}

to_matrix <- function(df, value, labs) {
  m <- tapply(df[[value]], list(factor(df$Age, levels = labs), df$Year), sum)
  m[labs, , drop = FALSE]
}

collapse_top <- function(M, cut = "100+") {
  top  <- switch(cut, "95+" = c("95-99","100-104","105-109","110+"),
                      "100+" = c("100-104","105-109","110+"))
  keep <- setdiff(rownames(M), top)
  rbind(M[keep, , drop = FALSE],
        matrix(colSums(M[top, , drop = FALSE]), 1, dimnames = list(cut, colnames(M))))
}

# 路徑相對於 repo 根目錄（工作目錄請設在 .Rproj 所在處）
load_data <- function(path_death = file.path("data", "五齡Death.xlsx"),
                      path_expo  = file.path("data", "5-age exposure.xls"),
                      sex = "Total", open_age = "100+") {
  d <- read_mort(path_death, age_labels)
  e <- read_mort(path_expo,  age_labels)
  Dx <- collapse_top(to_matrix(d, sex, age_labels), open_age)
  Ex <- collapse_top(to_matrix(e, sex, age_labels), open_age)
  stopifnot(all(Ex > 0))
  list(D = Dx, E = Ex, ages = rownames(Dx), years = as.numeric(colnames(Dx)))
}

## ============================ 2. 估計函數 ==================================
soft <- function(z, l) sign(z) * pmax(abs(z) - l, 0)

#' 懲罰 Poisson Lee-Carter
#' @param lambda 懲罰強度；建議以 lc_lambda_max() 正規化後用比例指定
#' @param beta0  收縮目標（長度 A）。預設 0 向量
#' @param w      各年齡的懲罰權重（長度 A）。自適應版可用 1/|beta_hat|
lc_poisson <- function(D, E, lambda = 0, beta0 = NULL, w = NULL,
                       maxit = 1000, tol = 1e-10, verbose = FALSE) {
  A <- nrow(D); Tn <- ncol(D)
  if (is.null(beta0)) beta0 <- rep(0, A)
  if (is.null(w))     w     <- rep(1, A)

  ## 起始值
  alpha <- log(rowSums(pmax(D, 0.5)) / rowSums(E))
  beta  <- rep(1/A, A)
  kappa <- seq(1, -1, length.out = Tn)

  loglik <- function(a, b, k) {
    eta <- outer(a, rep(1, Tn)) + outer(b, k)
    sum(D * eta - E * exp(eta))
  }
  pen <- function(b) lambda * sum(w * abs(b - beta0))

  obj_old <- -Inf
  for (it in seq_len(maxit)) {

    ## --- alpha：給定 beta,kappa 有封閉解（Poisson 的 offset 平移）---------
    mu    <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
    alpha <- alpha + log(rowSums(D) / rowSums(mu))

    ## --- kappa：逐年 Newton 一步，再置中 ---------------------------------
    mu    <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
    num   <- colSums((D - mu) * beta)
    den   <- colSums(mu * beta^2)
    kappa <- kappa + num / pmax(den, 1e-12)
    kappa <- kappa - mean(kappa)

    ## --- beta：proximal Newton（二次近似 + 軟閾值）+ sum(beta)=1 ----------
    mu   <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
    g    <- as.vector((D - mu) %*% kappa)          # 一階導
    Wd   <- pmax(as.vector(mu %*% kappa^2), 1e-12) # 負二階導（對角）
    btil <- beta + g / Wd                          # 未懲罰的 Newton 目標

    ## 給定乘子 nu 的封閉解；nu 用二分法使 sum(beta)=1
    bnu <- function(nu) beta0 + soft(btil - nu / Wd - beta0, lambda * w / Wd)
    lo <- -1e10; hi <- 1e10
    for (i in 1:300) {
      mid <- (lo + hi) / 2
      if (sum(bnu(mid)) > 1) lo <- mid else hi <- mid
    }
    bnew <- bnu((lo + hi) / 2)

    ## 步長折半確保目標函數不下降
    step <- 1; cur <- loglik(alpha, beta, kappa) - pen(beta)
    for (i in 1:40) {
      cand <- beta + step * (bnew - beta)
      if (abs(sum(cand)) > 1e-12) cand <- cand / sum(cand)
      if (is.finite(loglik(alpha, cand, kappa)) &&
          loglik(alpha, cand, kappa) - pen(cand) >= cur - 1e-12) break
      step <- step / 2
    }
    beta <- beta + step * (bnew - beta)

    ## --- 識別條件重整：sum(beta)=1 -----------------------------------------
    s <- sum(beta); beta <- beta / s; kappa <- kappa * s

    obj <- loglik(alpha, beta, kappa) - pen(beta)
    if (verbose && it %% 50 == 0) cat(sprintf("  it=%d obj=%.4f\n", it, obj))
    if (abs(obj - obj_old) < tol * max(1, abs(obj))) break
    obj_old <- obj
  }

  names(alpha) <- names(beta) <- rownames(D)
  names(kappa) <- colnames(D)
  mu  <- E * exp(outer(alpha, rep(1, Tn)) + outer(beta, kappa))
  dev <- 2 * sum(ifelse(D > 0, D * log(D / mu), 0) - (D - mu))

  list(alpha = alpha, beta = beta, kappa = kappa, fitted = mu,
       deviance = dev, loglik = loglik(alpha, beta, kappa),
       nfree = sum(abs(beta - beta0) > 1e-8), iter = it, lambda = lambda)
}

#' 使所有 beta_x 都被吸到 beta0 的最小 lambda（glmnet 式的 lambda_max）
lc_lambda_max <- function(D, E, beta0, w = NULL) {
  A <- nrow(D); Tn <- ncol(D)
  if (is.null(w)) w <- rep(1, A)
  f  <- lc_poisson(D, E, lambda = 0)
  mu <- f$fitted
  g  <- as.vector((D - mu) %*% f$kappa)
  Wd <- as.vector(mu %*% f$kappa^2)
  max(Wd * abs(f$beta + g / Wd - beta0) / w)
}

## ============================ 3. lambda 路徑 ===============================
lc_path <- function(D, E, beta0, w = NULL, nlam = 25, frac_min = 1e-4) {
  A <- nrow(D); Tn <- ncol(D); N <- A * Tn
  lmax <- lc_lambda_max(D, E, beta0, w)
  lams <- lmax * exp(seq(0, log(frac_min), length.out = nlam))
  res  <- lapply(lams, function(l) {
    f  <- lc_poisson(D, E, lambda = l, beta0 = beta0, w = w)
    df <- f$nfree + (Tn - 2) + A
    data.frame(lambda = l, frac = l / lmax, deviance = f$deviance,
               nfree = f$nfree, bic = f$deviance + log(N) * df,
               aic = f$deviance + 2 * df)
  })
  out <- do.call(rbind, res)
  attr(out, "lambda_max") <- lmax
  out
}

## ============================ 4. 執行範例 ==================================
if (sys.nframe() == 0) {

  dat <- load_data(sex = "Total")
  D <- dat$D; E <- dat$E; A <- nrow(D)

  ## (a) 未懲罰的 Poisson LC（基準）
  fit0 <- lc_poisson(D, E, lambda = 0)
  cat(sprintf("Poisson LC: deviance=%.1f, iter=%d, min(beta)=%.5f\n",
              fit0$deviance, fit0$iter, min(fit0$beta)))

  ## (b) 老師指定的形式：L1 朝 0 收縮 —— 對照組
  cat("\n--- L1 toward 0（對照組）---\n")
  for (l in c(0, 1e2, 1e3, 1e4, 1e5)) {
    f <- lc_poisson(D, E, lambda = l)
    cat(sprintf("lambda=%9.0f  sum|beta|=%.6f  #neg=%d  max|dbeta|=%.2e\n",
                l, sum(abs(f$beta)), sum(f$beta < 0),
                max(abs(f$beta - fit0$beta))))
  }
  cat("→ 估計值完全不變。原因見檔案末尾的恆等式。\n")

  ## (c) 朝參考曲線收縮：以 Total 為參考、估 Female
  datF <- load_data(sex = "Female")
  b_ref <- fit0$beta                       # 參考 = 全人口 Total
  pathF <- lc_path(datF$D, datF$E, beta0 = b_ref)
  print(pathF, digits = 5)

  best <- pathF[which.min(pathF$bic), ]
  cat(sprintf("\nBIC 最佳：lambda=%.4g (= %.4f * lambda_max)，自由參數 %d/%d\n",
              best$lambda, best$frac, best$nfree, A))

  fitF <- lc_poisson(datF$D, datF$E, lambda = best$lambda, beta0 = b_ref)
  snapped <- dat$ages[abs(fitF$beta - b_ref) < 1e-8]
  cat("被吸到參考曲線的年齡組：", paste(snapped, collapse = ", "), "\n")

  ## (d) 自適應權重版本：資料少的年齡罰得重
  w_ad <- 1 / sqrt(rowSums(datF$D)); w_ad <- w_ad / mean(w_ad)
  pathA <- lc_path(datF$D, datF$E, beta0 = b_ref, w = w_ad)
  cat(sprintf("自適應權重 BIC 最佳：lambda=%.4g，自由參數 %d\n",
              pathA$lambda[which.min(pathA$bic)],
              pathA$nfree[which.min(pathA$bic)]))

  ## --- 圖 ---
  op <- par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
  plot(pathF$frac, pathF$nfree, type = "b", log = "x", pch = 16,
       xlab = expression(lambda/lambda[max]), ylab = "自由 beta 個數",
       main = "稀疏度路徑")
  plot(pathF$frac, pathF$bic, type = "b", log = "x", pch = 16,
       xlab = expression(lambda/lambda[max]), ylab = "BIC", main = "BIC")
  abline(v = best$frac, col = 2, lty = 2)

  fF0 <- lc_poisson(datF$D, datF$E, lambda = 0)
  plot(seq_len(A), b_ref, type = "b", pch = 16, col = "grey50", xaxt = "n",
       ylim = range(c(b_ref, fF0$beta, fitF$beta)),
       xlab = "Age group", ylab = expression(beta[x]), main = "beta_x")
  lines(seq_len(A), fF0$beta,  type = "b", pch = 1,  col = "blue")
  lines(seq_len(A), fitF$beta, type = "b", pch = 17, col = "red")
  axis(1, at = seq_len(A), labels = dat$ages, las = 2, cex.axis = .6)
  legend("topright", c("參考(Total)", "Female 未懲罰", "Female 懲罰"),
         col = c("grey50", "blue", "red"), pch = c(16, 1, 17), bty = "n", cex = .8)
  par(op)

  write.csv(data.frame(Age = dat$ages, alpha = fitF$alpha,
                       beta_ref = b_ref, beta_unpen = fF0$beta,
                       beta_pen = fitF$beta),
            "output/tables/lc_poisson_lasso_params.csv", row.names = FALSE)
  write.csv(pathF, "output/tables/lc_poisson_lasso_path.csv", row.names = FALSE)
}

###############################################################################
# 為什麼 beta0 = 0 時懲罰項失效
#
# 在 LC 的識別條件 sum_x beta_x = 1 之下，
#
#     sum_x |beta_x| = sum_x beta_x + 2 * sum_x max(-beta_x, 0)
#                    = 1 + 2 * sum_x max(-beta_x, 0)
#
# 亦即 L1 懲罰項只罰 beta_x 的「負部」，對正值完全沒有作用。當 beta_hat 全為
# 正（大人口資料必然如此），sum|beta| 恆等於 1，懲罰項是常數，argmin 不變 ——
# 上面 (b) 的數值結果正是如此，lambda 加到 1e5 估計值一個位元都沒動。
#
# 這不代表 L1 毫無用處：在小人口下 beta_hat 會出現負值（表示該年齡死亡率隨時間
# 上升），此時 L1 的作用是「軟性非負限制」。若要讓 L1 產生年齡選擇的效果，
# 收縮目標必須改成非零的參考曲線 beta0（見 (c)），此時
# sum|beta - beta0| 在 sum(beta)=sum(beta0)=1 之下不是常數。
###############################################################################
