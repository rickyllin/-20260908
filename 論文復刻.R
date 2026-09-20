###############################################################################
# Rabbi & Mazzuco (2021) 全部分析步驟與圖表的復刻
#   資料：台灣 五齡組，預設 Female 2001-2024
#
# 需先 source("rabbi_mazzuco_replication.R") 與 source("patch_adjust_kappa.R")
#
# 對應關係
#   §3.1 Fig 1          -> fig1_smoothing()        三種修勻法比較 + 精確度表
#   §3.2 Fig 2          -> fig2_params()           a_x 與 b_x
#   §3.2 Fig 3 / Table1 -> fig3_kappa()            kappa_t 與 RW drift
#   §3.2 Fig 8          -> fig8_e0_edag()          e0 與 e0† 的關係
#   §3.3 Fig 4,5,6,7    -> fig4to7_accuracy()      樣本外準確度（面板）
#   §3.4 Table 2/Fig 9  -> fig9_e0_forecast()      2050 年 e0 與預測區間
#   §3.4 Fig 10,11      -> fig10_11_coverage()     區間預測涵蓋與誤判年數
#
# 未復刻：HU / HUR / HUW（需 demography 套件的函數型資料分解）
###############################################################################

source("rabbi_mazzuco_replication.R")
source("Patch adjust kappa.R")

SEX   <- "Female"; Y0 <- 2001; Y1 <- 2024
OUTDIR <- "rm_figs"; dir.create(OUTDIR, showWarnings = FALSE)

subset_years <- function(dat, y0, y1) {
  keep <- which(dat$years >= y0 & dat$years <= y1)
  list(D = dat$D[, keep, drop = FALSE], E = dat$E[, keep, drop = FALSE],
       ages = dat$ages, years = dat$years[keep])
}

## ===================== 0. 三種修勻法（同一組差分算子）======================
# L1 版 = smooth_l1_2d()（主檔已定義）
# L2 版：min ||y-z||^2 + lxx^2||Dxx z||^2 + lxt^2||Dxt z||^2 + ltt^2||Dtt z||^2
smooth_l2_2d <- function(logm, lxx = 5, lxt = 1, ltt = 5) {
  A <- nrow(logm); Tn <- ncol(logm); mid <- make_grid(A)$mid
  Da2 <- Matrix(D2_uneven(mid), sparse = TRUE)
  Dt2 <- Matrix(D2_uneven(seq_len(Tn) - 1), sparse = TRUE)
  Da1 <- Matrix(D1_uneven(mid), sparse = TRUE)
  Dt1 <- Matrix(D1_uneven(seq_len(Tn) - 1), sparse = TRUE)
  Rx <- kronecker(Da2, Diagonal(Tn)); Rt <- kronecker(Diagonal(A), Dt2)
  Rc <- kronecker(Da1, Dt1)
  P  <- lxx^2 * crossprod(Rx) + ltt^2 * crossprod(Rt) + lxt^2 * crossprod(Rc)
  z  <- solve(Diagonal(A * Tn) + P, as.vector(t(logm)))
  matrix(as.vector(z), A, Tn, byrow = TRUE, dimnames = dimnames(logm))
}
# 一維（僅年齡方向，逐年）：對應原文的 one-dimensional spline
smooth_l2_1d <- function(logm, lxx = 5) {
  A <- nrow(logm); mid <- make_grid(A)$mid
  Da2 <- Matrix(D2_uneven(mid), sparse = TRUE)
  M <- Diagonal(A) + lxx^2 * crossprod(Da2)
  out <- apply(logm, 2, function(y) as.vector(solve(M, y)))
  dimnames(out) <- dimnames(logm); out
}

rough_age  <- function(M) mean(abs(D2_uneven(make_grid(nrow(M))$mid) %*% M))
rough_time <- function(M) mean(abs(M %*% t(D2_uneven(seq_len(ncol(M)) - 1))))

## ===================== 1. 方法清單（含 BMS 與 MLCe†）=======================
# BMS（簡化版）：以 Poisson 離差配適 kappa（等同對齊各年齡的死亡分布），
#                並搜尋最佳配適起始年。非原文完整演算法，報告中請註明。
fit_bms <- function(D, E, min_len = 12) {
  Tn <- ncol(D)
  cand <- seq_len(Tn - min_len)
  score <- vapply(cand, function(s) {
    f <- lc_poisson(D[, s:Tn, drop = FALSE], E[, s:Tn, drop = FALSE], lambda = 0)
    f$deviance / (nrow(D) * (Tn - s + 1))
  }, numeric(1))
  s <- cand[which.min(score)]
  f <- lc_poisson(D[, s:Tn, drop = FALSE], E[, s:Tn, drop = FALSE], lambda = 0)
  list(a = f$alpha, b = f$beta, k = f$kappa, start = s)
}

fit_all <- function(D, E, lxx = 5, lxt = 1, ltt = 5, span = 25, ngrid = 120,
                    methods = c("LC","LCP","LM","BMS","LCedag","MLCedag")) {
  mobs <- D / E; logm <- log(mobs); Tn <- ncol(D)
  e0_obs <- apply(mobs, 2, function(m) life_table(m)$e0)
  ed_obs <- apply(mobs, 2, function(m) life_table(m)$edag)
  sm <- NULL
  get_sm <- function() { if (is.null(sm)) sm <<- smooth_l1_2d(logm, lxx, lxt, ltt); sm }
  
  out <- list()
  for (v in methods) {
    if (v == "LCP") {
      f <- lc_poisson(D, E, lambda = 0)
      out[[v]] <- list(a = f$alpha, b = f$beta, k = f$kappa, jumpoff = FALSE)
      
    } else if (v == "BMS") {
      f <- fit_bms(D, E)
      kk <- rep(NA_real_, Tn); kk[f$start:Tn] <- f$k
      out[[v]] <- list(a = f$a, b = f$b, k = kk, jumpoff = FALSE, start = f$start)
      
    } else if (v == "MLCedag") {
      S  <- get_sm()
      Ds <- E * exp(S)                       # 修勻後的期望死亡數
      f  <- lc_poisson(Ds, E, lambda = 0)
      kk <- adjust_kappa(f$alpha, f$beta, f$kappa,
                         function(m) life_table(m)$edag, ed_obs,
                         span = span, ngrid = ngrid, verbose = FALSE)
      out[[v]] <- list(a = f$alpha, b = f$beta, k = kk, jumpoff = TRUE)
      
    } else {
      M <- if (v == "LCedag") get_sm() else logm
      f <- lc_svd_fit(M); a <- f$a; b <- f$b; k <- f$k
      k <- switch(v,
                  LC = vapply(seq_len(Tn), function(j)
                    uniroot(function(kk) sum(E[, j] * exp(a + b * kk)) - sum(D[, j]),
                            c(k[j] - 60, k[j] + 60), extendInt = "yes",
                            tol = 1e-10)$root, numeric(1)),
                  LM = adjust_kappa(a, b, k, function(m) life_table(m)$e0, e0_obs,
                                    span = span, ngrid = ngrid, verbose = FALSE),
                  LCedag = adjust_kappa(a, b, k, function(m) life_table(m)$edag, ed_obs,
                                        span = span, ngrid = ngrid, verbose = FALSE))
      out[[v]] <- list(a = a, b = b, k = k,
                       jumpoff = v %in% c("LM", "LCedag"))
    }
    ## jump-off 修正：以最後一年實際值為起點
    o <- out[[v]]
    if (isTRUE(o$jumpoff) && !is.na(o$k[Tn]))
      out[[v]]$a <- log(mobs[, Tn]) - o$b * o$k[Tn]
    out[[v]]$smoothed <- v %in% c("LCedag", "MLCedag")
  }
  attr(out, "sm") <- sm
  out
}

drift_of <- function(k) { k <- k[!is.na(k)]; (k[length(k)] - k[1]) / (length(k) - 1) }
sigma_of <- function(k) { k <- k[!is.na(k)]; sd(diff(k) - drift_of(k)) }

## ===================== 2. Fig 1：修勻法比較 ================================
fig1_smoothing <- function(dat, lxx = 5, lxt = 1, ltt = 5, file = NULL) {
  logm <- log(dat$D / dat$E)
  S <- list(Observed = logm,
            `1D spline (L2, age)` = smooth_l2_1d(logm, lxx),
            `2D spline (L2)`      = smooth_l2_2d(logm, lxx, lxt, ltt),
            `LASSO (L1)`          = smooth_l1_2d(logm, lxx, lxt, ltt))
  mid <- make_grid(nrow(logm))$mid
  if (!is.null(file)) png(file, 1100, 950, res = 120)
  op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  cols <- rainbow(ncol(logm), start = 0, end = 0.75)
  for (nm in names(S))
    matplot(mid, S[[nm]], type = "l", lty = 1, col = cols,
            xlab = "Age", ylab = "log death rate", main = nm)
  par(op); if (!is.null(file)) dev.off()
  
  tab <- do.call(rbind, lapply(names(S)[-1], function(nm)
    data.frame(method = nm,
               MAE = mean(abs(S[[nm]] - logm)),
               MSE = mean((S[[nm]] - logm)^2),
               rough_age = rough_age(S[[nm]]),
               rough_time = rough_time(S[[nm]]))))
  tab <- rbind(data.frame(method = "Observed", MAE = 0, MSE = 0,
                          rough_age = rough_age(logm),
                          rough_time = rough_time(logm)), tab)
  list(surfaces = S, table = tab)
}

## ===================== 3. Fig 2 / Fig 3 / Table 1 ==========================
fig2_params <- function(dat, fits, file = NULL) {
  A <- nrow(dat$D); mid <- make_grid(A)$mid
  cols <- setNames(c("black","grey50","blue","darkgreen","red","orange"),
                   c("LC","LCP","LM","BMS","LCedag","MLCedag"))
  if (!is.null(file)) png(file, 1100, 500, res = 120)
  op <- par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(mid, fits$LC$a, type = "n", xlab = "Age", ylab = expression(a[x]),
       ylim = safe_range(sapply(fits, `[[`, "a")), main = expression(a[x]))
  for (v in names(fits)) lines(mid, fits[[v]]$a, col = cols[v], lwd = 2)
  plot(mid, fits$LC$b, type = "n", xlab = "Age", ylab = expression(b[x]),
       ylim = safe_range(sapply(fits, `[[`, "b")), main = expression(b[x]))
  for (v in names(fits)) lines(mid, fits[[v]]$b, col = cols[v], lwd = 2)
  abline(h = 0, lty = 3)
  legend("topright", names(fits), col = cols[names(fits)], lwd = 2,
         bty = "n", cex = .75)
  par(op); if (!is.null(file)) dev.off()
}

fig3_kappa <- function(dat, fits, h = 26, file = NULL) {
  yrs <- dat$years; last <- max(yrs)
  if (!is.null(file)) png(file, 1100, 800, res = 120)
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  for (v in names(fits)) {
    k <- fits[[v]]$k; ok <- !is.na(k)
    dr <- drift_of(k); sg <- sigma_of(k)
    kf <- k[max(which(ok))] + dr * seq_len(h)
    se <- sg * sqrt(seq_len(h)); yf <- last + seq_len(h)
    plot(c(yrs, yf), c(k, kf), type = "n", xlab = "Year",
         ylab = expression(kappa[t]), main = v)
    polygon(c(yf, rev(yf)), c(kf - 1.96*se, rev(kf + 1.96*se)),
            col = rgb(0,0,1,.15), border = NA)
    lines(yrs[ok], k[ok], lwd = 2); lines(yf, kf, col = 4, lwd = 2)
  }
  par(op); if (!is.null(file)) dev.off()
  
  data.frame(method = names(fits),
             drift = sapply(fits, function(f) drift_of(f$k)),
             sigma = sapply(fits, function(f) sigma_of(f$k)),
             k_lin_resid = sapply(fits, function(f) {
               k <- f$k[!is.na(f$k)]; sd(residuals(lm(k ~ seq_along(k)))) }),
             row.names = NULL)
}

## ===================== 4. Fig 8：e0 與 e0† 的關係 ==========================
fitted_rates <- function(f, Tn) {
  ok <- !is.na(f$k)
  exp(outer(f$a, rep(1, sum(ok))) + outer(f$b, f$k[ok]))
}

fig8_e0_edag <- function(dat, fits, which = c("LC","LCedag"), file = NULL) {
  mobs <- dat$D / dat$E
  obs <- data.frame(e0 = apply(mobs, 2, function(m) life_table(m)$e0),
                    ed = apply(mobs, 2, function(m) life_table(m)$edag))
  pans <- c(list(Observed = obs),
            setNames(lapply(which, function(v) {
              M <- fitted_rates(fits[[v]], ncol(mobs))
              data.frame(e0 = apply(M, 2, function(m) life_table(m)$e0),
                         ed = apply(M, 2, function(m) life_table(m)$edag))
            }), which))
  if (!is.null(file)) png(file, 1100, 420, res = 120)
  op <- par(mfrow = c(1, length(pans)), mar = c(4, 4, 3, 1))
  xl <- safe_range(sapply(pans, `[[`, "e0")); yl <- safe_range(sapply(pans, `[[`, "ed"))
  for (nm in names(pans)) {
    p <- pans[[nm]]
    plot(p$e0, p$ed, pch = 16, xlim = xl, ylim = yl,
         col = rainbow(nrow(p), start = 0, end = .75),
         xlab = expression(e[0]), ylab = expression(e[0]^"†"),
         main = sprintf("%s  (r = %.3f)", nm, cor(p$e0, p$ed)))
  }
  par(op); if (!is.null(file)) dev.off()
  data.frame(panel = names(pans),
             r = sapply(pans, function(p) cor(p$e0, p$ed)), row.names = NULL)
}

## ===================== 5. Fig 4-7：樣本外準確度（面板）=====================
oos_one <- function(D, E, h, methods, ...) {
  Tn <- ncol(D); fT <- Tn - h
  Df <- D[, 1:fT, drop = FALSE]; Ef <- E[, 1:fT, drop = FALSE]
  Dh <- D[, (fT+1):Tn, drop = FALSE]; Eh <- E[, (fT+1):Tn, drop = FALSE]
  lmh <- log(Dh / Eh)
  e0h <- apply(Dh/Eh, 2, function(m) life_table(m)$e0)
  edh <- apply(Dh/Eh, 2, function(m) life_table(m)$edag)
  fits <- fit_all(Df, Ef, methods = methods, ...)
  do.call(rbind, lapply(names(fits), function(v) {
    f <- fits[[v]]; ok <- !is.na(f$k)
    if (sum(ok) < 5) return(NULL)
    dr <- drift_of(f$k); kf <- f$k[max(which(ok))] + dr * seq_len(h)
    pred <- outer(f$a, rep(1, h)) + outer(f$b, kf)
    e0p <- apply(exp(pred), 2, function(m) life_table(m)$e0)
    edp <- apply(exp(pred), 2, function(m) life_table(m)$edag)
    err <- pred - lmh
    data.frame(method = v, MAE = mean(abs(err)), MSE = mean(err^2),
               MAE_h1 = mean(abs(err[,1])), MAE_hH = mean(abs(err[,h])),
               ME_e0 = mean(e0h - e0p), MAE_e0 = mean(abs(e0h - e0p)),
               MAE_edag = mean(abs(edh - edp)), drift = dr)
  }))
}

# 以「性別 x 配適起始年」取代原文的 20 個國家
fig4to7_accuracy <- function(h = 5, sexes = c("Female","Male","Total"),
                             starts = c(2001, 1991, 1981, 1971), y1 = 2024,
                             methods = c("LC","LCP","LM","BMS","LCedag","MLCedag"),
                             file = NULL) {
  res <- list()
  for (s in sexes) {
    base <- load_data(sex = s)
    for (y0 in starts) {
      dd <- subset_years(base, y0, y1)
      r <- oos_one(dd$D, dd$E, h, methods)
      if (!is.null(r)) { r$sex <- s; r$start <- y0; res[[length(res)+1]] <- r }
      cat(sprintf("  %s %d-%d 完成\n", s, y0, y1))
    }
  }
  out <- do.call(rbind, res)
  out$panel <- paste0(substr(out$sex,1,1), out$start)
  
  if (!is.null(file)) png(file, 1200, 950, res = 120)
  op <- par(mfrow = c(2, 2), mar = c(5, 4, 3, 1))
  for (m in c("MAE","MSE","ME_e0","MAE_edag")) {
    mm <- tapply(out[[m]], list(out$panel, out$method), identity)
    matplot(seq_len(ncol(mm)), t(mm), type = "p", pch = 16,
            col = rainbow(nrow(mm)), xaxt = "n", xlab = "", ylab = m,
            main = sprintf("%s (out-of-sample, h=%d)", m, h))
    axis(1, at = seq_len(ncol(mm)), labels = colnames(mm), las = 2, cex.axis = .8)
    if (m == "ME_e0") abline(h = 0, lty = 2, col = "grey40")
  }
  par(op); if (!is.null(file)) dev.off()
  out
}

## ===================== 6. Table 2 / Fig 9：e0 預測與區間 ===================
# 半參數 bootstrap：對 kappa 的 RW with drift 加擾動，乘上固定 b_x
sim_e0 <- function(f, h, nsim = 1000, seed = 1) {
  set.seed(seed)
  ok <- !is.na(f$k); dr <- drift_of(f$k); sg <- sigma_of(f$k)
  k0 <- f$k[max(which(ok))]
  paths <- k0 + t(apply(matrix(rnorm(nsim*h, dr, sg), nsim, h), 1, cumsum))
  apply(paths, 2, function(col)
    vapply(col, function(x) life_table(exp(f$a + f$b * x))$e0, numeric(1)))
}

fig9_e0_forecast <- function(dat, fits, to = 2050, nsim = 1000, file = NULL) {
  last <- max(dat$years); h <- to - last; yf <- last + seq_len(h)
  mobs <- dat$D / dat$E
  e0obs <- apply(mobs, 2, function(m) life_table(m)$e0)
  
  sims <- lapply(fits, sim_e0, h = h, nsim = nsim)
  if (!is.null(file)) png(file, 1200, 800, res = 120)
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
  for (v in names(fits)) {
    S <- sims[[v]]
    q <- apply(S, 2, quantile, c(.025,.1,.5,.9,.975))
    plot(c(dat$years, yf), c(e0obs, q[3,]), type = "n",
         ylim = safe_range(e0obs, q), xlab = "Year", ylab = expression(e[0]),
         main = v)
    polygon(c(yf, rev(yf)), c(q[1,], rev(q[5,])), col = rgb(1,0,0,.12), border = NA)
    polygon(c(yf, rev(yf)), c(q[2,], rev(q[4,])), col = rgb(0,0,1,.2),  border = NA)
    lines(dat$years, e0obs, lwd = 2); lines(yf, q[3,], col = 4, lwd = 2)
  }
  par(op); if (!is.null(file)) dev.off()
  
  data.frame(method = names(fits),
             e0_2050 = sapply(sims, function(S) median(S[, ncol(S)])),
             lo80 = sapply(sims, function(S) quantile(S[, ncol(S)], .1)),
             hi80 = sapply(sims, function(S) quantile(S[, ncol(S)], .9)),
             width80 = sapply(sims, function(S)
               diff(quantile(S[, ncol(S)], c(.1,.9)))), row.names = NULL)
}

## ===================== 7. Fig 10/11：區間預測的涵蓋 ========================
fig10_11_coverage <- function(D, E, years, h = 5, nsim = 2000,
                              methods = c("LC","LCP","LM","BMS","LCedag","MLCedag"),
                              file = NULL) {
  Tn <- ncol(D); fT <- Tn - h
  Dh <- D[, (fT+1):Tn, drop = FALSE]; Eh <- E[, (fT+1):Tn, drop = FALSE]
  e0h <- apply(Dh/Eh, 2, function(m) life_table(m)$e0)
  fits <- fit_all(D[,1:fT,drop=FALSE], E[,1:fT,drop=FALSE], methods = methods)
  
  res <- do.call(rbind, lapply(names(fits), function(v) {
    S <- sim_e0(fits[[v]], h, nsim)
    q <- apply(S, 2, quantile, c(.1,.9))
    inside <- e0h >= q[1,] & e0h <= q[2,]
    data.frame(method = v, covered = sum(inside), missed = h - sum(inside),
               emp_cov = mean(inside), CPD = abs(0.8 - mean(inside)),
               mean_width = mean(q[2,] - q[1,]))
  }))
  
  if (!is.null(file)) png(file, 1200, 500, res = 120)
  op <- par(mfrow = c(1, 2), mar = c(6, 4, 3, 1))
  barplot(res$missed, names.arg = res$method, las = 2, ylab = "誤判年數",
          main = sprintf("Interval forecast of e0 (h=%d)", h))
  barplot(res$CPD, names.arg = res$method, las = 2,
          ylab = "|0.8 - empirical coverage|", main = "Coverage probability deviation")
  par(op); if (!is.null(file)) dev.off()
  res
}

## ===================== 8. 執行全部 =========================================
if (sys.nframe() == 0) {
  
  dat <- subset_years(load_data(sex = SEX), Y0, Y1)
  cat(sprintf("== %s %d-%d：%d 齡組 x %d 年 ==\n", SEX, Y0, Y1,
              nrow(dat$D), ncol(dat$D)))
  
  cat("\n[Fig 1] 修勻法比較\n")
  f1 <- fig1_smoothing(dat, file = file.path(OUTDIR, "fig1_smoothing.png"))
  print(f1$table, digits = 4, row.names = FALSE)
  
  cat("\n[Fig 2/3, Table 1] 參數與 kappa\n")
  fits <- fit_all(dat$D, dat$E)
  fig2_params(dat, fits, file = file.path(OUTDIR, "fig2_params.png"))
  t1 <- fig3_kappa(dat, fits, file = file.path(OUTDIR, "fig3_kappa.png"))
  print(t1, digits = 4, row.names = FALSE)
  
  cat("\n[Fig 8] e0 與 e0† 的關係\n")
  print(fig8_e0_edag(dat, fits, file = file.path(OUTDIR, "fig8_e0_edag.png")),
        digits = 4, row.names = FALSE)
  
  cat("\n[Fig 4-7] 樣本外準確度（性別 x 起始年 面板；需數分鐘）\n")
  acc <- fig4to7_accuracy(h = 5, file = file.path(OUTDIR, "fig4to7_accuracy.png"))
  print(aggregate(cbind(MAE, MSE, ME_e0, MAE_edag) ~ method, acc, mean),
        digits = 4, row.names = FALSE)
  write.csv(acc, file.path(OUTDIR, "accuracy_panel.csv"), row.names = FALSE)
  
  cat("\n[Table 2 / Fig 9] 2050 年 e0 預測\n")
  t2 <- fig9_e0_forecast(dat, fits, file = file.path(OUTDIR, "fig9_e0_forecast.png"))
  print(t2, digits = 5, row.names = FALSE)
  
  cat("\n[Fig 10/11] 區間預測涵蓋\n")
  cov <- fig10_11_coverage(dat$D, dat$E, dat$years, h = 5,
                           file = file.path(OUTDIR, "fig10_coverage.png"))
  print(cov, digits = 4, row.names = FALSE)
  
  write.csv(f1$table, file.path(OUTDIR, "table_smoothing.csv"), row.names = FALSE)
  write.csv(t1, file.path(OUTDIR, "table1_drift.csv"), row.names = FALSE)
  write.csv(t2, file.path(OUTDIR, "table2_e0_2050.csv"), row.names = FALSE)
  write.csv(cov, file.path(OUTDIR, "table_coverage.csv"), row.names = FALSE)
  cat(sprintf("\n圖表已輸出到 %s/\n", OUTDIR))
}

###############################################################################
# 給老師報告時要一併說明的三點
#
# 1. Fig 4-7、10 的面板維度是「性別 x 配適起始年」，不是原文的 20 個國家。
#    原文用跨國變異來論證穩健性，本復刻改用跨子群與跨窗口變異。
#
# 2. 原文比較的一維／二維樣條，本檔以同一組差分算子的 L2 版本代替，
#    而非 MortalitySmooth / 約束迴歸樣條。好處是三種修勻法只差在範數，
#    是乾淨的受控比較；代價是與原文的數值不能逐項對照。
#
# 3. HU / HUR / HUW 未復刻（需 demography 的函數型多主成分分解）。
#
# 另：留出期 2020-2024 涵蓋 COVID，女性 e0 在 2020 達 84.30 後於 2022 降至
#     83.14，2024 才回到 83.94。所有方法的 ME(e0) 因此一致為負（高估存活
#     改善），這是真實的結構性衝擊而非模型缺陷。若要乾淨的比較，另外跑
#     Y1 = 2019 的版本當對照。
###############################################################################