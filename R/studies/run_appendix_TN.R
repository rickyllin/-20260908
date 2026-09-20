###############################################################################
# (A) 附錄：懲罰前後 beta_x 的逐年齡比較
# (B) 工程二：配適期長度 T x 人口規模 N 的二因子設計
# (C) Elastic Net 變體的比較（adaptive / relaxed / MCP）
#
# 輸出至 pen_figs/
###############################################################################
source("R/core/penalized_lc.R"); source("R/studies/mc_exposure.R")
FIGDIR <- "output/figures"; TABDIR <- "output/tables"
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
REPS <- 100; SEED <- 20260921

truth <- mc_truth("Female", 2001, 2024)
A <- length(truth$b); B0 <- rep(1/A, A); ages <- truth$ages

###############################################################################
# (A) beta_x 的逐年齡比較
###############################################################################
cat("[A] beta_x 逐年齡比較\n")

## --- A1. 實際資料（女性 2001-2024，單一性別約 1,190 萬人）---
dat <- load_data(sex = "Female"); kk <- which(dat$years >= 2001)
Dr <- dat$D[, kk]; Er <- dat$E[, kk]
lmax_r <- lc_pen_lambda_max(Dr, Er, B0)
fr <- list(
  unpen = lc_penalized(Dr, Er, 0),
  lasso = lc_penalized(Dr, Er, 0.10 * lmax_r, 1,   B0),
  enet  = lc_penalized(Dr, Er, 0.10 * lmax_r, 0.5, B0),
  ridge = lc_penalized(Dr, Er, 0.003 * lc_pen_lambda_max(Dr, Er, B0, NULL, 0), 0, B0))
tabA1 <- data.frame(
  age = ages, beta0 = B0,
  beta_unpen = fr$unpen$beta, beta_lasso = fr$lasso$beta,
  beta_enet  = fr$enet$beta,  beta_ridge = fr$ridge$beta)
tabA1$diff_lasso <- tabA1$beta_lasso - tabA1$beta_unpen
tabA1$diff_ridge <- tabA1$beta_ridge - tabA1$beta_unpen
tabA1$snapped    <- ifelse(abs(tabA1$beta_lasso - B0) < 1e-7, "是", "")
print(tabA1, digits = 4, row.names = FALSE)
write.csv(tabA1, file.path(TABDIR, "tableA1_beta_by_age_real.csv"), row.names = FALSE)

## --- A2. 模擬（有真值，可分解逐年齡的偏誤與標準差）---
set.seed(SEED)
simA <- function(N, frac) {
  E <- N * truth$w
  lam_all <- numeric(REPS)
  BU <- BP <- matrix(NA_real_, REPS, A)
  for (r in seq_len(REPS)) {
    D <- matrix(rpois(length(E), E * truth$m), nrow(E), ncol(E),
                dimnames = dimnames(truth$m))
    lm_ <- tryCatch(lc_pen_lambda_max(D, E, B0), error = function(e) NA)
    if (!is.finite(lm_)) next
    lam_all[r] <- lm_
    BU[r, ] <- tryCatch(lc_penalized(D, E, 0)$beta, error=function(e) rep(NA,A))
    BP[r, ] <- tryCatch(lc_penalized(D, E, frac*lm_, 1, B0)$beta,
                        error=function(e) rep(NA,A))
  }
  ok <- complete.cases(BU) & complete.cases(BP)
  BU <- BU[ok, , drop=FALSE]; BP <- BP[ok, , drop=FALSE]
  data.frame(N = N, age = ages, beta_true = truth$b,
             mean_unpen = colMeans(BU), sd_unpen = apply(BU, 2, sd),
             bias_unpen = colMeans(BU) - truth$b,
             mean_pen = colMeans(BP), sd_pen = apply(BP, 2, sd),
             bias_pen = colMeans(BP) - truth$b,
             snap_rate = colMeans(abs(BP - matrix(B0, nrow(BP), A, byrow=TRUE)) < 1e-7))
}
a2 <- rbind(simA(5e4, 0.20), simA(2e5, 0.03))
write.csv(a2, file.path(TABDIR, "tableA2_beta_by_age_sim.csv"), row.names = FALSE)
cat("\nN = 5e4（最適 frac = 0.20）：\n")
print(a2[a2$N==5e4, c("age","beta_true","mean_unpen","sd_unpen","mean_pen","sd_pen","snap_rate")],
      digits = 3, row.names = FALSE)

## --- 圖 ---
mid <- make_grid(A)$mid
png(file.path(FIGDIR, "figA_beta_by_age.png"), 1400, 900, res = 120)
op <- par(mfrow = c(2, 2), mar = c(4, 4.4, 3, 1))
# (1) 實際資料
matplot(mid, cbind(tabA1$beta_unpen, tabA1$beta_lasso, tabA1$beta_enet, tabA1$beta_ridge),
        type = "b", pch = c(16,17,15,18), lty = 1, col = c("black","red","blue","darkgreen"),
        xlab = "Age", ylab = expression(beta[x]),
        main = "Real data (Female 2001-2024)")
abline(h = 1/A, lty = 2, col = "grey40")
legend("topright", c("unpenalized","LASSO","Elastic Net","Ridge","1/A"),
       col = c("black","red","blue","darkgreen","grey40"),
       pch = c(16,17,15,18,NA), lty = c(1,1,1,1,2), bty = "n", cex = .7)
# (2)(3) 模擬：平均估計值
for (N in c(5e4, 2e5)) {
  s <- a2[a2$N == N, ]
  matplot(mid, cbind(s$beta_true, s$mean_unpen, s$mean_pen), type = "b",
          pch = c(NA,16,17), lty = c(1,1,1), lwd = c(2,1,1),
          col = c("grey30","black","red"), xlab = "Age", ylab = expression(beta[x]),
          main = sprintf("Simulation, N = %.0e", N))
  abline(h = 1/A, lty = 2, col = "grey60")
  legend("topright", c("truth","unpenalized (mean)","LASSO (mean)"),
         col = c("grey30","black","red"), lty = 1, pch = c(NA,16,17), bty="n", cex=.7)
}
# (4) 逐年齡的標準差比較
s <- a2[a2$N == 5e4, ]
plot(mid, s$sd_unpen, type = "b", pch = 16, log = "y",
     ylim = range(c(s$sd_unpen, pmax(s$sd_pen, 1e-6))),
     xlab = "Age", ylab = "SD of beta_x", main = "Per-age SD (N = 5e4)")
lines(mid, pmax(s$sd_pen, 1e-6), type = "b", pch = 17, col = "red")
legend("bottomleft", c("unpenalized","LASSO"), col = c("black","red"),
       pch = c(16,17), bty = "n", cex = .75)
par(op); dev.off()

###############################################################################
# (B) T x N 二因子設計
###############################################################################
cat("\n[B] T x N 二因子設計\n")
truthL <- mc_truth("Female", 1970, 2024)          # 55 年的長真值
AL <- length(truthL$b); B0L <- rep(1/AL, AL); TT <- length(truthL$k)

tn_cell <- function(T_, N, reps = REPS, frac = 0.1) {
  idx <- (TT - T_ + 1):TT                          # 取最近的 T 年
  E <- N * truthL$w[, idx, drop = FALSE]
  m <- truthL$m[, idx, drop = FALSE]
  kt <- truthL$k[idx]
  drift_true <- (kt[T_] - kt[1]) / (T_ - 1)
  su <- sp <- numeric(0); du <- dp <- numeric(0)
  for (r in seq_len(reps)) {
    D <- matrix(rpois(length(E), E * m), nrow(E), ncol(E), dimnames = dimnames(m))
    fu <- tryCatch(lc_penalized(D, E, 0), error = function(e) NULL)
    lm_ <- tryCatch(lc_pen_lambda_max(D, E, B0L), error = function(e) NA)
    fp <- if (is.finite(lm_))
      tryCatch(lc_penalized(D, E, frac*lm_, 1, B0L), error=function(e) NULL) else NULL
    if (!is.null(fu)) { su <- c(su, sum((fu$beta - truthL$b)^2))
      du <- c(du, (fu$kappa[T_] - fu$kappa[1])/(T_-1) - drift_true) }
    if (!is.null(fp)) { sp <- c(sp, sum((fp$beta - truthL$b)^2))
      dp <- c(dp, (fp$kappa[T_] - fp$kappa[1])/(T_-1) - drift_true) }
  }
  data.frame(T = T_, N = N, drift_true = drift_true, nok = length(su),
             sse_unpen = median(su), sse_pen = median(sp),
             drift_bias_unpen = median(du), drift_bias_pen = median(dp))
}
set.seed(SEED)
grid <- expand.grid(T_ = c(10, 15, 20, 30, 40, 55), N = c(2e4, 5e4, 2e5, 1e6))
tn <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  cat(sprintf("  T=%2d N=%.0e\n", grid$T_[i], grid$N[i]))
  tn_cell(grid$T_[i], grid$N[i])
}))
write.csv(tn, file.path(TABDIR, "tableB_TxN.csv"), row.names = FALSE)
cat("\nSSE(beta) 中位數（未懲罰）：\n")
print(reshape(tn[, c("T","N","sse_unpen")], idvar="T", timevar="N", direction="wide"),
      digits = 3, row.names = FALSE)

png(file.path(FIGDIR, "figB_TxN.png"), 1300, 500, res = 120)
op <- par(mfrow = c(1, 3), mar = c(4, 4.4, 3, 1))
Ns <- sort(unique(tn$N)); cols <- c("black","blue","darkgreen","red")
plot(NA, xlim=range(tn$T), ylim=range(tn$sse_unpen), log="y",
     xlab="Fitting window T (years)", ylab="median SSE(beta)",
     main="Unpenalized")
for (i in seq_along(Ns)) { s<-tn[tn$N==Ns[i],]; lines(s$T, s$sse_unpen, type="b", pch=16, col=cols[i]) }
legend("topright", sprintf("N=%.0e", Ns), col=cols, pch=16, lwd=1, bty="n", cex=.7)
plot(NA, xlim=range(tn$T), ylim=range(tn$sse_pen), log="y",
     xlab="Fitting window T (years)", ylab="median SSE(beta)",
     main="LASSO (frac=0.1)")
for (i in seq_along(Ns)) { s<-tn[tn$N==Ns[i],]; lines(s$T, s$sse_pen, type="b", pch=17, col=cols[i]) }
plot(NA, xlim=range(tn$T), ylim=range(tn$drift_bias_unpen),
     xlab="Fitting window T (years)", ylab="median drift bias", main="Drift bias (unpenalized)")
for (i in seq_along(Ns)) { s<-tn[tn$N==Ns[i],]; lines(s$T, s$drift_bias_unpen, type="b", pch=16, col=cols[i]) }
abline(h=0, lty=2)
par(op); dev.off()

###############################################################################
# (C) Elastic Net 變體
###############################################################################
cat("\n[C] Elastic Net 變體\n")
set.seed(SEED)
variants <- function(N, fracs = c(0.003,0.01,0.03,0.1,0.2,0.4)) {
  E <- N * truth$w
  acc <- array(NA_real_, c(REPS, length(fracs), 5, A))  # 5 個方法
  for (r in seq_len(REPS)) {
    D <- matrix(rpois(length(E), E*truth$m), nrow(E), ncol(E), dimnames=dimnames(truth$m))
    lm_ <- tryCatch(lc_pen_lambda_max(D,E,B0), error=function(e) NA)
    if (!is.finite(lm_)) next
    wad <- tryCatch(adaptive_w(D,E,B0), error=function(e) rep(1,A))
    for (j in seq_along(fracs)) {
      l <- fracs[j]*lm_
      g <- function(f) tryCatch(f$beta, error=function(e) rep(NA,A))
      acc[r,j,1,] <- g(tryCatch(lc_penalized(D,E,l,1,B0), error=function(e) NULL))
      acc[r,j,2,] <- g(tryCatch(lc_penalized(D,E,l,0.5,B0), error=function(e) NULL))
      acc[r,j,3,] <- g(tryCatch(lc_penalized(D,E,l,1,B0,w=wad), error=function(e) NULL))
      acc[r,j,4,] <- g(tryCatch(lc_relaxed(D,E,l,1,B0), error=function(e) NULL))
      acc[r,j,5,] <- g(tryCatch(lc_mcp(D,E,l,3,B0), error=function(e) NULL))
    }
    if (r %% 25 == 0) cat(sprintf("   rep %d/%d\n", r, REPS))
  }
  nm <- c("LASSO","Elastic Net","adaptive LASSO","relaxed LASSO","MCP")
  do.call(rbind, lapply(1:5, function(k) do.call(rbind, lapply(seq_along(fracs), function(j) {
    B <- acc[,j,k,]; B <- B[complete.cases(B), , drop=FALSE]
    if (!nrow(B)) return(NULL)
    data.frame(method=nm[k], frac=fracs[j], nrep=nrow(B),
               bias2=sum((colMeans(B)-truth$b)^2),
               variance=sum(apply(B,2,var)),
               mse=mean(rowSums(sweep(B,2,truth$b)^2)))
  }))))
}
vc <- variants(2e5)
write.csv(vc, file.path(TABDIR, "tableC_variants_grid.csv"), row.names = FALSE)
best <- do.call(rbind, lapply(split(vc, vc$method), function(s) s[which.min(s$mse),]))
best <- best[order(best$mse), ]
cat("\n各變體的最佳表現（N = 2e5）：\n"); print(best, digits=4, row.names=FALSE)
write.csv(best, file.path(TABDIR, "tableC_variants_best.csv"), row.names = FALSE)
cat("\n完成\n")
