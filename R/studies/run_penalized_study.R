###############################################################################
# 懲罰式 Lee-Carter 的實測：收斂、偏誤—變異數權衡、變異數收斂速度
#
# 對應老師指示的三項工作：
#   (1) 加入 beta_x 的 L1 懲罰後與直接估計比較
#   (2) 評估 penalized 造成的 bias 與 variance 縮減之間的 trade-off
#   (3) 評估不同人口暴露數之下變異數的收斂速度
#
# 輸出至 pen_figs/
###############################################################################

source("R/core/penalized_lc.R")
source("R/studies/mc_exposure.R")

FIGDIR <- "output/figures"; TABDIR <- "output/tables"
dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
REPS <- 100
truth <- mc_truth("Female", 2001, 2024)
A <- length(truth$b); B0 <- rep(1 / A, A)

cat(sprintf("真值：%d 齡組 x %d 年；beta 範圍 [%.4f, %.4f]\n",
            A, length(truth$k), min(truth$b), max(truth$b)))

## ========== 0. 懲罰項退化的示範（beta0 = 0 vs 1/A）========================
cat("\n[0] 收縮目標的選擇\n")
dat <- load_data(sex = "Female"); kk <- which(dat$years >= 2001)
Dr <- dat$D[, kk]; Er <- dat$E[, kk]
f_un <- lc_penalized(Dr, Er, lambda = 0)
deg <- do.call(rbind, lapply(c(0, 1e2, 1e4, 1e6), function(l) {
  fz <- lc_penalized(Dr, Er, lambda = l, beta0 = rep(0, A))
  fu <- lc_penalized(Dr, Er, lambda = l, beta0 = B0)
  data.frame(lambda = l,
             dev_beta0_zero = max(abs(fz$beta - f_un$beta)),
             nfree_zero = fz$nfree,
             dev_beta0_unif = max(abs(fu$beta - f_un$beta)),
             nfree_unif = fu$nfree)
}))
print(deg, digits = 3, row.names = FALSE)
write.csv(deg, file.path(TABDIR, "table0_shrinkage_target.csv"), row.names = FALSE)

## ========== 1. 迭代收斂 ====================================================
cat("\n[1] 迭代收斂診斷\n")
lam_tr <- c(0, 1e3, 1e4)
tr <- lapply(lam_tr, function(l)
  lc_penalized(Dr, Er, lambda = l, alpha_en = 1, beta0 = B0, track = TRUE))
names(tr) <- sprintf("lambda=%g", lam_tr)
conv_tab <- do.call(rbind, lapply(seq_along(tr), function(i) data.frame(
  lambda = lam_tr[i], iter = tr[[i]]$iter, converged = tr[[i]]$converged,
  final_obj = tr[[i]]$loglik, nfree = tr[[i]]$nfree,
  last_dbeta = tail(tr[[i]]$history$dbeta, 1))))
print(conv_tab, digits = 4, row.names = FALSE)
write.csv(conv_tab, file.path(TABDIR, "table1_convergence.csv"), row.names = FALSE)

png(file.path(FIGDIR, "fig1_convergence.png"), 1200, 500, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 3, 1))
cols <- c("black", "blue", "red")
gaps <- lapply(tr, function(x) pmax(max(x$history$obj) - x$history$obj, 1e-12))
plot(NA, xlim = c(1, max(sapply(tr, function(x) x$iter))),
     ylim = range(log10(unlist(gaps))),
     xlab = "Iteration", ylab = expression(log[10](f^"*" - f)),
     main = "Objective gap")
for (i in seq_along(tr)) lines(tr[[i]]$history$iter, log10(gaps[[i]]),
                               type = "b", pch = 16, cex = .5, col = cols[i])
legend("topright", names(tr), col = cols, lwd = 2, bty = "n", cex = .8)

plot(NA, xlim = c(1, max(sapply(tr, function(x) x$iter))),
     ylim = range(log10(pmax(unlist(lapply(tr, function(x) x$history$dbeta)), 1e-17))),
     xlab = "Iteration", ylab = expression(log[10]~max~abs(Delta~beta)),
     main = "Parameter change")
for (i in seq_along(tr)) lines(tr[[i]]$history$iter,
                               log10(pmax(tr[[i]]$history$dbeta, 1e-17)),
                               type = "b", pch = 16, cex = .5, col = cols[i])
abline(h = log10(1e-8), lty = 3)
par(op); dev.off()

## ========== 2. lambda 路徑與稀疏性 =========================================
cat("\n[2] lambda 路徑\n")
path_L1 <- lc_pen_path(Dr, Er, alpha_en = 1,   beta0 = B0, nlam = 18)
path_EN <- lc_pen_path(Dr, Er, alpha_en = 0.5, beta0 = B0, nlam = 18)
path_L2 <- lc_pen_path(Dr, Er, alpha_en = 0,   beta0 = B0, nlam = 18)
print(head(path_L1[, c("lambda","frac","deviance","nfree","bic","iter")], 6),
      digits = 4, row.names = FALSE)
write.csv(path_L1, file.path(TABDIR, "table2_path_lasso.csv"), row.names = FALSE)

## beta 的係數路徑
lams_p <- path_L1$lambda
bmat <- sapply(lams_p, function(l)
  lc_penalized(Dr, Er, lambda = l, alpha_en = 1, beta0 = B0)$beta)

png(file.path(FIGDIR, "fig2_path.png"), 1300, 480, res = 120)
op <- par(mfrow = c(1, 3), mar = c(4, 4.2, 3, 1))
matplot(log10(lams_p), t(bmat), type = "l", lty = 1, col = rainbow(A, end = .75),
        xlab = expression(log[10]~lambda), ylab = expression(beta[x]),
        main = "Coefficient path (LASSO)")
abline(h = 1 / A, lty = 2)
plot(log10(path_L1$lambda), path_L1$nfree, type = "b", pch = 16,
     xlab = expression(log[10]~lambda), ylab = "free coefficients",
     main = "Sparsity", ylim = c(0, A))
lines(log10(path_EN$lambda), path_EN$nfree, type = "b", pch = 17, col = "blue")
legend("bottomleft", c("LASSO", "Elastic Net (a=.5)"), col = c("black","blue"),
       pch = c(16,17), bty = "n", cex = .8)
plot(log10(path_L1$lambda), path_L1$bic, type = "b", pch = 16,
     xlab = expression(log[10]~lambda), ylab = "BIC", main = "BIC")
abline(v = log10(path_L1$lambda[which.min(path_L1$bic)]), col = 2, lty = 2)
par(op); dev.off()

## ========== 3. 偏誤—變異數權衡 =============================================
cat("\n[3] 偏誤—變異數權衡（reps =", REPS, "）\n")
# lambda 以 lambda_max 的比例指定：Hessian 尺度正比於 N，
# 絕對尺度的 lambda 在不同暴露數之間不可比。
FRACS <- c(0, 0.003, 0.01, 0.03, 0.1, 0.2, 0.4, 0.7, 1)
bv <- list()
for (N in c(5e4, 2e5)) {
  cat(sprintf(" N = %.0e\n", N))
  bv[[as.character(N)]] <- cbind(N = N,
    pen_bias_variance(truth, N, FRACS, alpha_en = 1, reps = REPS, beta0 = B0))
}
bv_all <- do.call(rbind, bv)
print(bv_all, digits = 4, row.names = FALSE)
write.csv(bv_all, file.path(TABDIR, "table3_bias_variance.csv"), row.names = FALSE)

png(file.path(FIGDIR, "fig3_bias_variance.png"), 1200, 500, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 3, 1))
for (N in c(5e4, 2e5)) {
  s <- bv_all[bv_all$N == N, ]
  x <- s$frac
  FLOOR <- 1e-7
  plot(x, s$mse, type = "b", pch = 16, log = "y",
       ylim = range(pmax(c(s$mse, s$variance, s$bias2), FLOOR)),
       xlab = expression(lambda/lambda[max]), ylab = "SSE(beta) component",
       main = sprintf("N = %.0e", N))
  lines(x, pmax(s$variance, FLOOR), type = "b", pch = 17, col = "blue")
  lines(x, pmax(s$bias2, FLOOR), type = "b", pch = 15, col = "red")
  abline(v = x[which.min(s$mse)], lty = 2)
  legend("bottomleft", c("MSE", "Variance", "Bias^2"),
         col = c("black","blue","red"), pch = c(16,17,15), bty = "n", cex = .75)
}
par(op); dev.off()

## ========== 4. 變異數的收斂速度 ============================================
cat("\n[4] 變異數收斂速度（reps =", REPS, "）\n")
Ns <- c(2e4, 5e4, 1e5, 2e5, 5e5, 1e6)
rate0 <- pen_variance_rate(truth, Ns, frac = 0,    reps = REPS, beta0 = B0)
FR <- 0.1
rate1 <- pen_variance_rate(truth, Ns, frac = FR, reps = REPS, beta0 = B0)
rate_tab <- rbind(cbind(method = "未懲罰", rate0),
                  cbind(method = sprintf("LASSO(frac=%.2f)", FR), rate1))
cat(sprintf("\n未懲罰   ：log-log 斜率 = %.3f (SE %.3f)\n",
            attr(rate0,"slope"), attr(rate0,"slope_se")))
cat(sprintf("LASSO    ：log-log 斜率 = %.3f (SE %.3f)\n",
            attr(rate1,"slope"), attr(rate1,"slope_se")))
print(rate_tab, digits = 4, row.names = FALSE)
write.csv(rate_tab, file.path(TABDIR, "table4_variance_rate.csv"), row.names = FALSE)

png(file.path(FIGDIR, "fig4_variance_rate.png"), 1200, 500, res = 120)
op <- par(mfrow = c(1, 2), mar = c(4, 4.2, 3, 1))
plot(rate0$N, rate0$variance, type = "b", pch = 16, log = "xy",
     ylim = range(c(rate0$variance, rate1$variance)),
     xlab = "Exposure N", ylab = "sum Var(beta)", main = "Variance convergence")
lines(rate1$N, rate1$variance, type = "b", pch = 17, col = "red")
lines(rate0$N, rate0$variance[1] * (rate0$N / rate0$N[1])^(-1), lty = 3)
legend("bottomleft", c("unpenalized","LASSO","slope -1 ref"),
       col = c("black","red","black"), lty = c(1,1,3), pch = c(16,17,NA),
       bty = "n", cex = .75)
plot(rate0$N, rate0$mse, type = "b", pch = 16, log = "xy",
     ylim = range(c(rate0$mse, rate1$mse)),
     xlab = "Exposure N", ylab = "MSE(beta)", main = "MSE")
lines(rate1$N, rate1$mse, type = "b", pch = 17, col = "red")
legend("bottomleft", c("unpenalized","LASSO"), col = c("black","red"),
       pch = c(16,17), bty = "n", cex = .75)
par(op); dev.off()

## ========== 5. LASSO vs Ridge vs Elastic Net ==============================
cat("\n[5] 三種懲罰的比較（N = 2e5, reps =", REPS, "）\n")
cmp <- do.call(rbind, lapply(c(1, 0.5, 0), function(a) {
  r <- pen_bias_variance(truth, 2e5, FRACS, alpha_en = a, reps = REPS,
                         beta0 = B0, verbose = FALSE)
  best <- r[which.min(r$mse), ]
  data.frame(penalty = c("1"="LASSO","0.5"="Elastic Net","0"="Ridge")[as.character(a)],
             alpha_en = a, best_frac = best$frac, bias2 = best$bias2,
             variance = best$variance, mse = best$mse, cor_med = best$cor_med)
}))
unp <- bv_all[bv_all$N == 2e5 & bv_all$frac == 0, ]
cmp <- rbind(data.frame(penalty = "未懲罰", alpha_en = NA, best_frac = 0,
                        bias2 = unp$bias2, variance = unp$variance,
                        mse = unp$mse, cor_med = unp$cor_med), cmp)
print(cmp, digits = 4, row.names = FALSE)
write.csv(cmp, file.path(TABDIR, "table5_penalty_comparison.csv"), row.names = FALSE)

cat(sprintf("\n圖已輸出至 %s/、表已輸出至 %s/\n", FIGDIR, TABDIR))
