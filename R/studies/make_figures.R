###############################################################################
# 產生報告使用的出版用圖
#
# 各 study 腳本負責跑模擬並輸出表（output/tables/），本檔則由這些表重繪圖，
# 統一套用 R/core/fig_axis_utils.R 的座標軸格式——把 R 預設的電腦式
# 科學記號（1e-03、5e+04）改為學術書寫的 10^{-3}、5 x 10^{4}。
#
# 用法（工作目錄為 repo 根目錄）：  Rscript R/studies/make_figures.R
###############################################################################

source("R/core/fig_axis_utils.R")
FIGDIR <- "output/figures"; TABDIR <- "output/tables"

## ---------- mc_plots.png ----------
source("R/studies/mc_exposure.R")
res <- read.csv("output/tables/mc_results.csv"); attr(res,"drift_true") <- -0.367781
ms <- unique(res$method); cols <- setNames(seq_along(ms), ms)
pick <- function(m, col) res[[col]][res$method == m]
Ns <- sort(unique(res$N))
png(file.path(FIGDIR, "mc_plots.png"), 1200, 950, res = 120)
op <- par(mfrow = c(2,2), mar = c(4,4.4,3,1))
pnl <- function(col, main, ylab, logy = TRUE, href = NULL) {
  yl <- range(res[[col]], na.rm = TRUE, finite = TRUE)
  if (!is.null(href)) yl <- range(c(yl, href))
  plot(Ns, pick(ms[1], col), type = "n", log = if (logy) "xy" else "x",
       ylim = yl, xlab = "Exposure N", ylab = ylab, main = main, xaxt = "n",
       yaxt = if (logy) "n" else "s")
  axis_pow(1, Ns, cex.axis = .72)
  if (logy) axis_pow(2, pow_ticks(res[[col]]))
  for (m in ms) lines(Ns, pick(m, col), type = "b", pch = 16, col = cols[m])
}
pnl("b_mse_med", "median SSE(beta)", "median SSE(beta)")
legend("bottomleft", ms, col = cols, lwd = 1, pch = 16, bty = "n", cex = .7)
pnl("a_bias", "mean bias(alpha)", "mean bias(alpha)", logy = FALSE); abline(h = 0, lty = 2)
dt <- attr(res,"drift_true")
pnl("drift_med", sprintf("drift (true = %.4f)", dt), "median drift", logy = FALSE, href = dt)
abline(h = dt, lty = 2, col = 2)
plot(Ns, pick(ms[1],"breakdown"), type="n", log="x", ylim=c(0,1), xaxt="n",
     xlab="Exposure N", ylab="P(cor(beta) < 0.8)", main="breakdown rate")
axis_pow(1, Ns, cex.axis = .72)
for (m in ms) lines(Ns, pick(m,"breakdown"), type="b", pch=16, col=cols[m])
par(op); dev.off()

## ---------- fig3_bias_variance.png ----------
bv <- read.csv(file.path(TABDIR,"table3_bias_variance.csv")); FLOOR <- 1e-7
png(file.path(FIGDIR,"fig3_bias_variance.png"), 1200, 500, res = 120)
op <- par(mfrow=c(1,2), mar=c(4,4.6,3,1))
for (N in sort(unique(bv$N))) {
  s <- bv[bv$N==N,]; x <- s$frac
  yv <- pmax(c(s$mse,s$variance,s$bias2), FLOOR); yl <- range(yv)
  plot(x, pmax(s$mse,FLOOR), type="b", pch=16, log="y", ylim=yl, yaxt="n",
       xlab=expression(lambda/lambda[max]), ylab="SSE(beta) component",
       main=bquote(N == .(signif(N/10^floor(log10(N)),2)) %*% 10^.(floor(log10(N)))))
  axis_pow(2, pow_ticks(yv, 6))
  lines(x, pmax(s$variance,FLOOR), type="b", pch=17, col="blue")
  lines(x, pmax(s$bias2,FLOOR),    type="b", pch=15, col="red")
  jb <- which.min(s$mse); abline(v=x[jb], lty=2)
  text(x[jb], yl[2], bquote(lambda/lambda[max] == .(x[jb])), pos=4, cex=.7)
  legend("bottomleft", c("MSE","Variance","Bias^2"), col=c("black","blue","red"),
         pch=c(16,17,15), bty="n", cex=.75)
}
par(op); dev.off()

## ---------- fig4_variance_rate.png ----------
r <- read.csv(file.path(TABDIR,"table4_variance_rate.csv"))
m <- unique(r$method); a <- r[r$method==m[1],]; b <- r[r$method==m[2],]
png(file.path(FIGDIR,"fig4_variance_rate.png"), 1200, 500, res=120)
op <- par(mfrow=c(1,2), mar=c(4,4.6,3,1))
yv <- c(a$variance,b$variance)
plot(a$N, a$variance, type="b", pch=16, log="xy", ylim=range(yv), xaxt="n", yaxt="n",
     xlab="Exposure N", ylab="sum Var(beta)", main="Variance convergence")
axis_pow(1, a$N, cex.axis = .72); axis_pow(2, pow_ticks(yv, 5))
lines(b$N, b$variance, type="b", pch=17, col="red")
lines(a$N, a$variance[1]*(a$N/a$N[1])^(-1), lty=3)
legend("bottomleft", c("unpenalized","LASSO","slope -1 ref"),
       col=c("black","red","black"), lty=c(1,1,3), pch=c(16,17,NA), bty="n", cex=.72)
yv2 <- c(a$mse,b$mse)
plot(a$N, a$mse, type="b", pch=16, log="xy", ylim=range(yv2), xaxt="n", yaxt="n",
     xlab="Exposure N", ylab="MSE(beta)", main="MSE (crossover)")
axis_pow(1, a$N, cex.axis = .72); axis_pow(2, pow_ticks(yv2, 5))
lines(b$N, b$mse, type="b", pch=17, col="red")
ix <- which(diff(sign(a$mse-b$mse))!=0)
if (length(ix)) abline(v=sqrt(a$N[ix]*a$N[ix+1]), lty=2, col="grey40")
legend("bottomleft", c("unpenalized","LASSO"), col=c("black","red"), pch=c(16,17), bty="n", cex=.75)
par(op); dev.off()

## ---------- figA_beta_by_age.png（右下 SD 面板為 log 軸）----------
source("R/core/penalized_lc.R")
t1 <- read.csv(file.path(TABDIR,"tableA1_beta_by_age_real.csv"))
a2 <- read.csv(file.path(TABDIR,"tableA2_beta_by_age_sim.csv"))
A <- nrow(t1); mid <- make_grid(A)$mid
png(file.path(FIGDIR,"figA_beta_by_age.png"), 1400, 900, res=120)
op <- par(mfrow=c(2,2), mar=c(4,4.6,3,1))
matplot(mid, cbind(t1$beta_unpen,t1$beta_lasso,t1$beta_enet,t1$beta_ridge), type="b",
        pch=c(16,17,15,18), lty=1, col=c("black","red","blue","darkgreen"),
        xlab="Age", ylab=expression(beta[x]), main="Real data (Female 2001-2024)")
abline(h=1/A, lty=2, col="grey40")
legend("topright", c("unpenalized","LASSO","Elastic Net","Ridge","1/A"),
       col=c("black","red","blue","darkgreen","grey40"), pch=c(16,17,15,18,NA),
       lty=c(1,1,1,1,2), bty="n", cex=.7)
for (N in sort(unique(a2$N))) {
  s <- a2[a2$N==N,]
  matplot(mid, cbind(s$beta_true,s$mean_unpen,s$mean_pen), type="b",
          pch=c(NA,16,17), lty=1, lwd=c(2,1,1), col=c("grey30","black","red"),
          xlab="Age", ylab=expression(beta[x]),
          main=bquote("Simulation," ~ N == .(signif(N/10^floor(log10(N)),2)) %*% 10^.(floor(log10(N)))))
  abline(h=1/A, lty=2, col="grey60")
  legend("topright", c("truth","unpenalized (mean)","LASSO (mean)"),
         col=c("grey30","black","red"), lty=1, pch=c(NA,16,17), bty="n", cex=.7)
}
s <- a2[a2$N==min(a2$N),]; yv <- pmax(c(s$sd_unpen,s$sd_pen), 1e-6)
plot(mid, s$sd_unpen, type="b", pch=16, log="y", ylim=range(yv), yaxt="n",
     xlab="Age", ylab="SD of beta_x",
     main=bquote("Per-age SD," ~ N == 5 %*% 10^4))
axis_pow(2, pow_ticks(yv, 5))
lines(mid, pmax(s$sd_pen,1e-6), type="b", pch=17, col="red")
legend("bottomleft", c("unpenalized","LASSO"), col=c("black","red"), pch=c(16,17), bty="n", cex=.75)
par(op); dev.off()
cat(sprintf("\n出版用圖已重繪至 %s/\n", FIGDIR))
