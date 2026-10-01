###############################################################################
# 緒論用圖：小人口下 Lee-Carter 估計的震盪與偏誤
#
#   目的。既有研究（王信忠等 2012；Wang et al. 2018；余清祥等 2021）指出
#   小人口的死亡率估計同時有「震盪」與「有方向的偏誤」兩個問題。
#   本圖以台灣女性資料重現之，三格分別呈現震盪的三個層次：
#     (a) 原始資料：單次實現的年齡別死亡率曲線，與全國的平滑曲線對照；
#         零死亡格另以空心點標出其被替代後的位置。
#     (b) 參數：beta_x 的多次重複軌跡，顯示形狀完全不穩定。
#     (c) 標的量：零歲平均餘命的重複分布，顯示震盪\textbf{與}系統性低估並存。
#
# 輸出：output/figures/figT0_oscillation.png
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/fig_axis_utils.R")

SEED <- 20260926; YEARS <- 2001:2024
dat <- load_data(sex="Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]; A<-nrow(D0); Tn<-ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage  <- E0[,Tn]/sum(E0[,Tn]); ages <- rownames(D0); xx <- seq_len(A)
e0_true <- life_table(mtrue[,Tn])$e0
cat(sprintf("真值 e0（配適末年）= %.3f\n", e0_true))

sim <- function(N, seed) {
  E <- outer(N*wage, rep(1,Tn)); dimnames(E) <- dimnames(D0)
  set.seed(seed); matrix(rpois(length(E), E*mtrue), A, dimnames=dimnames(E))
}
expo <- function(N) { E <- outer(N*wage, rep(1,Tn)); dimnames(E) <- dimnames(D0); E }

png_cjk("output/figures/figT0_oscillation.png", width=2250, height=800, res = 198)
par(mfrow=c(1,3), mar=c(5.2,4.6,3.4,1.0), mgp=c(2.7,0.8,0))

## ---- (a) 原始資料的震盪 ----
tsel <- Tn
plot(NA, xlim=c(1,A), ylim=c(-11,-0.5), xaxt="n", xlab="年齡組",
     ylab=expression(log~m[x]), main="(a) 年齡別死亡率：單次實現")
axis(1, at=xx, labels=ages, las=2, cex.axis=0.62)
abline(v=xx, col="grey93", lty=3)
cols <- c("#B03A2E", "#1F6FB2")
for (i in seq_along(c(1e4, 5e4))) {
  N <- c(1e4,5e4)[i]; D <- sim(N, SEED+9000+i); E <- expo(N)
  y <- log(pmax(D[,tsel],0.5)/E[,tsel])
  lines(xx, y, col=cols[i], lwd=1.7, type="b", pch=16, cex=0.55)
  z <- which(D[,tsel]==0)
  if (length(z)) points(z, y[z], pch=21, bg="white", col=cols[i], cex=1.25, lwd=1.7)
}
lines(xx, log(mtrue[,tsel]), col="black", lwd=2.6)
legend("bottomright", bty="n", cex=0.8,
  legend=c("全國（真值）", "N = 10,000", "N = 50,000", "零死亡格（以 0.5 替代）"),
  col=c("black", cols, "grey30"), lwd=c(2.6,1.7,1.7,NA),
  pch=c(NA,16,16,21), pt.bg=c(NA,NA,NA,"white"))

## ---- (b) beta_x 的重複軌跡 ----
R20 <- 20; N <- 5e4; E <- expo(N)
Bm <- matrix(NA_real_, R20, A)
for (r in seq_len(R20)) {
  set.seed(SEED+1000*2+r)
  D <- matrix(rpois(length(E), E*mtrue), A, dimnames=dimnames(E))
  Bm[r,] <- lc_svd_fit(log(pmax(D,0.5)/E))$b
}
plot(NA, xlim=c(1,A), ylim=c(-0.18,0.38), xaxt="n", xlab="年齡組",
     ylab=expression(hat(beta)[x]), main="(b) 標準 LC 的 beta：20 次重複（N = 50,000）")
axis(1, at=xx, labels=ages, las=2, cex.axis=0.62)
abline(h=0, col="grey40"); abline(v=xx, col="grey93", lty=3)
for (r in seq_len(R20)) lines(xx, Bm[r,], col=grDevices::rgb(0.65,0.25,0.18,0.42), lwd=1.0)
lines(xx, truth$b, col="black", lwd=2.8); points(xx, truth$b, pch=16, cex=0.7)
legend("topleft", bty="n", cex=0.8, legend=c("真值","各次重複"),
       col=c("black", grDevices::rgb(0.65,0.25,0.18,0.7)), lwd=c(2.8,1.0))

## ---- (c) e0 的重複分布 ----
NS <- c(1e4, 5e4, 2e5); R100 <- 100
E0s <- matrix(NA_real_, R100, length(NS))
for (i in seq_along(NS)) {
  N <- NS[i]; E <- expo(N)
  for (r in seq_len(R100)) {
    set.seed(SEED+1000*i+r)
    D <- matrix(rpois(length(E), E*mtrue), A, dimnames=dimnames(E))
    f <- lc_svd_fit(log(pmax(D,0.5)/E))
    mh <- exp(outer(f$a,rep(1,Tn))+outer(f$b,f$k))
    E0s[r,i] <- tryCatch(life_table(mh[,Tn])$e0, error=function(e) NA_real_)
  }
}
plot(NA, xlim=c(0.5,length(NS)+0.5), ylim=range(c(E0s,e0_true), na.rm=TRUE),
     xaxt="n", xlab="人口規模", ylab="零歲平均餘命（歲）",
     main="(c) 標準 LC 的平均餘命：100 次重複")
axis(1, at=seq_along(NS), labels=format(NS, big.mark=",", scientific=FALSE))
abline(h=e0_true, col="black", lwd=2.4, lty=2)
for (i in seq_along(NS)) {
  v <- E0s[,i]; v <- v[is.finite(v)]
  points(jitter(rep(i,length(v)), amount=0.14), v,
         pch=16, cex=0.55, col=grDevices::rgb(0.12,0.44,0.70,0.45))
  segments(i-0.26, median(v), i+0.26, median(v), col="#B03A2E", lwd=3)
  cat(sprintf("N=%.0e  e0 中位 %.2f（真值 %.2f，偏誤 %+.2f）  5-95%%範圍 %.2f 歲\n",
      NS[i], median(v), e0_true, median(v)-e0_true,
      diff(quantile(v, c(0.05,0.95)))))
}
legend("bottomright", bty="n", cex=0.8,
       legend=c("真值","各次重複","中位數"),
       col=c("black", grDevices::rgb(0.12,0.44,0.70,0.7), "#B03A2E"),
       lwd=c(2.4,NA,3), lty=c(2,NA,1), pch=c(NA,16,NA))
dev.off()
cat("\n已輸出 figT0_oscillation.png\n")
