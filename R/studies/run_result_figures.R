###############################################################################
# 實證分析用圖：各估計量的比較，以及對平均餘命的後果
#
#   figT3：三項指標（alpha 最大偏誤、SSE(beta)、|漂移偏誤|）隨人口規模的變化，
#          實線為卜瓦松模擬、虛線為實際觀測死亡數的二項抽薄。兩組並列即可看出
#          結論是否因模型誤設而改變。縱軸取對數，因各估計量相差達兩個數量級。
#   figT4：零歲平均餘命的重複分布（各估計量，N = 5e4），
#          把參數層次的改善翻譯成使用者真正關心的量。
#
# 輸出：output/figures/figT3_estimator_compare.png
#       output/figures/figT4_e0_by_estimator.png
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/fig_axis_utils.R")

## ---------- figT3：由既有表格作圖 ----------
sim  <- read.csv("output/tables/tableT_analytic_bias.csv", check.names = FALSE)
thin <- read.csv("output/tables/tableW_real_thinning.csv", check.names = FALSE)
map <- c("A  標準 LC"="標準 LC", "A 標準 LC"="標準 LC",
         "B  卜瓦松 MLE"="卜瓦松 MLE", "B 卜瓦松 MLE"="卜瓦松 MLE",
         "C  Firth"="Firth", "C Firth"="Firth",
         "D  解析校正 shrink=1.0"="中心化", "D 中心化"="中心化",
         "G  解析校正 + 加權 SVD"="中心化＋加權", "E 中心化 + 加權"="中心化＋加權")
sim$g <- map[sim$估計量]; thin$g <- map[thin$估計量]
sim <- sim[!is.na(sim$g), ]; thin <- thin[!is.na(thin$g), ]
KEEP <- c("標準 LC","卜瓦松 MLE","Firth","中心化","中心化＋加權")
COL  <- c("#B03A2E","#7D3C98","#D68910","#1F6FB2","#1E8449"); names(COL) <- KEEP
NS <- c(1e4, 5e4, 2e5)
METR <- list(list(col="alpha最大", lab=expression(paste(alpha," 最大偏誤"))),
             list(col="SSE_beta", lab=expression(paste("SSE(", beta, ")"))),
             list(col="漂移偏誤", lab="|漂移項偏誤|"))

png_cjk("output/figures/figT3_estimator_compare.png", width=2250, height=820, res = 198)
par(mfrow=c(1,3), mar=c(5.0,4.8,3.4,1.0), mgp=c(3.0,0.8,0))
for (mm in METR) {
  vals <- c(abs(sim[[mm$col]]), abs(thin[[mm$col]]))
  vals <- vals[is.finite(vals) & vals > 0]
  plot(NA, xlim=range(log10(NS)), ylim=range(vals), log="y", xaxt="n",
       xlab="人口規模", ylab=mm$lab, main=mm$lab)
  axis(1, at=log10(NS), labels=format(NS, big.mark=",", scientific=FALSE))
  abline(v=log10(NS), col="grey93", lty=3)
  for (g in KEEP) {
    a <- sim[sim$g==g, ]; a <- a[order(a$N), ]
    b <- thin[thin$g==g, ]; b <- b[order(b$N), ]
    lines(log10(a$N), abs(a[[mm$col]]), col=COL[g], lwd=2.4, type="b", pch=16, cex=0.9)
    if (nrow(b)) lines(log10(b$N), abs(b[[mm$col]]), col=COL[g], lwd=1.6, lty=2,
                       type="b", pch=1, cex=0.9)
  }
  if (identical(mm$col, "alpha最大"))
    legend("bottomleft", bty="n", cex=0.74, legend=KEEP, col=COL[KEEP], lwd=2.4, pch=16)
  if (identical(mm$col, "SSE_beta"))
    legend("bottomleft", bty="n", cex=0.74,
           legend=c("卜瓦松模擬（實線、實心）","真實資料抽薄（虛線、空心）"),
           col="grey30", lwd=c(2.4,1.6), lty=c(1,2), pch=c(16,1))
}
dev.off()

## ---------- figT4：平均餘命 ----------
SEED <- 20260926; REPS <- 100; YEARS <- 2001:2024
dat <- load_data(sex="Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]; A<-nrow(D0); Tn<-ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage <- E0[,Tn]/sum(E0[,Tn]); e0_true <- life_table(mtrue[,Tn])$e0
BF <- make_logbias(0.5)
fit_w <- function(D,E) {
  g <- lc_analytic(D,E,bfun=BF,shrink=1)
  mu <- E*exp(outer(g$a,rep(1,ncol(D)))+outer(g$b,g$k))
  lmc <- log(pmax(D,0.5)/E) - BF(mu); W <- mu/mean(mu)
  am <- rowSums(W*lmc)/rowSums(W); Z <- (lmc-am)*sqrt(W)
  sv <- svd(Z); u1<-sv$u[,1]; v1<-sv$v[,1]
  if (sum(u1)<0){u1<--u1; v1<--v1}
  bb<-u1/sum(u1); kk<-sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
  kk<-kk-mean(kk); sb<-sum(bb); list(a=am,b=bb/sb,k=kk*sb)
}
FIT <- list("標準 LC"=function(D,E) lc_svd_fit(log(pmax(D,0.5)/E)),
            "卜瓦松 MLE"=function(D,E) lc_poisson_firth(D,E,firth=FALSE),
            "Firth"=function(D,E) lc_poisson_firth(D,E,firth=TRUE),
            "中心化"=function(D,E) lc_analytic(D,E,bfun=BF,shrink=1),
            "中心化＋加權"=fit_w)
N <- 5e4; E <- outer(N*wage, rep(1,Tn)); dimnames(E) <- dimnames(D0)
EE <- matrix(NA_real_, REPS, length(FIT)); colnames(EE) <- names(FIT)
for (r in seq_len(REPS)) {
  set.seed(SEED+1000*2+r)
  D <- matrix(rpois(length(E), E*mtrue), A, dimnames=dimnames(E))
  for (j in seq_along(FIT)) {
    f <- tryCatch(FIT[[j]](D,E), error=function(e) NULL)
    if (is.null(f)||!all(is.finite(f$a))||!all(is.finite(f$b))) next
    mh <- exp(outer(f$a,rep(1,Tn))+outer(f$b,f$k))
    EE[r,j] <- tryCatch(life_table(mh[,Tn])$e0, error=function(e) NA_real_)
  }
  if (r%%25==0) cat(sprintf("  e0 rep %d/%d\n", r, REPS))
}
png_cjk("output/figures/figT4_e0_by_estimator.png", width=1700, height=820, res = 198)
par(mar=c(6.4,4.8,3.4,1.0), mgp=c(3.0,0.8,0))
K <- length(FIT)
plot(NA, xlim=c(0.5,K+0.5), ylim=range(c(EE,e0_true), na.rm=TRUE), xaxt="n",
     xlab="", ylab="零歲平均餘命（歲）",
     main=sprintf("配適末年的平均餘命（N = 50,000，重複 100 次；真值 %.2f 歲）", e0_true))
axis(1, at=seq_len(K), labels=names(FIT), las=2, cex.axis=0.85)
abline(h=e0_true, col="black", lwd=2.4, lty=2)
for (j in seq_len(K)) {
  v <- EE[,j]; v <- v[is.finite(v)]; if (!length(v)) next
  points(jitter(rep(j,length(v)), amount=0.13), v, pch=16, cex=0.6,
         col=grDevices::rgb(t(col2rgb(COL[min(j,length(COL))]))/255, alpha=0.45))
  q <- quantile(v, c(0.05,0.5,0.95))
  segments(j, q[1], j, q[3], col=COL[min(j,length(COL))], lwd=2.0)
  segments(j-0.24, q[2], j+0.24, q[2], col=COL[min(j,length(COL))], lwd=3.4)
  text(j, max(EE, na.rm=TRUE), sprintf("%+.2f", q[2]-e0_true), cex=0.82,
       col=COL[min(j,length(COL))], font=2)
}
legend("bottomright", bty="n", cex=0.8, legend=c("真值","中位數與 5–95%"),
       col=c("black","grey30"), lwd=c(2.4,2.4), lty=c(2,1))
dev.off()
cat("\n=== e0 中位偏誤（N=5e4）===\n")
for (j in seq_len(K)) cat(sprintf("  %-14s %+.3f 歲（5-95%% 全距 %.2f）\n",
  names(FIT)[j], median(EE[,j],na.rm=TRUE)-e0_true,
  diff(quantile(EE[,j], c(0.05,0.95), na.rm=TRUE))))
cat("\n已輸出 figT3 / figT4\n")
