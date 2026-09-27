###############################################################################
# 方法的極限：人口規模小到什麼程度時中心化失效
#
#   可事先寫下的失效機制。中心化須代入 mu_hat，而 mu_hat 由配適值算出。
#   若某年齡 x 的各年幾乎全為零死亡，則 ell_{x,t} = log(c0/E) 對所有 t 成立，
#   故 alpha_hat_x 被釘在 log(c0/E) 附近，於是
#       mu_hat_{x,t} = E exp(alpha_hat + beta_hat kappa_hat) ≈ c0,
#   而 b(c0) 是一個有限的小數（c0=0.5 時 b(0.5)=0.342）。
#   但真實的偏誤為 log(c0/mu)，在 mu=0.02 時高達 3.2。
#   亦即\textbf{校正在最需要它的時候會自己關掉}，因為它所依據的 mu_hat
#   已無資訊可用。據此預測：存在一個人口規模，低於該規模後中心化的
#   改善幅度迅速消失。
#
#   本檔掃描 N = 2000 至 1e6，回報三項：alpha 最大偏誤、零歲平均餘命的
#   偏誤、以及「整列全零」的年齡數（失效機制的直接指標）。
#
# 輸出：output/tables/tableX_size_limit.csv
#       output/figures/figT5_size_limit.png
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/fig_axis_utils.R")

SEED <- 20260926; REPS <- 100; YEARS <- 2001:2024
NS <- c(2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 1e6)
dat <- load_data(sex="Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]; A<-nrow(D0); Tn<-ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage <- E0[,Tn]/sum(E0[,Tn]); e0_true <- life_table(mtrue[,Tn])$e0
BF <- make_logbias(0.5)
cat(sprintf("真值 e0 = %.3f\n", e0_true))

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
            "Firth"=function(D,E) lc_poisson_firth(D,E,firth=TRUE),
            "中心化＋加權"=fit_w)
K <- length(FIT); res <- list()

for (N in NS) {
  E <- outer(N*wage, rep(1,Tn)); dimnames(E) <- dimnames(D0); mu <- E*mtrue
  aM <- array(NA_real_,c(REPS,A,K)); e0M <- matrix(NA_real_,REPS,K)
  bad <- numeric(REPS); corr_used <- numeric(REPS)
  for (r in seq_len(REPS)) {
    set.seed(SEED + 300*which(NS==N) + r)
    D <- matrix(rpois(length(E), mu), A, dimnames=dimnames(E))
    bad[r] <- sum(rowSums(D)==0)                      # 整列全零的年齡數
    ## 校正實際用到的量 vs 真正該扣的量（失效機制的直接度量）
    g <- tryCatch(lc_analytic(D,E,bfun=BF,shrink=1), error=function(e) NULL)
    if (!is.null(g)) {
      mh <- E*exp(outer(g$a,rep(1,Tn))+outer(g$b,g$k))
      corr_used[r] <- max(abs(rowMeans(BF(mh))))      # 實際扣的最大值
    }
    for (j in seq_len(K)) {
      f <- tryCatch(FIT[[j]](D,E), error=function(e) NULL)
      if (is.null(f)||!all(is.finite(f$a))||!all(is.finite(f$b))) next
      aM[r,,j] <- f$a-truth$a
      mh <- exp(outer(f$a,rep(1,Tn))+outer(f$b,f$k))
      e0M[r,j] <- tryCatch(life_table(mh[,Tn])$e0, error=function(e) NA_real_)
    }
  }
  need <- max(abs(rowMeans(BF(mu))))                  # 真正該扣的最大值
  res[[as.character(N)]] <- data.frame(
    N=N, 全零年齡數=round(mean(bad),1),
    該扣最大=round(need,3), 實扣最大=round(median(corr_used),3),
    校正達成率=round(median(corr_used)/need,3),
    setNames(as.data.frame(lapply(seq_len(K), function(j)
      round(max(abs(apply(aM[,,j],2,median,na.rm=TRUE))),4))), paste0("a最大_",names(FIT))),
    setNames(as.data.frame(lapply(seq_len(K), function(j)
      round(median(e0M[,j],na.rm=TRUE)-e0_true,3))), paste0("e0偏誤_",names(FIT))),
    check.names=FALSE)
  cat(sprintf("N=%7.0f 全零年齡 %4.1f  該扣 %.2f 實扣 %.2f (%.0f%%) | a最大 %6.3f/%6.3f/%6.3f | e0偏誤 %+6.2f/%+6.2f/%+6.2f\n",
    N, mean(bad), need, median(corr_used), 100*median(corr_used)/need,
    max(abs(apply(aM[,,1],2,median,na.rm=TRUE))), max(abs(apply(aM[,,2],2,median,na.rm=TRUE))),
    max(abs(apply(aM[,,3],2,median,na.rm=TRUE))),
    median(e0M[,1],na.rm=TRUE)-e0_true, median(e0M[,2],na.rm=TRUE)-e0_true,
    median(e0M[,3],na.rm=TRUE)-e0_true))
}
tab <- do.call(rbind,res)
dir.create("output/tables",recursive=TRUE,showWarnings=FALSE)
write.csv(tab,"output/tables/tableX_size_limit.csv",row.names=FALSE)

## ---- 圖 ----
COL <- c("#B03A2E","#D68910","#1E8449")
png_cjk("output/figures/figT5_size_limit.png", width=2250, height=820, res=150)
par(mfrow=c(1,3), mar=c(5.0,4.8,3.4,1.0), mgp=c(3.0,0.8,0))
plot(NA, xlim=range(log10(NS)), ylim=range(tab[,grep("^a最大",names(tab))]), log="y",
     xaxt="n", xlab="人口規模", ylab=expression(paste(alpha," 最大偏誤")),
     main=expression(paste(alpha," 最大偏誤")))
axis(1,at=log10(NS),labels=format(NS,big.mark=",",scientific=FALSE),cex.axis=0.8,las=2)
abline(v=log10(NS),col="grey93",lty=3)
for (j in seq_len(K)) lines(log10(tab$N), tab[[paste0("a最大_",names(FIT)[j])]],
  col=COL[j], lwd=2.4, type="b", pch=16)
legend("bottomleft", bty="n", cex=0.8, legend=names(FIT), col=COL, lwd=2.4, pch=16)

plot(NA, xlim=range(log10(NS)), ylim=range(tab[,grep("^e0偏誤",names(tab))]),
     xaxt="n", xlab="人口規模", ylab="平均餘命偏誤（歲）", main="零歲平均餘命偏誤")
axis(1,at=log10(NS),labels=format(NS,big.mark=",",scientific=FALSE),cex.axis=0.8,las=2)
abline(h=0,col="black",lwd=1.6,lty=2); abline(v=log10(NS),col="grey93",lty=3)
for (j in seq_len(K)) lines(log10(tab$N), tab[[paste0("e0偏誤_",names(FIT)[j])]],
  col=COL[j], lwd=2.4, type="b", pch=16)

plot(log10(tab$N), 100*tab$校正達成率, type="b", pch=16, lwd=2.4, col="#1F6FB2",
     ylim=c(0,110), xaxt="n", xlab="人口規模", ylab="校正達成率（%）",
     main="校正實際扣除的量 / 應扣的量")
axis(1,at=log10(NS),labels=format(NS,big.mark=",",scientific=FALSE),cex.axis=0.8,las=2)
abline(h=100,col="grey40",lty=2); abline(v=log10(NS),col="grey93",lty=3)
par(new=TRUE)
plot(log10(tab$N), tab$全零年齡數, type="b", pch=1, lty=2, col="#B03A2E",
     axes=FALSE, xlab="", ylab="", ylim=c(0,max(tab$全零年齡數)*1.15))
axis(4, col="#B03A2E", col.axis="#B03A2E"); mtext("整列全零的年齡數", 4, 2.4, col="#B03A2E")
legend("right", bty="n", cex=0.78, legend=c("校正達成率（左軸）","全零年齡數（右軸）"),
       col=c("#1F6FB2","#B03A2E"), lwd=c(2.4,1.6), lty=c(1,2), pch=c(16,1))
dev.off()
cat("\n已輸出 tableX / figT5\n")
