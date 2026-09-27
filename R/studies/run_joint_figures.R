###############################################################################
# 實證分析用圖：alpha 與 beta 能否同時改善
#   figT6：alpha--beta 的取捨圖（Pareto 圖）。橫軸為 SSE(beta)、縱軸為 alpha
#          最大偏誤，皆取對數；愈靠左下愈好。同時改善即「往左下移動」。
#   figT7：僅校正 alpha ＋ 加權 與標準 LC、Firth 的逐年齡帶（alpha 與 beta）。
# 輸出：figT6_pareto.png / figT7_joint_band.png
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/fig_axis_utils.R")

SEED <- 20260926; REPS <- 100; NS <- c(1e4,5e4,2e5); YEARS <- 2001:2024
dat <- load_data(sex="Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]; A<-nrow(D0); Tn<-ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage <- E0[,Tn]/sum(E0[,Tn]); ages <- rownames(D0); xx <- seq_len(A)
BF <- make_logbias(0.5)
fit_cellw <- function(D,E){ g<-lc_analytic(D,E,bfun=BF,shrink=1)
  mu<-E*exp(outer(g$a,rep(1,ncol(D)))+outer(g$b,g$k)); lmc<-log(pmax(D,0.5)/E)-BF(mu)
  W<-mu/mean(mu); am<-rowSums(W*lmc)/rowSums(W); Z<-(lmc-am)*sqrt(W)
  sv<-svd(Z); u1<-sv$u[,1]; v1<-sv$v[,1]; if(sum(u1)<0){u1<--u1;v1<--v1}
  bb<-u1/sum(u1); kk<-sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
  kk<-kk-mean(kk); sb<-sum(bb); list(a=am,b=bb/sb,k=kk*sb) }
SPEC <- list(
  list(id="標準 LC",            col="#B03A2E", pch=16, f=function(D,E) lc_svd_fit(log(pmax(D,0.5)/E))),
  list(id="Firth",              col="#D68910", pch=17, f=function(D,E) lc_poisson_firth(D,E,firth=TRUE)),
  list(id="逐格校正",           col="#7D3C98", pch=15, f=function(D,E) lc_analytic(D,E,bfun=BF,shrink=1)),
  list(id="逐格校正＋加權",     col="#1F6FB2", pch=18, f=fit_cellw),
  list(id="僅校正 α",           col="#117A65", pch=1,  f=function(D,E) lc_alpha_only(D,E,bfun=BF,weight=FALSE)),
  list(id="僅校正 α＋加權",     col="#1E8449", pch=19, f=function(D,E) lc_alpha_only(D,E,bfun=BF,weight=TRUE))
)
K <- length(SPEC)
AA <- BB <- vector("list", length(NS))
for (i in seq_along(NS)) {
  N<-NS[i]; E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mtrue
  aA<-array(NA_real_,c(REPS,A,K)); bA<-array(NA_real_,c(REPS,A,K))
  for (r in seq_len(REPS)) {
    set.seed(SEED+1000*i+r)
    D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    for (j in seq_len(K)) {
      f<-tryCatch(SPEC[[j]]$f(D,E),error=function(e) NULL)
      if(is.null(f)||!all(is.finite(f$a))||!all(is.finite(f$b))) next
      aA[r,,j]<-f$a-truth$a; bA[r,,j]<-f$b
    }
    if(r%%25==0) cat(sprintf("  N=%.0e rep %d/%d\n",N,r,REPS))
  }
  AA[[i]]<-aA; BB[[i]]<-bA
}
amax <- function(M) max(abs(apply(M,2,median,na.rm=TRUE)))
bsse <- function(M) median(apply(M,1,function(v) sum((v-truth$b)^2)/sum(truth$b^2)),na.rm=TRUE)

## ---------- figT6：Pareto 圖 ----------
png_cjk("output/figures/figT6_pareto.png", width=2250, height=820, res=150)
par(mfrow=c(1,3), mar=c(5.0,5.0,3.4,1.0), mgp=c(3.0,0.8,0))
for (i in seq_along(NS)) {
  ax <- sapply(seq_len(K), function(j) amax(AA[[i]][,,j]))
  bx <- sapply(seq_len(K), function(j) bsse(BB[[i]][,,j]))
  plot(bx, ax, log="xy", type="n",
       xlim=range(bx)*c(0.75,1.5), ylim=range(ax)*c(0.6,1.9),
       xlab=expression(paste("SSE(", beta, ")  →  愈小愈好")),
       ylab=expression(paste(alpha, " 最大偏誤  →  愈小愈好")),
       main=sprintf("N = %s", format(NS[i], big.mark=",", scientific=FALSE)))
  grid(col="grey92", lty=3)
  ## 由標準 LC 指向兩個改良方向的箭頭
  arrows(bx[1], ax[1], bx[1], ax[5], col="grey55", lwd=1.6, length=0.09, lty=2)
  arrows(bx[1], ax[1], bx[4], ax[1], col="grey55", lwd=1.6, length=0.09, lty=2)
  arrows(bx[1], ax[1], bx[6], ax[6], col="black",  lwd=2.2, length=0.11)
  for (j in seq_len(K)) points(bx[j], ax[j], col=SPEC[[j]]$col, pch=SPEC[[j]]$pch, cex=1.9, lwd=2)
  text(bx, ax, sapply(SPEC,`[[`,"id"), pos=c(4,4,2,3,4,1), cex=0.72, col=sapply(SPEC,`[[`,"col"))
  if (i==1) legend("topright", bty="n", cex=0.7,
    legend=c("只改善 α（縱向）","只改善 β（橫向）","兩者同時（對角）"),
    col=c("grey55","grey55","black"), lwd=c(1.6,1.6,2.2), lty=c(2,2,1))
}
dev.off()

## ---------- figT7：逐年齡帶 ----------
band <- function(M) list(lo=apply(M,2,quantile,0.05,na.rm=TRUE),
                         md=apply(M,2,median,na.rm=TRUE),
                         hi=apply(M,2,quantile,0.95,na.rm=TRUE))
tp <- function(col,a=0.20){v<-col2rgb(col)/255; grDevices::rgb(v[1],v[2],v[3],a)}
i <- 2   # N = 5e4
use <- c(1,2,6)
png_cjk("output/figures/figT7_joint_band.png", width=2250, height=880, res=150)
par(mfrow=c(1,2), mar=c(5.4,4.8,3.4,1.0), mgp=c(2.9,0.8,0))
## alpha
bd <- lapply(use, function(j) band(AA[[i]][,,j]))
allv <- unlist(lapply(bd,function(b) c(b$lo,b$hi)))
yl <- range(c(quantile(allv,c(0.03,0.97),na.rm=TRUE),0)); yl <- yl+c(-1,1)*0.12*diff(yl)
plot(NA, xlim=c(1,A), ylim=yl, xaxt="n", xlab="年齡組",
     ylab=expression(hat(alpha)[x]-alpha[x]), main="(a) α 的偏誤（N = 50,000）")
axis(1,at=xx,labels=ages,las=2,cex.axis=0.62); abline(h=0,col="grey35",lwd=1.2)
abline(v=xx,col="grey94",lty=3)
for (m in seq_along(use)) { j<-use[m]; b<-bd[[m]]; cl<-SPEC[[j]]$col
  polygon(c(xx,rev(xx)),c(pmax(b$lo,yl[1]),rev(pmin(b$hi,yl[2]))),col=tp(cl),border=NA)
  lines(xx,pmax(pmin(b$lo,yl[2]),yl[1]),col=cl,lwd=1.0,lty=2)
  lines(xx,pmax(pmin(b$hi,yl[2]),yl[1]),col=cl,lwd=1.0,lty=2)
  lines(xx,b$md,col=cl,lwd=2.4) }
legend("topright",bty="n",cex=0.76,legend=sapply(SPEC[use],`[[`,"id"),
       col=sapply(SPEC[use],`[[`,"col"),lwd=2.4)
## beta
bd <- lapply(use, function(j) band(BB[[i]][,,j]))
rowv <- unlist(lapply(bd,function(b) c(b$lo,b$hi)))
yl <- range(c(quantile(rowv,c(0.12,0.88),na.rm=TRUE),truth$b)); yl<-yl+c(-1,1)*0.12*diff(yl)
plot(NA, xlim=c(1,A), ylim=yl, xaxt="n", xlab="年齡組",
     ylab=expression(hat(beta)[x]), main="(b) β 的估計（N = 50,000）")
axis(1,at=xx,labels=ages,las=2,cex.axis=0.62); abline(h=0,col="grey35",lwd=1.2)
abline(v=xx,col="grey94",lty=3)
for (m in seq_along(use)) { j<-use[m]; b<-bd[[m]]; cl<-SPEC[[j]]$col
  polygon(c(xx,rev(xx)),c(pmax(b$lo,yl[1]),rev(pmin(b$hi,yl[2]))),col=tp(cl),border=NA)
  lines(xx,pmax(pmin(b$lo,yl[2]),yl[1]),col=cl,lwd=1.0,lty=2)
  lines(xx,pmax(pmin(b$hi,yl[2]),yl[1]),col=cl,lwd=1.0,lty=2)
  lines(xx,b$md,col=cl,lwd=2.2) }
lines(xx,truth$b,col="black",lwd=2.8); points(xx,truth$b,pch=16,cex=0.65)
legend("topleft",bty="n",cex=0.76,legend=c("真值",sapply(SPEC[use],`[[`,"id")),
       col=c("black",sapply(SPEC[use],`[[`,"col")),lwd=c(2.8,2.4,2.4,2.4))
dev.off()
cat("\n=== N=5e4 的 alpha/beta 同時比較 ===\n")
for (j in seq_len(K)) cat(sprintf("  %-18s a最大 %.4f   SSE(b) %.4f\n",
  SPEC[[j]]$id, amax(AA[[2]][,,j]), bsse(BB[[2]][,,j])))
cat("\n已輸出 figT6 / figT7\n")
