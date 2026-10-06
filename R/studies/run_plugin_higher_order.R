###############################################################################
# 插入式誤差的二項展開（進度報告_1012 第伍節之一）
#
#   本文扣除 b(mu_hat) 而非 b(mu)，其差可展開為
#       E[b(mu_hat)] - b(mu) = b'(mu)(E[mu_hat]-mu) + b''(mu)Var(mu_hat)/2 + ...
#   兩個導數皆由差分恆等式 d/dmu E[g(D)] = E[g(D+1)] - E[g(D)] 精確算出，
#   施用兩次得 b''(mu) = E[g(D+2)] - 2E[g(D+1)] + E[g(D)] + 1/mu^2。
#
#   結論分兩層：展開抓到了正確的結構（逐年齡相關 0.817 與 0.899），
#   但兩項在幼年組幾乎互相抵銷（人數五萬、1-4 歲：-0.1132 與 +0.0754，
#   兩項和 -0.0378，而實際只有 +0.0006），淨值落在三階，
#   故直接施以二階校正在數值上不穩定。這同時解釋為何只修勻 mu_hat 的
#   版本效果不如預期：修勻壓低二階項卻給 mu_hat 添上平滑偏誤而放大一階項，
#   而 |b'| 在小 mu 處甚大（mu=0.1 時 -9.31），後者佔上風。
#
# 輸出：終端表格（報告表 6）
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
source("R/core/analytic_bias_lc.R")
c0<-0.5; REPS<-200; SEED<-20261006
g<-function(d) log(pmax(d,c0))
Eg<-function(mu,s=0) sapply(mu,function(m){K<-max(60,ceiling(m+12*sqrt(m)+12));d<-0:K
  p<-dpois(d,m); sum(p*g(d+s))/sum(p)})
bfun <-function(mu) Eg(mu,0)-log(mu)
b1fun<-function(mu) Eg(mu,1)-Eg(mu,0)-1/mu
b2fun<-function(mu) Eg(mu,2)-2*Eg(mu,1)+Eg(mu,0)+1/mu^2
dat<-load_data(sex="Female"); keep<-which(dat$years %in% 2001:2024)
D0<-dat$D[,keep,drop=FALSE]; E0<-dat$E[,keep,drop=FALSE]
A<-nrow(D0); Tn<-ncol(D0)
truth<-lc_poisson_firth(D0,E0,firth=FALSE)
mtrue<-exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage<-E0[,Tn]/sum(E0[,Tn]); BF<-make_logbias(c0)
for (N in c(5e4, 2e5)) {
  E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mtrue
  MU<-array(NA_real_,c(REPS,A,Tn))
  for(r in seq_len(REPS)){ set.seed(SEED+13000+r)
    D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    f<-lc_alpha_only(D,E,bfun=BF,weight=TRUE)
    MU[r,,]<-E*exp(outer(f$a,rep(1,Tn))+outer(f$b,f$k)) }
  mhat<-apply(MU,c(2,3),mean); vhat<-apply(MU,c(2,3),var)
  bt<-matrix(bfun(as.vector(mu)),A); d1<-matrix(b1fun(as.vector(mu)),A)
  d2<-matrix(b2fun(as.vector(mu)),A)
  bias_mu<-mhat-mu
  Ebhat<-apply(MU,c(2,3),function(v) mean(bfun(v)))
  t1<-d1*bias_mu; t2<-0.5*d2*vhat
  cat(sprintf("\n===== N=%.0e（重複 %d 次）=====\n", N, REPS))
  cat(sprintf("%-8s %9s %9s %9s %9s %9s %9s\n","年齡","mu","E[b(muh)]-b(mu)","一階項","二階項","兩項和","mu偏誤%%"))
  for(x in c(1,2,3,4,5,8,12,16,20,22)) cat(sprintf("%-8s %9.2f %12.4f %12.4f %9.4f %9.4f %8.1f\n",
    rownames(D0)[x], mean(mu[x,]), mean(Ebhat[x,]-bt[x,]), mean(t1[x,]), mean(t2[x,]),
    mean(t1[x,]+t2[x,]), 100*mean(bias_mu[x,]/mu[x,])))
  obs<-rowMeans(Ebhat-bt); pred<-rowMeans(t1+t2)
  cat(sprintf("  逐年齡相關 %.4f｜迴歸斜率 %.3f｜一階項佔 %.0f%%、二階項佔 %.0f%%（絕對值和）\n",
    cor(obs,pred), coef(lm(obs~pred))[2],
    100*sum(abs(rowMeans(t1)))/sum(abs(rowMeans(t1))+abs(rowMeans(t2))),
    100*sum(abs(rowMeans(t2)))/sum(abs(rowMeans(t1))+abs(rowMeans(t2)))))
  cat(sprintf("  過度扣除總量（列平均後）：中位 %+.4f、最大 %+.4f\n",
    median(obs), obs[which.max(abs(obs))]))
}
