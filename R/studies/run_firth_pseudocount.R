###############################################################################
# 核對：Firth pseudo-count 的雙線性項在數值上是否要緊
# （1008策略文件_可行性評估.md 第三節；對應該文件一 §8.1 的方法三）
#
#   Lee-Carter 是雙線性而非線性預測子，故 Firth 的 pseudo-count 分解為
#       a_{x,t} = h_{x,t}/2 + c_{x,t}，  c_{x,t} = mu * (I_c^-)_{beta_x, kappa_t}，
#   而第二項可為負。現行實作（fh_firth_kalman.R 的 lc_leverage/2）只做第一項。
#
#   兩項量測。其一，|c|/(h/2) 的中位數為 0.006 至 0.011，但最大值在人數一萬
#   時達 1.398 且落在 1-4 歲，且有 249/528 格的 c < 0。其二，實際改動估計量
#   的效果：人數一萬時 alpha 最大偏誤 0.2993 -> 0.3161、SSE(beta)
#   18.105 -> 16.904（方向相反），人數五萬時三者到小數四位幾乎相同。
#
#   結論：分解在理論上正確、在數值上不要緊，建議只在報告中交代而不改實作。
#   又 D + h/2 + c 在 300 組樣本中從未為負，故文件所提的 cellwise damping
#   在本文的資料上從未啟動。
#
#   另順帶查出 sum(h/2) = 33.000 而非 p/2 = 34——因參數化有兩個冗餘方向、
#   rank(I) = p - 2。1012 V1 原寫「總量恆為 p/2」須更正。
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
library(MASS)
dat<-load_data(sex="Female"); keep<-which(dat$years %in% 2001:2024)
D0<-dat$D[,keep,drop=FALSE]; E0<-dat$E[,keep,drop=FALSE]
A<-nrow(D0); Tn<-ncol(D0)
tr<-lc_poisson_firth(D0,E0,firth=FALSE); a0<-tr$a;b0<-tr$b;k0<-tr$k
mt<-exp(outer(a0,rep(1,Tn))+outer(b0,k0)); wage<-E0[,Tn]/sum(E0[,Tn])
## 完整的 Firth pseudo-count：a = h/2 + c，c_{x,t} = mu * (I^-)[beta_x, kappa_t]
parts <- function(a,b,k,E) {
  n<-A*Tn; p<-2*A+Tn
  mu<-as.vector(t(E*exp(outer(a,rep(1,Tn))+outer(b,k))))   # 列優先，與 lc_leverage 一致
  xi<-rep(seq_len(A),each=Tn); ti<-rep(seq_len(Tn),times=A)
  J<-matrix(0,n,p)
  J[cbind(seq_len(n),xi)]<-1
  J[cbind(seq_len(n),A+xi)]<-k[ti]
  J[cbind(seq_len(n),2*A+ti)]<-b[xi]
  Jw<-J*sqrt(mu); Ii<-ginv(crossprod(Jw))
  h<-rowSums((Jw%*%Ii)*Jw)
  cc<-mu*Ii[cbind(A+xi, 2*A+ti)]       # 雙線性項（beta_x, kappa_t 交叉區塊）
  list(h=matrix(h,A,Tn,byrow=TRUE), c=matrix(cc,A,Tn,byrow=TRUE))
}
cat("### Firth pseudo-count 的兩項：h/2（現行實作）與 c（雙線性項，現行實作所無）\n\n")
for (N in c(1e4, 5e4, 2e5)) {
  E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mt
  set.seed(20260926+11001); D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
  f<-lc_poisson_firth(D,E,firth=TRUE)
  P<-parts(f$a,f$b,f$k,E)
  h2<-P$h/2; cc<-P$c
  cat(sprintf("  N=%.0e｜sum(h/2)=%.3f（應為 p/2=%.1f）｜sum(c)=%+.4f\n",
      N, sum(h2), (2*A+Tn)/2, sum(cc)))
  cat(sprintf("        |c|/|h/2| 的中位 %.4f、最大 %.4f｜c<0 的格數 %d/%d\n",
      median(abs(cc)/h2), max(abs(cc)/h2), sum(cc<0), length(cc)))
  cat(sprintf("        D + h/2 + c 為負的格數 %d｜其中 D=0 的 %d\n",
      sum(D + h2 + cc < 0), sum(D==0 & (D + h2 + cc) < 0)))
  o<-order(-abs(cc))[1:3]
  cat(sprintf("        |c| 最大的三格：%s\n", paste(sprintf("%s/%s: h/2=%.3f c=%+.3f",
      rownames(D0)[(o-1)%%A+1], colnames(D0)[(o-1)%/%A+1], h2[o], cc[o]), collapse="； ")))
}

## ---- 改動估計量的效果 ----------------------------------------------------
SEED<-20260926; REPS<-100
dat<-load_data(sex="Female"); keep<-which(dat$years %in% 2001:2024)
D0<-dat$D[,keep,drop=FALSE]; E0<-dat$E[,keep,drop=FALSE]
A<-nrow(D0); Tn<-ncol(D0)
tr<-lc_poisson_firth(D0,E0,firth=FALSE); a0<-tr$a;b0<-tr$b;k0<-tr$k
mt<-exp(outer(a0,rep(1,Tn))+outer(b0,k0)); wage<-E0[,Tn]/sum(E0[,Tn])
e0t<-life_table(mt[,Tn])$e0
pc <- function(a,b,k,E,mode) {           # pseudo-count 的三種取法
  n<-A*Tn; p<-2*A+Tn
  mu<-as.vector(t(E*exp(outer(a,rep(1,Tn))+outer(b,k))))
  xi<-rep(seq_len(A),each=Tn); ti<-rep(seq_len(Tn),times=A)
  J<-matrix(0,n,p); J[cbind(seq_len(n),xi)]<-1
  J[cbind(seq_len(n),A+xi)]<-k[ti]; J[cbind(seq_len(n),2*A+ti)]<-b[xi]
  Jw<-J*sqrt(mu); Ii<-ginv(crossprod(Jw))
  h<-rowSums((Jw%*%Ii)*Jw)/2
  cc<-mu*Ii[cbind(A+xi,2*A+ti)]
  v <- switch(mode, "h"=h, "full"=h+cc, "damp"=h+pmax(cc,-h))  # damp：保 pseudo-count 非負
  matrix(v,A,Tn,byrow=TRUE)
}
fit <- function(D,E,mode,maxit=400,tol=1e-10) {
  a<-log(pmax(rowSums(D),0.5)/rowSums(E)); b<-rep(1/A,A); k<-seq(1,-1,length.out=Tn)
  for (it in seq_len(maxit)) {
    h<-pc(a,b,k,E,mode)
    mu<-E*exp(outer(a,rep(1,Tn))+outer(b,k))
    a<-a+rowSums(D-mu+h)/pmax(rowSums(mu),1e-12)
    mu<-E*exp(outer(a,rep(1,Tn))+outer(b,k))
    k<-k+colSums((D-mu+h)*b)/pmax(colSums(mu*b^2),1e-12); k<-k-mean(k)
    mu<-E*exp(outer(a,rep(1,Tn))+outer(b,k))
    bn<-b+rowSums((D-mu+h)*matrix(k,A,Tn,byrow=TRUE))/
          pmax(rowSums(mu*matrix(k^2,A,Tn,byrow=TRUE)),1e-12)
    s<-sum(bn); if(abs(s)<1e-12) break
    nb<-bn/s; nk<-k*s
    del<-max(max(abs(nb-b)),max(abs(nk-k))/max(1,max(abs(k)))); b<-nb;k<-nk
    if(del<tol) break }
  list(a=a,b=b,k=k) }
LAB<-c("h/2（現行實作）","h/2 + c（完整 Firth）","h/2 + max(c,-h/2)（阻尼）")
cat("### Firth 的 pseudo-count：省略雙線性項是否要緊\n")
for (N in c(1e4,5e4)) {
  E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mt
  K<-3; aM<-array(NA_real_,c(REPS,A,K)); eS<-bS<-matrix(NA_real_,REPS,K)
  for (r in seq_len(REPS)) { set.seed(SEED+11000+r)
    D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    for (j in 1:K) { f<-tryCatch(fit(D,E,c("h","full","damp")[j]),error=function(e)NULL)
      if(is.null(f)||!all(is.finite(f$a))) next
      aM[r,,j]<-f$a-a0; bS[r,j]<-sum((f$b-b0)^2)/sum(b0^2)
      eS[r,j]<-tryCatch(life_table(exp(f$a+f$b*f$k[Tn]))$e0,error=function(e)NA) } }
  cat(sprintf("\n  N=%.0e\n  %-28s %9s %9s %9s %9s\n",N,"pseudo-count","a中位","a最大","SSEb","e0偏誤"))
  for (j in 1:K) { da<-apply(aM[,,j],2,median,na.rm=TRUE)
    cat(sprintf("  %-28s %9.4f %9.4f %9.4f %+9.3f\n", LAB[j],
      median(abs(da)),max(abs(da)),median(bS[,j],na.rm=TRUE),
      median(eS[,j],na.rm=TRUE)-e0t)) }
}
