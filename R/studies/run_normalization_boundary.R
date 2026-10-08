###############################################################################
# 核對：改用 ||beta||_2 = 1 可消掉該邊界
# （1008策略文件_可行性評估.md 第二節；對應該文件二 §10.2 的【猜想】）
#
#   該猜想未經證明，而它是四份文件裡報酬最高的一項——若成立則是一行改動。
#
#   本程式把同一條退化路徑上的點改以 ||beta||_2 = 1 表示
#   （c = ||beta||_2、bhat = beta/c、khat = c*kappa，eta 不變），
#   再比較兩種正規化下的 log det I。
#
#   結果：sum(beta)=1 下 ||beta|| 由 0.34 漲到 26050 而 log det I 發散
#   （每十倍 +4.606 = 2(T-X-1) log 10）；||beta||_2=1 下 ||kappa|| 收斂到
#   1.0213、log det I 收斂到 161.617。亦即那條逃逸路徑在後一種正規化下
#   只是一條收斂路徑，猜想成立。
#
#   理由是幾何的：sum(beta)=1 不限制 ||beta||（正負相消即可），而
#   ||beta||_2=1 把 beta 釘在緊緻的單位球上，故「beta 放大而 kappa 縮小」
#   這條使 eta 保持有界的出路被封住。
#
#   代價：兩種正規化描述同一個模型，故 eta、mu、m、e0 與所有可識別量不變，
#   改變的只有參數層的報告尺度；但 beta_hat 的標準誤會隨正規化而變
#   （1012 V1 第貳節之五所引 Cui and Xu 的論點）。
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
## 同一條退化路徑，在兩種正規化下的座標與 log det I
onb0 <- function(n){M<-diag(n)-1/n; qr.Q(qr(M))[,1:(n-1),drop=FALSE]}
onbS <- function(b){            # 單位球在 b 處的切空間（與 b 正交），維度 X-1
  X<-length(b); Q<-qr.Q(qr(cbind(b, diag(X))))[,2:X,drop=FALSE]; Q }
logdetI <- function(a,b,k,E,Bb,Bk){
  X<-length(a); T<-length(k); n<-X*T
  mu<-as.vector(E*exp(outer(a,rep(1,T))+outer(b,k)))
  xi<-rep(seq_len(X),times=T); ti<-rep(seq_len(T),each=X)
  q<-X+ncol(Bb)+ncol(Bk); J<-matrix(0,n,q)
  for(j in seq_len(X)) J[,j]<-as.numeric(xi==j)
  for(j in seq_len(ncol(Bb))) J[,X+j]<-Bb[xi,j]*k[ti]
  for(j in seq_len(ncol(Bk))) J[,X+ncol(Bb)+j]<-b[xi]*Bk[ti,j]
  sv<-svd(J*sqrt(mu))$d; 2*sum(log(sv)) }
set.seed(7); X<-22; T<-24
a<-seq(-9,-1,length.out=X); bt<-rnorm(X); bt<-bt-mean(bt); bt<-bt/sum(abs(bt))
kt<-rnorm(T); kt<-kt-mean(kt); E<-matrix(1e4/X,X,T)
Bk<-onb0(T)
cat("### 同一退化路徑在兩種正規化下的行為（X=22、T=24）\n\n")
cat(sprintf("%8s %12s %12s | %12s %12s %12s\n","s",
  "sum=1: |b|","sum=1: logdetI","L2=1: |k|","L2=1: logdetI","eta 最大變動"))
e_ref <- outer(a,rep(1,T))+outer(bt,kt)
for (s in 10^(0:5)) {
  b1<-s*bt+1/X; k1<-kt/s                       # sum(beta)=1 的座標
  ld1<-logdetI(a,b1,k1,E,onb0(X),Bk)
  cc<-sqrt(sum(b1^2)); b2<-b1/cc; k2<-k1*cc    # 同一點改以 ||beta||_2=1 表示
  ld2<-logdetI(a,b2,k2,E,onbS(b2),Bk)
  eta<-outer(a,rep(1,T))+outer(b1,k1)
  cat(sprintf("%8.0e %12.3f %12.4f | %12.4f %12.4f %12.2e\n",
      s, sqrt(sum(b1^2)), ld1, sqrt(sum(k2^2)), ld2, max(abs(eta-e_ref))))
}
cat("\n  （兩種正規化描述的是同一個模型點，故 eta 相同；差別只在座標與資訊矩陣。）\n")
