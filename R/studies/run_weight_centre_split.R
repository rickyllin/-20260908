###############################################################################
# 權重與中心化的分離，以及分離式的估計量設計
# （進度報告_理論初探_擬概似與變異安定轉換 第肆節）
#
#   命題 2 的分工預測：權重修正變異函數故作用於期望死亡數大的年齡（beta、
#   kappa、漂移）；中心化修正轉換的偏誤故作用於期望死亡數小的年齡（alpha）。
#   本程式以 2x2 檢驗，結果顯示分離是精確的而非近似的：
#     只中心化時 SSE(beta)、漂移、cor(kappa)、kappa 全距四項到小數第四位
#       與標準 Lee-Carter 相同——設計使然，中心化只在 SVD 之後把 bbar 扣在
#       alpha 上，依建構不觸及 beta 與 kappa。
#     只加權時 alpha 最大偏誤幾乎不動（0.8658 -> 0.8660）。
#
#   據此修改設計：加權只用於秩一結構，alpha 維持等權列平均（閉式
#   mean_t b 本來正是對等權列平均成立的）。五萬時 alpha 中位偏誤
#   0.0193 -> 0.0067、最大 0.0768 -> 0.0391，而 SSE(beta) 只差 1.5%。
#   一萬時幾乎全面較優（e0 -1.317 -> -0.649）。
#
#   五萬時 e0 由 -0.041 變為 +0.273 不是精度變差，而是失去一個相消：
#   依 run_bias_attribution.R，-0.041 是 alpha(-0.326)、beta(+1.492)、
#   kappa(-0.603) 三項相消的結果。
#
#   另檢驗一項推測而未成立：現行中心化項取 b 的算術平均而列平均是加權的，
#   兩者不一致；改為一致後差異可忽略（0.0193->0.0203、0.0442->0.0437），
#   因權重在同一年齡內只隨時間趨勢變動、列內近乎等權。alpha 因加權而變差
#   的真正來源應是權重由資料估出、與該列的 l 相關所生的比值偏誤，尚未推導。
#
# 輸出：終端三組表格（報告表 5、表 6）
###############################################################################
## ---- 一、2x2：權重與中心化分開施作 --------------------------------------
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R"); source("R/core/analytic_bias_lc.R")
SEED<-20261006; REPS<-100; c0<-0.5
dat<-load_data(sex="Female"); keep<-which(dat$years %in% 2001:2024)
D0<-dat$D[,keep,drop=FALSE]; E0<-dat$E[,keep,drop=FALSE]
A<-nrow(D0); Tn<-ncol(D0); ages<-rownames(D0)
tr<-lc_poisson_firth(D0,E0,firth=FALSE); a0<-tr$a;b0<-tr$b;k0<-tr$k
mt<-exp(outer(a0,rep(1,Tn))+outer(b0,k0)); wage<-E0[,Tn]/sum(E0[,Tn])
BF<-make_logbias(c0); NOB<-function(m) 0*m
e0t<-life_table(mt[,Tn])$e0
LAB<-c("w=1、不中心化（標準）","w=mu、不中心化（只加權）",
       "w=1、中心化（只中心化）","w=mu、中心化（兩者）")
cat("### 加權與中心化的 2x2：各自修正哪一個參數\n")
cat("    命題 2 的預測：加權修正變異函數故作用於 mu 大的年齡（漂移、kappa）；\n")
cat("    中心化修正轉換的偏誤故作用於 mu 小的年齡（alpha）。\n")
for (N in c(1e4,5e4)) {
  E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mt
  K<-4; aM<-array(NA_real_,c(REPS,A,K)); bS<-eS<-dr<-ck<-rg<-matrix(NA_real_,REPS,K)
  for (r in seq_len(REPS)) { set.seed(SEED+17000+r)
    D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    fits<-list(
      lc_alpha_only(D,E,bfun=NOB,weight=FALSE),
      lc_alpha_only(D,E,bfun=NOB,weight=TRUE),
      lc_alpha_only(D,E,bfun=BF, weight=FALSE),
      lc_alpha_only(D,E,bfun=BF, weight=TRUE))
    for (j in 1:K) { f<-fits[[j]]
      aM[r,,j]<-f$a-a0; bS[r,j]<-sum((f$b-b0)^2)/sum(b0^2)
      dr[r,j]<-(f$k[Tn]-f$k[1])/(Tn-1)-(k0[Tn]-k0[1])/(Tn-1)
      ck[r,j]<-cor(f$k,k0); rg[r,j]<-diff(range(f$k))
      eS[r,j]<-tryCatch(life_table(exp(f$a+f$b*f$k[Tn]))$e0,error=function(e)NA) } }
  cat(sprintf("\n===== N = %.0e（真值 kappa 全距 %.3f）=====\n", N, diff(range(k0))))
  cat(sprintf("  %-26s %8s %8s %8s %8s %8s %8s\n","估計量",
      "a 中位","a 最大","SSEb","漂移偏誤","cor(k)","k 全距"))
  for (j in 1:4) {
    da<-apply(aM[,,j],2,median)
    cat(sprintf("  %-26s %8.4f %8.4f %8.4f %+8.4f %8.3f %8.3f\n", LAB[j],
        median(abs(da)), max(abs(da)), median(bS[,j],na.rm=TRUE),
        median(dr[,j],na.rm=TRUE), median(ck[,j],na.rm=TRUE), median(rg[,j],na.rm=TRUE)))
  }
  cat("  e0 偏誤：")
  for (j in 1:4) cat(sprintf(" %s %+.3f｜", c("標準","只加權","只中心化","兩者")[j],
      median(eS[,j],na.rm=TRUE)-e0t))
  cat("\n")
}

## ---- 二、中心化項的權重是否須一致（結論：可忽略）------------------------
#' lc_alpha_only 的修正版：中心化項與列平均採同一組權重
lc_alpha_cons <- function(D,E,c0=0.5,weight=TRUE,maxit=20,tol=1e-9,bfun=NULL,
                          consistent=TRUE) {
  if (is.null(bfun)) bfun <- make_logbias(c0)
  A<-nrow(D); Tn<-ncol(D); lm_<-log(pmax(D,c0)/E)
  f<-lc_svd_fit(lm_); a<-f$a;b<-f$b;k<-f$k
  for (it in seq_len(maxit)) {
    mu<-E*exp(outer(a,rep(1,Tn))+outer(b,k)); Bm<-bfun(mu)
    if (weight) {
      W<-mu/mean(mu)
      am<-rowSums(W*lm_)/rowSums(W)
      bbar <- if (consistent) rowSums(W*Bm)/rowSums(W) else rowMeans(Bm)
      Z<-(lm_-am)*sqrt(W); sv<-svd(Z); u1<-sv$u[,1]; v1<-sv$v[,1]
      if(sum(u1)<0){u1<--u1;v1<--v1}
      b2<-u1/sum(u1); k2<-sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
      k2<-k2-mean(k2); sb<-sum(b2); b2<-b2/sb; k2<-k2*sb
    } else {
      f2<-lc_svd_fit(lm_); am<-f2$a;b2<-f2$b;k2<-f2$k; bbar<-rowMeans(Bm)
    }
    a2<-am-bbar
    d<-max(abs(a2-a),abs(b2-b),abs(k2-k)/max(1,max(abs(k)))); a<-a2;b<-b2;k<-k2
    if(d<tol) break }
  list(a=a,b=b,k=k) }
dat<-load_data(sex="Female"); keep<-which(dat$years %in% 2001:2024)
D0<-dat$D[,keep,drop=FALSE]; E0<-dat$E[,keep,drop=FALSE]
A<-nrow(D0); Tn<-ncol(D0)
tr<-lc_poisson_firth(D0,E0,firth=FALSE); a0<-tr$a;b0<-tr$b;k0<-tr$k
mt<-exp(outer(a0,rep(1,Tn))+outer(b0,k0)); wage<-E0[,Tn]/sum(E0[,Tn])
BF<-make_logbias(c0); e0t<-life_table(mt[,Tn])$e0
LAB<-c("只中心化（w=1，一致）","加權＋算術平均的 b（現行）",
       "加權＋加權平均的 b（修正）")
cat("### 中心化項的權重是否須與列平均一致\n")
cat("    E[alpha_hat] = sum_t w (log m + b) / sum_t w，故應扣加權平均的 b。\n")
for (N in c(1e4,5e4)) {
  E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mt
  K<-3; aM<-array(NA_real_,c(REPS,A,K)); bS<-eS<-matrix(NA_real_,REPS,K)
  for (r in seq_len(REPS)) { set.seed(SEED+17000+r)
    D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    fits<-list(lc_alpha_cons(D,E,bfun=BF,weight=FALSE),
               lc_alpha_cons(D,E,bfun=BF,weight=TRUE,consistent=FALSE),
               lc_alpha_cons(D,E,bfun=BF,weight=TRUE,consistent=TRUE))
    for (j in 1:K){f<-fits[[j]]; aM[r,,j]<-f$a-a0
      bS[r,j]<-sum((f$b-b0)^2)/sum(b0^2)
      eS[r,j]<-tryCatch(life_table(exp(f$a+f$b*f$k[Tn]))$e0,error=function(e)NA)} }
  cat(sprintf("\n  N=%.0e\n  %-30s %9s %9s %9s %9s\n", N,
      "估計量","a 中位","a 最大","SSEb","e0 偏誤"))
  for (j in 1:K) { da<-apply(aM[,,j],2,median)
    cat(sprintf("  %-30s %9.4f %9.4f %9.4f %+9.3f\n", LAB[j],
        median(abs(da)), max(abs(da)), median(bS[,j],na.rm=TRUE),
        median(eS[,j],na.rm=TRUE)-e0t)) }
}

## ---- 三、分離式設計 ------------------------------------------------------
#' alpha 用等權列平均（閉式所適用者），beta 與 kappa 用加權奇異值分解
lc_split <- function(D,E,bfun,maxit=20,tol=1e-9) {
  A<-nrow(D); Tn<-ncol(D); lm_<-log(pmax(D,c0)/E)
  f<-lc_svd_fit(lm_); a<-f$a;b<-f$b;k<-f$k
  for (it in seq_len(maxit)) {
    mu<-E*exp(outer(a,rep(1,Tn))+outer(b,k))
    am <- rowMeans(lm_)                       # 等權：閉式正是對此成立
    bbar <- rowMeans(bfun(mu))
    W  <- mu/mean(mu)                         # 加權只用於秩一結構
    Z  <- (lm_-am)*sqrt(W)
    sv<-svd(Z); u1<-sv$u[,1]; v1<-sv$v[,1]
    if(sum(u1)<0){u1<--u1;v1<--v1}
    b2<-u1/sum(u1); k2<-sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
    k2<-k2-mean(k2); sb<-sum(b2); b2<-b2/sb; k2<-k2*sb
    a2<-am-bbar
    d<-max(abs(a2-a),abs(b2-b),abs(k2-k)/max(1,max(abs(k)))); a<-a2;b<-b2;k<-k2
    if(d<tol) break }
  list(a=a,b=b,k=k) }
dat<-load_data(sex="Female"); keep<-which(dat$years %in% 2001:2024)
D0<-dat$D[,keep,drop=FALSE]; E0<-dat$E[,keep,drop=FALSE]
A<-nrow(D0); Tn<-ncol(D0)
tr<-lc_poisson_firth(D0,E0,firth=FALSE); a0<-tr$a;b0<-tr$b;k0<-tr$k
mt<-exp(outer(a0,rep(1,Tn))+outer(b0,k0)); wage<-E0[,Tn]/sum(E0[,Tn])
BF<-make_logbias(c0); e0t<-life_table(mt[,Tn])$e0
LAB<-c("標準 Lee-Carter","只中心化","只加權","現行（兩者同時）",
       "分離式：alpha 等權、beta/kappa 加權")
cat("### 分離式：把加權只用在秩一結構，alpha 維持等權列平均\n")
for (N in c(1e4,5e4)) {
  E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mt
  K<-5; aM<-array(NA_real_,c(REPS,A,K)); bS<-eS<-dr<-ck<-matrix(NA_real_,REPS,K)
  for (r in seq_len(REPS)) { set.seed(SEED+17000+r)
    D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    NOB<-function(m) 0*m
    fits<-list(lc_alpha_only(D,E,bfun=NOB,weight=FALSE),
               lc_alpha_only(D,E,bfun=BF, weight=FALSE),
               lc_alpha_only(D,E,bfun=NOB,weight=TRUE),
               lc_alpha_only(D,E,bfun=BF, weight=TRUE),
               lc_split(D,E,bfun=BF))
    for (j in 1:K){f<-fits[[j]]; aM[r,,j]<-f$a-a0
      bS[r,j]<-sum((f$b-b0)^2)/sum(b0^2); ck[r,j]<-cor(f$k,k0)
      dr[r,j]<-(f$k[Tn]-f$k[1])/(Tn-1)-(k0[Tn]-k0[1])/(Tn-1)
      eS[r,j]<-tryCatch(life_table(exp(f$a+f$b*f$k[Tn]))$e0,error=function(e)NA)} }
  cat(sprintf("\n  N=%.0e\n  %-36s %8s %8s %8s %8s %8s %8s\n", N,
      "估計量","a 中位","a 最大","SSEb","漂移","cor(k)","e0 偏誤"))
  for (j in 1:K) { da<-apply(aM[,,j],2,median)
    cat(sprintf("  %-36s %8.4f %8.4f %8.4f %+8.4f %8.3f %+8.3f\n", LAB[j],
        median(abs(da)), max(abs(da)), median(bS[,j],na.rm=TRUE),
        median(dr[,j],na.rm=TRUE), median(ck[,j],na.rm=TRUE),
        median(eS[,j],na.rm=TRUE)-e0t)) }
}
