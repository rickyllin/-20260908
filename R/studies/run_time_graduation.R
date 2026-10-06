###############################################################################
# 沿時間方向的修勻（分支：修勻與貝氏混合｜進度報告_1012 第參節）
#
#   命題（報告第參節）：對任一凸的 J，loss (r-y)'W(r-y) + J(Delta^z r) 的解
#   滿足 1'W r = 1'W y（z >= 1），且 z >= 2 時另有 t'W r = t'W y。
#   故 alpha_hat（即 W 加權列平均）精確不變，沿時間修勻只作用於 beta 與
#   kappa 的形狀，與中心化校正完全不互相干擾——這與年齡方向須輸送校正量
#   （run_graduation_transport.R）形成對比。
#
#   本程式驗證該不變性，並量出修勻對 SSE(beta)、漂移與 kappa 粗糙度的效果。
#   kappa 粗糙度以 sum_t (Delta^2 kappa_t)^2 衡量，真值為 20.57；
#   lambda 應校到粗糙度對上真值（一萬約 0.3、五萬約 0.1）。
#
# 輸出：終端表格（報告表 3）
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
SEED <- 20261006; REPS <- 100
dat <- load_data(sex="Female"); keep <- which(dat$years %in% 2001:2024)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]
A <- nrow(D0); Tn <- ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage <- E0[,Tn]/sum(E0[,Tn]); BF <- make_logbias(0.5)
tt <- seq_len(Tn)
drift_end <- function(k) (k[Tn]-k[1])/(Tn-1)                  # 本文現用
drift_wls <- function(k,w) sum(w*(tt-sum(w*tt)/sum(w))*k)/sum(w*(tt-sum(w*tt)/sum(w))^2)
d_end_true <- drift_end(truth$k); d_wls_true <- drift_wls(truth$k, rep(1,Tn))

#' 沿時間的 Whittaker 修勻，逐年齡施作，權重取該年齡各年的 mu_hat
time_smooth <- function(L, W, lam, z = 2) {
  Dm <- diag(ncol(L)); for(i in seq_len(z)) Dm <- diff(Dm); P <- crossprod(Dm)
  out <- L
  for (x in seq_len(nrow(L))) {
    w <- W[x,]; out[x,] <- solve(diag(w)+lam*mean(w)*P, w*L[x,])
  }
  out
}
lc_time <- function(D,E,bfun,lam,z=2,maxit=20,tol=1e-9){
  lm_<-log(pmax(D,0.5)/E); f<-lc_svd_fit(lm_); a<-f$a;b<-f$b;k<-f$k
  for(it in seq_len(maxit)){
    mu<-E*exp(outer(a,rep(1,Tn))+outer(b,k)); W<-mu/mean(mu)
    lt <- if (is.null(lam)) lm_ else time_smooth(lm_, W, lam, z)
    bb<-rowMeans(bfun(mu))
    am<-rowSums(W*lt)/rowSums(W); Z<-(lt-am)*sqrt(W)
    sv<-svd(Z);u1<-sv$u[,1];v1<-sv$v[,1]; if(sum(u1)<0){u1<--u1;v1<--v1}
    b2<-u1/sum(u1); k2<-sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
    k2<-k2-mean(k2); sb<-sum(b2); b2<-b2/sb; k2<-k2*sb; a2<-am-bb
    d<-max(abs(a2-a),abs(b2-b),abs(k2-k)/max(1,max(abs(k)))); a<-a2;b<-b2;k<-k2
    if(d<tol) break}
  list(a=a,b=b,k=k,W=W)}

LAM<-c(NA,0.1,0.4,10)
for (N in c(1e4,5e4)) {
  E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mtrue
  out<-array(NA_real_,c(REPS,length(LAM),5)); rough<-matrix(NA_real_,REPS,length(LAM))
  for(r in seq_len(REPS)){
    set.seed(SEED+13000+r)
    D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    for(j in seq_along(LAM)){
      f<-tryCatch(lc_time(D,E,BF,if(is.na(LAM[j])) NULL else LAM[j]),error=function(e)NULL)
      if(is.null(f)) next
      out[r,j,1]<-max(abs(f$a-truth$a))
      out[r,j,2]<-median(abs(f$a-truth$a))
      out[r,j,3]<-sum((f$b-truth$b)^2)/sum(truth$b^2)
      out[r,j,4]<-drift_end(f$k)-d_end_true
      out[r,j,5]<-drift_wls(f$k,rep(1,Tn))-d_wls_true
      rough[r,j]<-sum(diff(diff(f$k))^2)
    }}
  cat(sprintf("\n===== N = %.0e （真值 二階差分平方和 %.4f）=====\n", N, sum(diff(diff(truth$k))^2)))
  cat(sprintf("%-10s %10s %10s %10s %12s %12s %10s\n","λ","α最大","α中位","SSEβ","漂移(端點)","漂移(WLS)","κ粗糙度"))
  for(j in seq_along(LAM))
    cat(sprintf("%-10s %10.4f %10.4f %10.4f %12.4f %12.4f %10.4f\n",
      ifelse(is.na(LAM[j]),"不修勻",LAM[j]),
      median(out[,j,1],na.rm=TRUE), median(out[,j,2],na.rm=TRUE),
      median(out[,j,3],na.rm=TRUE), median(out[,j,4],na.rm=TRUE),
      median(out[,j,5],na.rm=TRUE), median(rough[,j],na.rm=TRUE)))
}
