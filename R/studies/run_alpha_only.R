###############################################################################
# alpha 與 beta 能否同時改善？
#   比較三種設計：逐格校正、只校正 alpha、只校正 alpha + 加權。
# 輸出：output/tables/tableY_alpha_only.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
SEED <- 20260926; REPS <- 100; NS <- c(1e4,5e4,2e5); YEARS <- 2001:2024
dat <- load_data(sex="Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]; A<-nrow(D0); Tn<-ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage <- E0[,Tn]/sum(E0[,Tn]); dtrue<-(truth$k[Tn]-truth$k[1])/(Tn-1)
e0_true <- life_table(mtrue[,Tn])$e0; BF <- make_logbias(0.5)
fit_cellw <- function(D,E){ g<-lc_analytic(D,E,bfun=BF,shrink=1)
  mu<-E*exp(outer(g$a,rep(1,ncol(D)))+outer(g$b,g$k)); lmc<-log(pmax(D,0.5)/E)-BF(mu)
  W<-mu/mean(mu); am<-rowSums(W*lmc)/rowSums(W); Z<-(lmc-am)*sqrt(W)
  sv<-svd(Z); u1<-sv$u[,1]; v1<-sv$v[,1]; if(sum(u1)<0){u1<--u1;v1<--v1}
  bb<-u1/sum(u1); kk<-sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
  kk<-kk-mean(kk); sb<-sum(bb); list(a=am,b=bb/sb,k=kk*sb) }
FIT <- list("標準 LC"=function(D,E) lc_svd_fit(log(pmax(D,0.5)/E)),
            "Firth"=function(D,E) lc_poisson_firth(D,E,firth=TRUE),
            "逐格校正"=function(D,E) lc_analytic(D,E,bfun=BF,shrink=1),
            "逐格校正＋加權"=fit_cellw,
            "僅校正 alpha"=function(D,E) lc_alpha_only(D,E,bfun=BF,weight=FALSE),
            "僅校正 alpha＋加權"=function(D,E) lc_alpha_only(D,E,bfun=BF,weight=TRUE))
K<-length(FIT); res<-list()
for (N in NS) {
  E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mtrue
  aM<-array(NA_real_,c(REPS,A,K)); bS<-dS<-eS<-matrix(NA_real_,REPS,K)
  for (r in seq_len(REPS)) {
    set.seed(SEED+1000*which(NS==N)+r)
    D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    for (j in seq_len(K)) {
      f<-tryCatch(FIT[[j]](D,E),error=function(e) NULL)
      if(is.null(f)||!all(is.finite(f$a))||!all(is.finite(f$b))) next
      aM[r,,j]<-f$a-truth$a; bS[r,j]<-sum((f$b-truth$b)^2)/sum(truth$b^2)
      dS[r,j]<-(f$k[Tn]-f$k[1])/(Tn-1)-dtrue
      mh<-exp(outer(f$a,rep(1,Tn))+outer(f$b,f$k))
      eS[r,j]<-tryCatch(life_table(mh[,Tn])$e0,error=function(e) NA_real_)
    }
    if(r%%25==0) cat(sprintf("  N=%.0e rep %d/%d\n",N,r,REPS))
  }
  res[[as.character(N)]]<-data.frame(N=N,估計量=names(FIT),
    alpha中位=round(sapply(seq_len(K),function(j) median(abs(apply(aM[,,j],2,median,na.rm=TRUE)))),4),
    alpha最大=round(sapply(seq_len(K),function(j) max(abs(apply(aM[,,j],2,median,na.rm=TRUE)))),4),
    SSE_beta=round(apply(bS,2,median,na.rm=TRUE),4),
    漂移偏誤=round(apply(dS,2,median,na.rm=TRUE),4),
    e0偏誤=round(apply(eS,2,median,na.rm=TRUE)-e0_true,3),stringsAsFactors=FALSE)
  cat(sprintf("\n===== N = %.0e =====\n",N)); print(res[[as.character(N)]][,-1],row.names=FALSE)
}
tab<-do.call(rbind,res); write.csv(tab,"output/tables/tableY_alpha_only.csv",row.names=FALSE)
cat("\n已輸出 tableY\n")
