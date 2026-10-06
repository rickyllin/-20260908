###############################################################################
# 計數尺度各法與對數尺度各法的對照（理論初探_擬概似與變異安定轉換 第伍、陸節）
#
#   兩件待辦：以 Anscombe 尺度配適，以及補上有界影響擬概似的對照。
#   兩者都不取對數，故依命題 1 其估計方程無偏；問題在解的存在性與效率。
#
#   Anscombe 分三個版本，以分離三項實作選擇各自的貢獻：
#     近似平均 2 sqrt(mu+3/8) 對精確平均 E[h(D)]；
#     單位權重對精確變異 Var[h(D)] 的倒數。
#   變異安定在小 mu 處失效（mu=0.05 時精確變異 0.061 而名目為 1），
#   故單位權重在最需要資訊的年齡低估資訊量。
#
#   有界影響依 Cantoni and Ronchetti (2001)：以 Huber 截斷 Pearson 殘差，
#   並扣 a(mu) = E[psi_c(r)] 以保 Fisher 一致性。該修正項以截斷級數精確算出
#   （mu=0.05 時為 -0.147，不可略）。
#
#   解法器已與既有的 lc_poisson_firth 核對：取 u = D - mu 時兩者的
#   alpha 最大偏誤與 SSE(beta) 到小數第四位相同。
#
# 輸出：output/tables/tableQ1_count_scale.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/count_scale_lc.R")
SEED <- 20261006; REPS <- 100; NS <- c(1e4, 5e4); c0 <- 0.5
dat <- load_data(sex="Female"); keep <- which(dat$years %in% 2001:2024)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]
A <- nrow(D0); Tn <- ncol(D0)
tr <- lc_poisson_firth(D0,E0,firth=FALSE); a0<-tr$a; b0<-tr$b; k0<-tr$k
mt <- exp(outer(a0,rep(1,Tn))+outer(b0,k0)); wage <- E0[,Tn]/sum(E0[,Tn])
BF <- make_logbias(c0); e0t <- life_table(mt[,Tn])$e0

## alpha 等權、beta 與 kappa 加權的分離式（見第肆節之三）
lc_split <- function(D,E,bfun,maxit=20,tol=1e-9) {
  A<-nrow(D); Tn<-ncol(D); lm_<-log(pmax(D,c0)/E)
  f<-lc_svd_fit(lm_); a<-f$a;b<-f$b;k<-f$k
  for (it in seq_len(maxit)) {
    mu<-E*exp(outer(a,rep(1,Tn))+outer(b,k))
    am<-rowMeans(lm_); bbar<-rowMeans(bfun(mu)); W<-mu/mean(mu)
    Z<-(lm_-am)*sqrt(W); sv<-svd(Z); u1<-sv$u[,1]; v1<-sv$v[,1]
    if(sum(u1)<0){u1<--u1;v1<--v1}
    b2<-u1/sum(u1); k2<-sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
    k2<-k2-mean(k2); sb<-sum(b2); b2<-b2/sb; k2<-k2*sb
    a2<-am-bbar
    d<-max(abs(a2-a),abs(b2-b),abs(k2-k)/max(1,max(abs(k)))); a<-a2;b<-b2;k<-k2
    if(d<tol) break }
  list(a=a,b=b,k=k) }

UA1 <- u_anscombe(FALSE,TRUE); UA2 <- u_anscombe(TRUE,TRUE)
UA3 <- u_anscombe(TRUE,FALSE); UH <- u_huber(1.345); UH2 <- u_huber(3)
LAB <- c("標準 Lee-Carter","中心化＋加權","分離式","Firth","Poisson 擬分數",
         "Anscombe 近似平均、單位權","Anscombe 精確平均、單位權",
         "Anscombe 精確平均、精確變異","有界影響 Huber c=1.345",
         "有界影響 Huber c=3")
res <- list()
for (ni in seq_along(NS)) {
  N <- NS[ni]
  E <- outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu <- E*mt
  K <- length(LAB)
  aM <- array(NA_real_,c(REPS,A,K)); bS<-eS<-dr<-matrix(NA_real_,REPS,K)
  fl <- rep(0,K)
  for (r in seq_len(REPS)) {
    set.seed(SEED+17000+r)
    D <- matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    fits <- vector("list",K)
    fits[[1]] <- lc_svd_fit(log(pmax(D,c0)/E))
    fits[[2]] <- tryCatch(lc_alpha_only(D,E,bfun=BF,weight=TRUE),error=function(e)NULL)
    fits[[3]] <- tryCatch(lc_split(D,E,bfun=BF),error=function(e)NULL)
    fits[[4]] <- tryCatch(lc_poisson_firth(D,E,firth=TRUE),error=function(e)NULL)
    fits[[5]] <- tryCatch(lc_count_fit(D,E,u_poisson()),error=function(e)NULL)
    for (j in 6:10) {
      uu <- list(UA1,UA2,UA3,UH,UH2)[[j-5]]
      fits[[j]] <- tryCatch(lc_count_fit(D,E,uu),error=function(e)NULL)
    }
    for (j in seq_len(K)) {
      f <- fits[[j]]
      if (is.null(f)||!all(is.finite(f$a))||!all(is.finite(f$b))) next
      if (!is.null(f$floor_hit) && f$floor_hit) fl[j] <- fl[j]+1
      aM[r,,j]<-f$a-a0; bS[r,j]<-sum((f$b-b0)^2)/sum(b0^2)
      dr[r,j]<-(f$k[Tn]-f$k[1])/(Tn-1)-(k0[Tn]-k0[1])/(Tn-1)
      eS[r,j]<-tryCatch(life_table(exp(f$a+f$b*f$k[Tn]))$e0,error=function(e)NA)
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }
  res[[ni]] <- data.frame(N=N, 估計量=LAB,
    alpha中位=round(sapply(seq_len(K),function(j) median(abs(apply(aM[,,j],2,median,na.rm=TRUE)))),4),
    alpha最大=round(sapply(seq_len(K),function(j) max(abs(apply(aM[,,j],2,median,na.rm=TRUE)))),4),
    SSE_beta=round(apply(bS,2,median,na.rm=TRUE),4),
    漂移偏誤=round(apply(dr,2,median,na.rm=TRUE),4),
    e0偏誤=round(apply(eS,2,median,na.rm=TRUE)-e0t,3),
    e0標準差=round(apply(eS,2,sd,na.rm=TRUE),3),
    下界生效次數=fl, stringsAsFactors=FALSE)
  cat(sprintf("\n===== N = %.0e =====\n", N)); print(res[[ni]][,-1], row.names=FALSE)
}
tab <- do.call(rbind,res)
write.csv(tab,"output/tables/tableQ1_count_scale.csv",row.names=FALSE)
cat("\n已輸出 tableQ1_count_scale.csv\n")
