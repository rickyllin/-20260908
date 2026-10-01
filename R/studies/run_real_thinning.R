###############################################################################
# 真實資料分析：以卜瓦松抽薄產生小人口
#
#   與既有模擬設計的差別。既有設計取 D ~ Poisson(E * m_fitted)，
#   亦即\textbf{假設 LC 模型完全正確}；小人口的資料因而不含任何不合模型的
#   結構。本檔改以\textbf{實際觀測到的}全國死亡數抽薄：
#       D_small_{x,t} ~ Poisson( p * D_obs_{x,t} ),   E_small = p * E_obs
#   其中 p = N / 全國人數。如此小人口繼承台灣死亡率真實的逐年不規則性
#   （戰後嬰兒潮、事故、疫情等所有不合秩一雙線性結構的部分），
#   而非由模型生成。真值仍取全國資料的 LC 配適。
#
#   這是對推導的一項更嚴格的檢驗：式 b(mu;c0) 的推導依賴卜瓦松假設，
#   而抽薄後 E[D_small] = p * D_obs，其中 D_obs 本身偏離 LC 配適值。
#   若校正在此仍有效，表示它對「模型誤設」有一定的穩健性；
#   若失效，則表示其效力依賴 LC 模型正確——兩種結果都是可寫的。
#
#   註：抽薄方式依 0929 會議意見由二項改為卜瓦松。兩者的期望值相同
#   （皆為 p * D_obs），但二項抽薄有 D_small <= D_obs 的上界，變異數為
#   p(1-p) D_obs 而非 p D_obs，在 p 很小時兩者差異雖小，卜瓦松抽薄卻
#   與本文推導所假設的分配完全一致，不必再討論上界的影響。
#   被打破的不是卜瓦松假設，而是「mu 落在秩一雙線性結構上」這個假設。
#
# 輸出：output/tables/tableW_real_thinning.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")

SEED <- 20260926; REPS <- 100; NS <- c(1e4, 5e4, 2e5); YEARS <- 2001:2024
dat <- load_data(sex="Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]; A<-nrow(D0); Tn<-ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
dtrue <- (truth$k[Tn]-truth$k[1])/(Tn-1)
Npop  <- sum(E0[,Tn])
BF <- make_logbias(0.5)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
e0_true <- life_table(mtrue[,Tn])$e0
cat(sprintf("全國末年人數 = %.0f；真值 drift = %.4f\n", Npop, dtrue))
cat(sprintf("觀測值對配適值的偏離：log(D_obs/E) - (a+b k) 的標準差 = %.4f\n",
  sd(log(pmax(D0,0.5)/E0) - (outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k)))))

SPEC <- list(
  list(id="A 標準 LC",        f=function(D,E) lc_svd_fit(log(pmax(D,0.5)/E))),
  list(id="B 卜瓦松 MLE",     f=function(D,E) lc_poisson_firth(D,E,firth=FALSE)),
  list(id="C Firth",          f=function(D,E) lc_poisson_firth(D,E,firth=TRUE)),
  list(id="D 中心化",         f=function(D,E) lc_analytic(D,E,bfun=BF,shrink=1)),
  list(id="E 中心化 + 加權",  f=function(D,E) {
      g <- lc_analytic(D,E,bfun=BF,shrink=1)
      mu <- E*exp(outer(g$a,rep(1,ncol(D)))+outer(g$b,g$k))
      lmc <- log(pmax(D,0.5)/E) - BF(mu); W <- mu/mean(mu)
      am <- rowSums(W*lmc)/rowSums(W); Z <- (lmc-am)*sqrt(W)
      sv <- svd(Z); u1<-sv$u[,1]; v1<-sv$v[,1]
      if (sum(u1)<0){u1<--u1; v1<--v1}
      bb<-u1/sum(u1); kk<-sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
      kk<-kk-mean(kk); sb<-sum(bb); list(a=am,b=bb/sb,k=kk*sb)})
)
ESTS <- vapply(SPEC,`[[`,character(1),"id"); J<-length(SPEC)
ok <- function(f) all(is.finite(f$a))&&all(is.finite(f$b))&&all(is.finite(f$k))&&max(abs(f$k))<1e3

res <- list()
for (N in NS) {
  p <- N/Npop
  E <- p*E0; dimnames(E) <- dimnames(D0)
  aM <- array(NA_real_,c(REPS,A,J)); bS<-dS<-cS<-eS<-matrix(NA_real_,REPS,J)
  div<-integer(J); zc<-numeric(REPS)
  for (r in seq_len(REPS)) {
    set.seed(SEED + 7000*which(NS==N) + r)
    D <- matrix(rpois(length(D0), p * D0), A, dimnames=dimnames(D0))
    zc[r] <- sum(D==0)
    for (j in seq_len(J)) {
      f <- tryCatch(SPEC[[j]]$f(D,E), error=function(e) NULL)
      if (is.null(f)||!ok(f)) { div[j]<-div[j]+1L; next }
      aM[r,,j] <- f$a-truth$a
      bS[r,j] <- sum((f$b-truth$b)^2)/sum(truth$b^2)
      cS[r,j] <- if (sd(f$b)<1e-12) NA_real_ else cor(f$b,truth$b)
      dS[r,j] <- (f$k[Tn]-f$k[1])/(Tn-1)-dtrue
      mh <- exp(outer(f$a,rep(1,Tn))+outer(f$b,f$k))
      eS[r,j] <- tryCatch(life_table(mh[,Tn])$e0, error=function(e) NA_real_)
    }
    if (r%%25==0) cat(sprintf("  N=%.0e rep %d/%d\n",N,r,REPS))
  }
  res[[as.character(N)]] <- data.frame(N=N, 估計量=ESTS, 平均零格=round(mean(zc),1),
    alpha中位=round(sapply(seq_len(J),function(j) median(abs(apply(aM[,,j],2,median,na.rm=TRUE)))),4),
    alpha最大=round(sapply(seq_len(J),function(j) max(abs(apply(aM[,,j],2,median,na.rm=TRUE)))),4),
    SSE_beta=round(apply(bS,2,median,na.rm=TRUE),4),
    cor_beta=round(apply(cS,2,median,na.rm=TRUE),3),
    漂移偏誤=round(apply(dS,2,median,na.rm=TRUE),4),
    e0偏誤=round(apply(eS,2,median,na.rm=TRUE)-e0_true,3),
    e0四分位距=round(apply(eS,2,function(v) IQR(v,na.rm=TRUE)),3),
    發散率=round(div/REPS,3), stringsAsFactors=FALSE)
  cat(sprintf("\n===== N = %.0e（p=%.5f，平均零格 %.1f/%d）=====\n", N, p, mean(zc), A*Tn))
  print(res[[as.character(N)]][,-1], row.names=FALSE)
}
tab <- do.call(rbind,res)
dir.create("output/tables",recursive=TRUE,showWarnings=FALSE)
write.csv(tab,"output/tables/tableW_real_thinning.csv",row.names=FALSE)
cat("\n已輸出 tableW_real_thinning.csv\n")
