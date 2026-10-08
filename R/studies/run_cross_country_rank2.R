###############################################################################
# 跨國的適用性評估（三）：受控的偏離秩一
#
#   對應審稿意見 2.10：「Add at least one scenario with a controlled departure
#   from rank one, for example a second age-time component of increasing
#   magnitude. Plot estimator performance against the magnitude of
#   misspecification.」並回答其所提的問題——偏離到何種程度時插入式校正
#   不再改善目標量。
#
#   設計。以各國全國資料的第二個奇異成分為偏離的方向（而非人造的方向，
#   如此偏離的形狀是該國真實資料所呈現者），真值取
#       log m = alpha_x + beta_x kappa_t + delta * u2_x v2_t * d2，
#   delta = 0 即秩一。delta 掃 0、0.25、0.5、1.0、2.0。
#   delta = 1 時第二成分的量即該國全國資料所實際呈現的量。
#
#   取臺灣與美國兩國，因其秩一佔比與 kappa 線性度分居兩端
#   （臺灣 0.896/0.982、美國 0.623/0.711）。
#
# 輸出：output/tables/tableX3_rank2.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/hmd_load.R")
SEED <- 20261008; REPS <- 100; N <- 5e4; c0 <- 0.5; BF <- make_logbias(c0)
DELTAS <- c(0, 0.25, 0.5, 1, 2)

lc_split <- function(D, E, bfun, maxit = 20, tol = 1e-9) {
  Tn <- ncol(D); lm_ <- log(pmax(D, c0)/E)
  f <- lc_svd_fit(lm_); a<-f$a;b<-f$b;k<-f$k
  for (it in seq_len(maxit)) {
    mu <- E*exp(outer(a,rep(1,Tn))+outer(b,k))
    am<-rowMeans(lm_); bb<-rowMeans(bfun(mu)); W<-mu/mean(mu)
    Z<-(lm_-am)*sqrt(W); sv<-svd(Z); u1<-sv$u[,1]; v1<-sv$v[,1]
    if(sum(u1)<0){u1<--u1;v1<--v1}
    b2<-u1/sum(u1); k2<-sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
    k2<-k2-mean(k2); sb<-sum(b2); b2<-b2/sb; k2<-k2*sb; a2<-am-bb
    d<-max(abs(a2-a),abs(b2-b),abs(k2-k)/max(1,max(abs(k)))); a<-a2;b<-b2;k<-k2
    if(d<tol) break }
  list(a=a,b=b,k=k) }

rows <- list()
for (nm in c("Taiwan","US")) {
  x <- get_country(nm); A<-nrow(x$D); Tn<-ncol(x$D)
  tr <- lc_poisson_firth(round(x$D), x$E, firth=FALSE)
  a0<-tr$a; b0<-tr$b; k0<-tr$k
  L1 <- outer(a0,rep(1,Tn)) + outer(b0,k0)
  ## 由實際資料取第二成分
  Lobs <- log(pmax(x$D,c0)/x$E); R <- Lobs - L1
  sv <- svd(R); u2<-sv$u[,1]; v2<-sv$v[,1]; d2<-sv$d[1]
  w <- x$E[,Tn]/sum(x$E[,Tn]); E <- outer(N*w,rep(1,Tn)); dimnames(E)<-dimnames(x$D)
  for (de in DELTAS) {
    Ltrue <- L1 + de * d2 * outer(u2, v2)
    mtrue <- exp(Ltrue); mu <- E*mtrue
    e0t <- life_table(mtrue[,Tn])$e0
    ## 真值的 alpha 定義為該尺度的列平均（偏離秩一時 alpha 不再唯一，
    ## 故以列平均為目標，與閉式所刻畫的對象一致）
    aT <- rowMeans(Ltrue)
    aM <- array(NA_real_,c(REPS,A,2)); eS <- matrix(NA_real_,REPS,2)
    for (r in seq_len(REPS)) {
      set.seed(SEED+13000+r)
      D <- matrix(rpois(length(E), mu), A, dimnames=dimnames(E))
      fits <- list(lc_svd_fit(log(pmax(D,c0)/E)),
                   tryCatch(lc_split(D,E,BF), error=function(e) NULL))
      for (j in 1:2) { f<-fits[[j]]
        if (is.null(f)||!all(is.finite(f$a))) next
        aM[r,,j]<-f$a-aT
        eS[r,j]<-tryCatch(life_table(exp(f$a+f$b*f$k[Tn]))$e0,error=function(e)NA) }
    }
    o1<-apply(aM[,,1],2,median,na.rm=TRUE); o2<-apply(aM[,,2],2,median,na.rm=TRUE)
    rows[[length(rows)+1]] <- data.frame(國家=nm, delta=de,
      第二成分佔比 = round((de*d2)^2/((de*d2)^2+sum((L1-rowMeans(L1))^2)),4),
      標準_a最大=round(max(abs(o1)),4), 分離式_a最大=round(max(abs(o2)),4),
      改善倍數=round(max(abs(o1))/max(abs(o2)),2),
      標準_e0=round(median(eS[,1],na.rm=TRUE)-e0t,3),
      分離式_e0=round(median(eS[,2],na.rm=TRUE)-e0t,3),
      stringsAsFactors=FALSE)
  }
  cat(sprintf("  %s 完成\n", nm))
}
tab <- do.call(rbind, rows)
write.csv(tab,"output/tables/tableX3_rank2.csv",row.names=FALSE)
cat("\n### 偏離秩一時校正是否仍改善（人數五萬）\n"); print(tab, row.names=FALSE)
