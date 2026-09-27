###############################################################################
# 推導的變異數最適權重 vs 冪次權重
#
#   受檢驗的推導（見 R/core/mquantile_lc.R 檔頭）：
#     w_x ∝ mu_x ( 2 Phi(c sqrt(mu_x)) - 1 ) / beta_{c sqrt(mu_x)}
#   其兩個極限為 sqrt(mu)（檢查函數）與 mu（最小平方），內插變數為 c sqrt(mu)。
#
#   預測：推導權重應優於 sqrt(mu) 與 mu 兩個固定冪次，因為同一筆資料的
#   各年齡 c sqrt(mu) 橫跨兩個數量級，沒有單一冪次能同時服務兩端。
#   若推導權重僅與 mu 相當，則表示本問題中 c sqrt(mu) 大多落在「最小平方」
#   一側，推導雖正確但無實務差異——此亦為可寫的結果。
#
#   種子與 run_eda_to_correction.R 完全相同。
# 輸出：output/tables/tableU_derived_weight.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/mquantile_lc.R")

SEED <- 20260926; REPS <- 100; NS <- c(1e4, 5e4, 2e5); YEARS <- 2001:2024
dat <- load_data(sex="Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]; A<-nrow(D0); Tn<-ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage <- E0[,Tn]/sum(E0[,Tn]); drift_true <- (truth$k[Tn]-truth$k[1])/(Tn-1)
cat(sprintf("真值 drift = %.4f\n", drift_true))

SPEC <- list()
for (cc in c(0.7, 1.345)) for (wm in c("none","sqrt","mu","opt"))
  SPEC[[length(SPEC)+1]] <- list(
    id = sprintf("c=%-5s w=%s", cc, wm), c = cc, wm = wm)
ESTS <- vapply(SPEC, `[[`, character(1), "id"); J <- length(SPEC)

fit1 <- function(D,E,sp) switch(sp$wm,
  none = lc_mquantile(D,E,0.5,k_c=sp$c,wmode="none"),
  sqrt = lc_mquantile(D,E,0.5,k_c=sp$c,wmode="pow",wpow=0.5),
  mu   = lc_mquantile(D,E,0.5,k_c=sp$c,wmode="mu"),
  opt  = lc_mquantile(D,E,0.5,k_c=sp$c,wmode="opt"))
ok <- function(f) all(is.finite(f$a))&&all(is.finite(f$b))&&all(is.finite(f$k))&&max(abs(f$k))<1e3

res <- list()
for (N in NS) {
  E <- outer(N*wage, rep(1,Tn)); dimnames(E) <- dimnames(D0); mu_exp <- E*mtrue
  aM <- array(NA_real_, c(REPS,A,J)); bS<-dS<-cS<-matrix(NA_real_,REPS,J); div<-integer(J)
  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000*which(NS==N) + r)
    D <- matrix(rpois(length(E), mu_exp), A, dimnames=dimnames(E))
    for (j in seq_len(J)) {
      f <- tryCatch(fit1(D,E,SPEC[[j]]), error=function(e) NULL)
      if (is.null(f)||!ok(f)) { div[j]<-div[j]+1L; next }
      aM[r,,j] <- f$a-truth$a
      bS[r,j] <- sum((f$b-truth$b)^2)/sum(truth$b^2)
      cS[r,j] <- if (sd(f$b)<1e-12) NA_real_ else cor(f$b,truth$b)
      dS[r,j] <- (f$k[Tn]-f$k[1])/(Tn-1)-drift_true
    }
    if (r%%25==0) cat(sprintf("  N=%.0e rep %d/%d\n",N,r,REPS))
  }
  res[[as.character(N)]] <- data.frame(N=N, 估計量=ESTS,
    alpha中位=round(sapply(seq_len(J),function(j) median(abs(apply(aM[,,j],2,median,na.rm=TRUE)))),4),
    alpha最大=round(sapply(seq_len(J),function(j) max(abs(apply(aM[,,j],2,median,na.rm=TRUE)))),4),
    SSE_beta=round(apply(bS,2,median,na.rm=TRUE),4),
    cor_beta=round(apply(cS,2,median,na.rm=TRUE),3),
    漂移偏誤=round(apply(dS,2,median,na.rm=TRUE),4),
    發散率=round(div/REPS,3), stringsAsFactors=FALSE)
  cat(sprintf("\n===== N = %.0e =====\n",N)); print(res[[as.character(N)]][,-1], row.names=FALSE)
}
tab <- do.call(rbind,res)
dir.create("output/tables",recursive=TRUE,showWarnings=FALSE)
write.csv(tab,"output/tables/tableU_derived_weight.csv",row.names=FALSE)
cat("\n已輸出 tableU_derived_weight.csv\n")
