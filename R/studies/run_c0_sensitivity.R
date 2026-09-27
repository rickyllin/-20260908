###############################################################################
# 零格替代值 c0 的敏感度：校正能否消掉這個任意的選擇
#
#   問題。本文已指出偏誤的變號點 mu* 隨 c0 移動（0.1/0.3/0.5/1.0 ->
#   0.134/0.533/0.923/1.509），亦即「哪些年齡被高估」是分析者選出來的。
#   計量經濟學對 log(y+c) 的批評（Chen & Roth 2024）指出這類估計量所估的
#   東西隨 c 而變，不是良定義的參數。
#
#   可檢驗的預測。b(mu;c0) 的定義中已含 c0，故扣除它之後，
#   估計量對 c0 的選擇應\textbf{不再敏感}——校正把那個任意的選擇消掉了。
#   若成立，這是支持本文作法（保留 c0 並扣除其後果，而非調整 c0）的
#   直接證據；若不成立，則須說明殘餘的敏感度有多大。
#
# 輸出：output/tables/tableZ_c0_sensitivity.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
SEED <- 20260926; REPS <- 100; NS <- c(1e4, 5e4); YEARS <- 2001:2024
C0S <- c(0.1, 0.3, 0.5, 1.0)
dat <- load_data(sex="Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]; A<-nrow(D0); Tn<-ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage <- E0[,Tn]/sum(E0[,Tn]); e0_true <- life_table(mtrue[,Tn])$e0
BF <- lapply(C0S, make_logbias)          # 各 c0 一張查表
names(BF) <- as.character(C0S)

res <- list()
for (N in NS) {
  E <- outer(N*wage, rep(1,Tn)); dimnames(E) <- dimnames(D0); mu <- E*mtrue
  for (ci in seq_along(C0S)) {
    c0 <- C0S[ci]
    aU <- aC <- matrix(NA_real_, REPS, A); eU <- eC <- numeric(REPS)
    for (r in seq_len(REPS)) {
      set.seed(SEED + 1000*which(NS==N) + r)
      D <- matrix(rpois(length(E), mu), A, dimnames=dimnames(E))
      ## 未校正：僅更換 c0
      f0 <- lc_svd_fit(log(pmax(D,c0)/E))
      aU[r,] <- f0$a - truth$a
      m0 <- exp(outer(f0$a,rep(1,Tn))+outer(f0$b,f0$k))
      eU[r] <- tryCatch(life_table(m0[,Tn])$e0, error=function(e) NA_real_)
      ## 校正（僅校正 alpha ＋ 加權），b 依同一個 c0 計算
      f1 <- tryCatch(lc_alpha_only(D,E,c0=c0,bfun=BF[[ci]],weight=TRUE),
                     error=function(e) NULL)
      if (!is.null(f1)) {
        aC[r,] <- f1$a - truth$a
        m1 <- exp(outer(f1$a,rep(1,Tn))+outer(f1$b,f1$k))
        eC[r] <- tryCatch(life_table(m1[,Tn])$e0, error=function(e) NA_real_)
      }
    }
    res[[length(res)+1]] <- data.frame(N=N, c0=c0,
      未校正_a中位=round(median(abs(apply(aU,2,median,na.rm=TRUE))),4),
      未校正_a最大=round(max(abs(apply(aU,2,median,na.rm=TRUE))),4),
      未校正_e0偏誤=round(median(eU,na.rm=TRUE)-e0_true,3),
      校正後_a中位=round(median(abs(apply(aC,2,median,na.rm=TRUE))),4),
      校正後_a最大=round(max(abs(apply(aC,2,median,na.rm=TRUE))),4),
      校正後_e0偏誤=round(median(eC,na.rm=TRUE)-e0_true,3))
    cat(sprintf("N=%.0e c0=%.1f | 未校正 a中位 %.4f a最大 %.4f e0 %+.3f | 校正後 a中位 %.4f a最大 %.4f e0 %+.3f\n",
      N, c0, res[[length(res)]]$未校正_a中位, res[[length(res)]]$未校正_a最大,
      res[[length(res)]]$未校正_e0偏誤, res[[length(res)]]$校正後_a中位,
      res[[length(res)]]$校正後_a最大, res[[length(res)]]$校正後_e0偏誤))
  }
}
tab <- do.call(rbind,res)
write.csv(tab,"output/tables/tableZ_c0_sensitivity.csv",row.names=FALSE)
cat("\n=== 對 c0 的敏感度（同一 N 內，跨 c0 的全距）===\n")
for (N in NS) {
  d <- tab[tab$N==N,]
  cat(sprintf("N=%.0e  未校正：a最大 全距 %.4f、e0 全距 %.3f 歲\n", N,
    diff(range(d$未校正_a最大)), diff(range(d$未校正_e0偏誤))))
  cat(sprintf("        校正後：a最大 全距 %.4f、e0 全距 %.3f 歲   <-- 愈小表示校正消掉了 c0 的任意性\n",
    diff(range(d$校正後_a最大)), diff(range(d$校正後_e0偏誤))))
}
cat("\n已輸出 tableZ_c0_sensitivity.csv\n")
