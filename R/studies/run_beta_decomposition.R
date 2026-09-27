###############################################################################
# beta 的偏誤分解：確定性成分（可由中心化移除）與隨機成分（須靠加權）
#
#   推導。beta_hat 來自中心化矩陣 Z 的 SVD，其中 Z_{x,t} = ell_{x,t} - mean_t ell_{x,t}。
#   取期望（並用 sum_t kappa_t = 0）：
#       E[Z_{x,t}] = beta_x kappa_t + ( b(mu_{x,t}) - bbar_x ),   bbar_x = mean_t b(mu_{x,t})
#   故期望矩陣\textbf{不是}秩一訊號，而多了一個確定性的擾動
#       Delta_{x,t} = b(mu_{x,t}) - bbar_x .
#   若 b(mu_{x,t}) 在列內不隨 t 變動則 Delta = 0，beta 不受對數偏誤影響；
#   但死亡率隨年下降使 mu_{x,t} 隨 t 變動，而 b 為非線性，故 Delta != 0。
#
#   本檔把 beta 的總誤差分解為兩部分，且分解是精確的（不需一階近似）：
#     (甲) 確定性成分：對\textbf{無噪音}的期望矩陣 beta kappa' + Delta 作 SVD，
#          所得 beta_hat 與真值之差。此即 Delta 單獨造成的誤差。
#     (乙) 隨機成分：總誤差減去 (甲)。來源是異質變異噪音旋轉主奇異向量
#          （Zhang-Cai-Wu 的 HeteroPCA 機制），與 Delta 無關。
#   中心化移除的是 (甲)；加權處理的是 (乙)。
#
#   可檢驗的預測：若 (乙) 遠大於 (甲)，則中心化單獨不會改善 beta——
#   這正是實測所見（中心化使 SSE(beta) 變差）。本檔給出兩者的相對大小。
#
# 輸出：output/tables/tableV_beta_decomp.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")

NS <- c(1e4, 5e4, 2e5, 1e6); YEARS <- 2001:2024
dat <- load_data(sex="Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]; A<-nrow(D0); Tn<-ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage <- E0[,Tn]/sum(E0[,Tn]); dtrue <- (truth$k[Tn]-truth$k[1])/(Tn-1)
BF <- make_logbias(0.5)
SEED <- 20260926; REPS <- 100

out <- list()
for (N in NS) {
  E <- outer(N*wage, rep(1,Tn)); dimnames(E) <- dimnames(D0)
  mu <- E*mtrue
  b  <- BF(mu)
  bbar <- rowMeans(b)
  Delta <- b - bbar                                   # 確定性擾動
  ## --- 訊號與擾動的相對大小 ---
  Sig <- outer(truth$b, truth$k)                      # beta kappa'
  s1  <- svd(Sig)$d[1]
  nD  <- sqrt(sum(Delta^2)); nS <- sqrt(sum(Sig^2))
  ## --- (甲) 確定性成分：對無噪音期望矩陣作 SVD ---
  Zdet <- Sig + Delta
  sv <- svd(Zdet); u1 <- sv$u[,1]; v1 <- sv$v[,1]
  if (sum(u1) < 0) { u1 <- -u1; v1 <- -v1 }
  bdet <- u1/sum(u1); kdet <- sv$d[1]*v1*sum(u1); kdet <- kdet - mean(kdet)
  sse_det <- sum((bdet-truth$b)^2)/sum(truth$b^2)
  dr_det  <- (kdet[Tn]-kdet[1])/(Tn-1) - dtrue
  a_det   <- max(abs(bbar))                           # alpha 的確定性偏誤 = |bbar|
  ## --- 總誤差（模擬） ---
  sse <- dr <- amax <- numeric(REPS)
  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000*which(NS==N) + r)
    D <- matrix(rpois(length(E), mu), A, dimnames=dimnames(E))
    f <- lc_svd_fit(log(pmax(D,0.5)/E))
    sse[r] <- sum((f$b-truth$b)^2)/sum(truth$b^2)
    dr[r]  <- (f$k[Tn]-f$k[1])/(Tn-1) - dtrue
    amax[r]<- max(abs(f$a - truth$a))
  }
  out[[as.character(N)]] <- data.frame(
    N=N,
    `訊號 s1`       = round(s1,4),
    `||Delta||_F`   = round(nD,4),
    `||Delta||/||Sig||` = round(nD/nS,4),
    `alpha 確定性`  = round(a_det,4),
    `alpha 實測`    = round(median(amax),4),
    `SSE(b) 確定性(甲)` = round(sse_det,5),
    `SSE(b) 總(模擬)`   = round(median(sse),5),
    `甲占總比`      = round(sse_det/median(sse),4),
    `漂移 確定性(甲)`   = round(dr_det,5),
    `漂移 總(模擬)`     = round(median(dr),5),
    check.names=FALSE)
  cat(sprintf("N=%.0e  ||Delta||/||Sig||=%.4f  SSE(b): 確定性 %.5f / 總 %.5f (占 %.2f%%)  漂移: 確定性 %+.5f / 總 %+.5f\n",
      N, nD/nS, sse_det, median(sse), 100*sse_det/median(sse), dr_det, median(dr)))
}
tab <- do.call(rbind,out)
dir.create("output/tables",recursive=TRUE,showWarnings=FALSE)
write.csv(tab,"output/tables/tableV_beta_decomp.csv",row.names=FALSE)
cat("\n=== 完整表 ===\n"); print(t(tab))
cat("\n已輸出 tableV_beta_decomp.csv\n")
