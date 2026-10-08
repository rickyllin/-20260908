###############################################################################
# 邊界配適的系統性回報（B-6：Pawel et al. 2026, TAS 80(1): 31-48）
#
#   該文檢視 482 篇模擬研究，只有 23% 提到不收斂或演算法失敗、14% 說明如何
#   處理。本文的模擬在小人口下幼年組必然出現整列全零，此時 Poisson 最大概似
#   的 alpha_x 落在邊界（對數概似沿 alpha_x -> -無窮 單調而漸近），
#   V2 只零星提過一例（人數一萬、5-9 歲的 -13.03）。此處改為系統性回報：
#   一百次重複中有多少次、哪些年齡出現，以及如何處理。
#
# 輸出：output/tables/tableB6_boundary.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
SEED <- 20260926; REPS <- 100; NS <- c(2e3, 5e3, 1e4, 5e4, 2e5); c0 <- 0.5
dat <- load_data(sex="Female"); keep <- which(dat$years %in% 2001:2024)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]
A <- nrow(D0); Tn <- ncol(D0); ages <- rownames(D0)
tr <- lc_poisson_firth(D0,E0,firth=FALSE)
mt <- exp(outer(tr$a,rep(1,Tn))+outer(tr$b,tr$k)); wage <- E0[,Tn]/sum(E0[,Tn])
THRESH <- -12     # alpha_x 低於此值視為落在邊界（真值最低者約 -9.0）
cat(sprintf("真值 alpha 的最小值為 %.3f，故以 %g 為邊界判準\n\n", min(tr$a), THRESH))
rows <- list()
for (N in NS) {
  E <- outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu <- E*mt
  zero_row <- rep(0,A); bnd_mle <- matrix(0,REPS,A); nz <- 0
  for (r in seq_len(REPS)) {
    set.seed(SEED+11000+r)
    D <- matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    zr <- rowSums(D)==0; zero_row <- zero_row + zr; nz <- nz + sum(zr)
    f <- tryCatch(lc_poisson_firth(D,E,firth=FALSE), error=function(e) NULL)
    if (!is.null(f)) bnd_mle[r,] <- as.numeric(f$a < THRESH | !is.finite(f$a))
  }
  nrep_any <- sum(apply(bnd_mle,1,function(v) any(v>0)))
  rows[[length(rows)+1]] <- data.frame(
    N=N, 整列全零年齡數_平均=round(nz/REPS,2),
    有邊界配適的重複次數=nrep_any,
    邊界配適的年齡=paste(ages[colSums(bnd_mle)>0], collapse=", "),
    stringsAsFactors=FALSE)
  cat(sprintf("N=%7.0f｜整列全零 %.2f 個年齡／次｜100 次中 %3d 次出現邊界配適｜年齡：%s\n",
      N, nz/REPS, nrep_any,
      ifelse(any(colSums(bnd_mle)>0), paste(ages[colSums(bnd_mle)>0],collapse=", "), "—")))
}
tab <- do.call(rbind, rows)
write.csv(tab,"output/tables/tableB6_boundary.csv",row.names=FALSE)
cat("\n已輸出 tableB6_boundary.csv\n")
