###############################################################################
# 跨國的適用性評估（二）：閉式與校正是否逐國成立
#
#   對應審稿意見 2.10（內部驗證與外部效度須分開標示）與其 minor 的
#   「多個年齡結構」一項。此處的設計刻意只做內部驗證的跨國版本：
#   以各國全國資料的 Poisson 配適為真值，按該國自己的年齡結構縮放，
#   故資料恰好落在秩一結構上。其目的是檢驗閉式是否逐國成立，
#   不是主張跨國的實務優越性。
#
#   秩一適足性另單獨回報，因為它逐國不同而閉式不依賴它——
#   這一區分正是審稿意見 2.10 所要求的。
#
# 輸出：output/tables/tableX2_cross_country_validate.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/hmd_load.R")
SEED <- 20261008; REPS <- 100; NS <- c(1e4, 5e4); c0 <- 0.5
BF <- make_logbias(c0)

## alpha 等權、beta 與 kappa 加權的分離式（1012 V1 第柒節之三）
lc_split <- function(D, E, bfun, maxit = 20, tol = 1e-9) {
  Tn <- ncol(D); lm_ <- log(pmax(D, c0) / E)
  f <- lc_svd_fit(lm_); a <- f$a; b <- f$b; k <- f$k
  for (it in seq_len(maxit)) {
    mu <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    am <- rowMeans(lm_); bb <- rowMeans(bfun(mu)); W <- mu / mean(mu)
    Z <- (lm_ - am) * sqrt(W); sv <- svd(Z); u1 <- sv$u[,1]; v1 <- sv$v[,1]
    if (sum(u1) < 0) { u1 <- -u1; v1 <- -v1 }
    b2 <- u1 / sum(u1); k2 <- sv$d[1] * v1 * sum(u1) / sqrt(pmax(colMeans(W), 1e-12))
    k2 <- k2 - mean(k2); sb <- sum(b2); b2 <- b2 / sb; k2 <- k2 * sb
    a2 <- am - bb
    d <- max(abs(a2-a), abs(b2-b), abs(k2-k)/max(1,max(abs(k)))); a<-a2;b<-b2;k<-k2
    if (d < tol) break }
  list(a = a, b = b, k = k) }

rows <- list()
for (nm in ALL_COUNTRIES) {
  x <- get_country(nm); A <- nrow(x$D); Tn <- ncol(x$D)
  tr <- lc_poisson_firth(round(x$D), x$E, firth = FALSE)
  a0 <- tr$a; b0 <- tr$b; k0 <- tr$k
  mt <- exp(outer(a0, rep(1,Tn)) + outer(b0, k0))
  w  <- x$E[, Tn] / sum(x$E[, Tn])
  e0t <- life_table(mt[, Tn])$e0
  ## 秩一適足性（全國資料，抽樣誤差可忽略）
  Lc <- log(x$D / x$E); Lc[!is.finite(Lc)] <- NA
  Lc <- sweep(Lc, 1, rowMeans(Lc, na.rm = TRUE)); Lc[is.na(Lc)] <- 0
  sv1 <- svd(Lc)$d; r1 <- sv1[1]^2 / sum(sv1^2)
  ## kappa 的線性度（漂移假設的檢查）
  tt <- seq_len(Tn); lin <- summary(lm(k0 ~ tt))$r.squared

  for (N in NS) {
    E <- outer(N*w, rep(1,Tn)); dimnames(E) <- dimnames(x$D); mu <- E*mt
    pred <- rowMeans(BF(mu))                       # 閉式預測的 alpha 偏誤
    aM <- array(NA_real_, c(REPS, A, 3)); eS <- bS <- matrix(NA_real_, REPS, 3)
    for (r in seq_len(REPS)) {
      set.seed(SEED + 13000 + r)
      D <- matrix(rpois(length(E), mu), A, dimnames = dimnames(E))
      fits <- list(lc_svd_fit(log(pmax(D,c0)/E)),
                   tryCatch(lc_split(D,E,BF), error=function(e) NULL),
                   tryCatch(lc_poisson_firth(D,E,firth=TRUE), error=function(e) NULL))
      for (j in 1:3) { f <- fits[[j]]
        if (is.null(f) || !all(is.finite(f$a))) next
        aM[r,,j] <- f$a - a0; bS[r,j] <- sum((f$b-b0)^2)/sum(b0^2)
        eS[r,j] <- tryCatch(life_table(exp(f$a+f$b*f$k[Tn]))$e0, error=function(e) NA) }
    }
    obs <- apply(aM[,,1], 2, median, na.rm = TRUE)
    rows[[length(rows)+1]] <- data.frame(國家=nm, N=N,
      零格佔比 = round(mean(mu < c0), 3),
      閉式相關 = round(cor(obs, pred), 4),
      閉式斜率 = round(coef(lm(obs ~ pred))[2], 3),
      標準_a最大 = round(max(abs(obs)), 4),
      分離式_a最大 = round(max(abs(apply(aM[,,2],2,median,na.rm=TRUE))), 4),
      Firth_a最大 = round(max(abs(apply(aM[,,3],2,median,na.rm=TRUE))), 4),
      標準_e0 = round(median(eS[,1],na.rm=TRUE)-e0t, 3),
      分離式_e0 = round(median(eS[,2],na.rm=TRUE)-e0t, 3),
      Firth_e0 = round(median(eS[,3],na.rm=TRUE)-e0t, 3),
      秩一佔比 = round(r1,4), kappa線性R2 = round(lin,4),
      stringsAsFactors = FALSE)
  }
  cat(sprintf("  %s 完成\n", nm))
}
tab <- do.call(rbind, rows)
write.csv(tab, "output/tables/tableX2_cross_country_validate.csv", row.names=FALSE)
cat("\n### 閉式的跨國驗證與校正的成績\n")
print(tab[, c("國家","N","零格佔比","閉式相關","閉式斜率",
              "標準_a最大","分離式_a最大","Firth_a最大")], row.names=FALSE)
cat("\n### 平均餘命與結構診斷\n")
print(tab[, c("國家","N","標準_e0","分離式_e0","Firth_e0","秩一佔比","kappa線性R2")],
      row.names=FALSE)
