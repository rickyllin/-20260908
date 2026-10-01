###############################################################################
# 人口年齡結構是否影響偏誤（回應 0929 會議：Whether population structure
# influences the bias）
#
#   本文的閉式 b(mu; c0) 只透過期望死亡數 mu_{x,t} = E_{x,t} m_{x,t} 進入偏誤，
#   而 E_{x,t} = N * w_x 中的 w_x 正是年齡結構。閉式因此對本問題給出一個
#   可檢驗的預測：在總人數 N 固定下改變 w_x，標準 LC 的逐年齡偏誤應隨之
#   改變，且改變的方向與幅度完全由 b(mu) 決定；而扣除該量之後，
#   剩餘偏誤應與結構無關。這比「偏誤是否隨 N 變化」更嚴格，
#   因為此處 N 不變，變的只有資訊在年齡間的分配。
#
#   三組結構皆取自臺灣女性的實際曝露數：
#     年輕結構 = 2001 年（65 歲以上占 8.5%）
#     基準結構 = 2024 年（20.3%）
#     高齡結構 = 以 2024 年結構作指數傾斜，使 65 歲以上占約三成,
#                對應人口外流的鄉鎮市區。
#
# 輸出：output/tables/tableM2_popstructure.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")

SEED <- 20260926; REPS <- 100; NS <- c(1e4, 5e4); YEARS <- 2001:2024
dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0); ages <- rownames(D0)
truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
e0_true <- life_table(mtrue[, Tn])$e0
BF <- make_logbias(0.5)

w_young <- E0[, 1]  / sum(E0[, 1])
w_base  <- E0[, Tn] / sum(E0[, Tn])
tilt    <- exp(0.085 * (seq_len(A) - 1))
w_old   <- w_base * tilt; w_old <- w_old / sum(w_old)
i65 <- 15:A
STR <- list("年輕（2001）" = w_young, "基準（2024）" = w_base, "高齡（傾斜）" = w_old)
cat(sprintf("65+ 占比：%s\n",
    paste(sprintf("%s=%.3f", names(STR), sapply(STR, function(w) sum(w[i65]))),
          collapse = "  ")))

FIT <- list("標準 LC"      = function(D, E) lc_svd_fit(log(pmax(D, 0.5) / E)),
            "Firth"        = function(D, E) lc_poisson_firth(D, E, firth = TRUE),
            "中心化＋加權" = function(D, E) lc_alpha_only(D, E, bfun = BF, weight = TRUE))

res <- list(); prof <- list()
for (ni in seq_along(NS)) for (si in seq_along(STR)) {
  N <- NS[ni]; w <- STR[[si]]
  E  <- outer(N * w, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu <- E * mtrue
  bpred <- rowMeans(BF(mu))                      # 閉式預測的逐年齡 alpha 偏誤
  K <- length(FIT)
  aM <- array(NA_real_, c(REPS, A, K)); bS <- eS <- matrix(NA_real_, REPS, K)
  for (r in seq_len(REPS)) {
    set.seed(SEED + 23000 * ni + 311 * si + r)
    D <- matrix(rpois(length(E), mu), A, dimnames = dimnames(E))
    for (j in seq_len(K)) {
      f <- tryCatch(FIT[[j]](D, E), error = function(e) NULL)
      if (is.null(f) || !all(is.finite(f$a)) || !all(is.finite(f$b))) next
      aM[r, , j] <- f$a - truth$a
      bS[r, j] <- sum((f$b - truth$b)^2) / sum(truth$b^2)
      mh <- exp(outer(f$a, rep(1, Tn)) + outer(f$b, f$k))
      eS[r, j] <- tryCatch(life_table(mh[, Tn])$e0, error = function(e) NA_real_)
    }
  }
  aMed <- sapply(seq_len(K), function(j) apply(aM[, , j], 2, median, na.rm = TRUE))
  cr <- cor(bpred, aMed[, 1]); sl <- unname(coef(lm(aMed[, 1] ~ bpred))[2])
  res[[length(res) + 1]] <- data.frame(
    N = N, 結構 = names(STR)[si], 老年占比 = round(sum(w[i65]), 3),
    估計量 = names(FIT),
    alpha中位 = round(apply(abs(aMed), 2, median), 4),
    alpha最大 = round(apply(abs(aMed), 2, max), 4),
    SSE_beta  = round(apply(bS, 2, median, na.rm = TRUE), 4),
    e0偏誤    = round(apply(eS, 2, median, na.rm = TRUE) - e0_true, 3),
    e0標準差  = round(apply(eS, 2, sd, na.rm = TRUE), 3),
    閉式相關  = round(cr, 3), 閉式斜率 = round(sl, 3),
    stringsAsFactors = FALSE)
  prof[[length(prof) + 1]] <- data.frame(
    N = N, 結構 = names(STR)[si], 年齡 = ages, 曝露占比 = round(w, 5),
    閉式預測 = round(bpred, 4), 標準LC實測 = round(aMed[, 1], 4),
    中心化實測 = round(aMed[, 3], 4), stringsAsFactors = FALSE)
  cat(sprintf("N=%.0e %s：閉式相關 %.3f、斜率 %.3f\n", N, names(STR)[si], cr, sl))
}
tab <- do.call(rbind, res)
write.csv(tab, "output/tables/tableM2_popstructure.csv", row.names = FALSE)
write.csv(do.call(rbind, prof), "output/tables/tableM2_popstructure_byage.csv",
          row.names = FALSE)
print(tab, row.names = FALSE)
cat("\n已輸出 tableM2_popstructure.csv\n")
