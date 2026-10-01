###############################################################################
# 與參考母體修勻法的比較（回應 0929 會議：「目前方法與 Partial SMR 孰優」）
#
#   設計。Partial SMR（Lee 2003）與 Whittaker ratio 都需要一個人數較多的
#   參考母體，其表現取決於參考母體與目標小區域的死亡率型態是否相近。
#   因此比較不能只在「兩者型態相同」的一個情境下進行，否則等於讓對照組
#   先看到答案。本文沿用 Yue, Wang and Wang (2019, Figure 2) 的七種比值情境
#   s_x = m^small_x / m^ref_x，涵蓋型態相同（0.8、1.0、1.2）與型態不同
#   （遞增、遞減、V 型、倒 V 型）兩類。
#
#   小區域的死亡數在各情境下完全相同（亂數種子只隨 N 與重複次數變動），
#   故三個不使用參考母體的估計量（標準 LC、Firth、本文中心化加權）
#   在各情境下的結果一致，可作為橫向比較的基準線。
#
# 輸出：output/tables/tableM1_smr_compare.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/partial_smr.R")

SEED <- 20260926; REPS <- 100; NS <- c(1e4, 5e4, 2e5); NREF <- 2e6
YEARS <- 2001:2024

dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0)
truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn])
e0_true <- life_table(mtrue[, Tn])$e0
BF <- make_logbias(0.5)
SCN <- ratio_scenarios(A)

## 不使用參考母體的三個估計量 --------------------------------------------------
FIT0 <- list(
  "標準 LC"        = function(D, E) lc_svd_fit(log(pmax(D, 0.5) / E)),
  "Firth"          = function(D, E) lc_poisson_firth(D, E, firth = TRUE),
  "中心化＋加權"   = function(D, E) lc_alpha_only(D, E, bfun = BF, weight = TRUE))

metric <- function(f) {
  if (is.null(f) || !all(is.finite(f$a)) || !all(is.finite(f$b)))
    return(c(NA_real_, NA_real_, NA_real_))
  mh <- exp(outer(f$a, rep(1, Tn)) + outer(f$b, f$k))
  c(sum((f$b - truth$b)^2) / sum(truth$b^2),
    tryCatch(life_table(mh[, Tn])$e0, error = function(e) NA_real_),
    NA_real_)
}

res <- list()
for (ni in seq_along(NS)) {
  N <- NS[ni]
  E   <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  Eref <- outer(NREF * wage, rep(1, Tn)); dimnames(Eref) <- dimnames(D0)
  mu  <- E * mtrue

  nm0 <- names(FIT0); nmS <- c("Partial SMR", "Whittaker ratio", "零格以 SMR 填補")
  lab <- c(nm0, as.vector(outer(nmS, names(SCN), paste, sep = " ｜ ")))
  K   <- length(lab)
  aM  <- array(NA_real_, c(REPS, A, K))
  bS  <- eS <- matrix(NA_real_, REPS, K)

  for (r in seq_len(REPS)) {
    set.seed(SEED + 11000 * ni + r)
    D <- matrix(rpois(length(E), mu), A, dimnames = dimnames(E))

    j <- 0L
    for (nm in nm0) {
      j <- j + 1L
      f <- tryCatch(FIT0[[nm]](D, E), error = function(e) NULL)
      if (!is.null(f) && all(is.finite(f$a)) && all(is.finite(f$b))) {
        aM[r, , j] <- f$a - truth$a
        bS[r, j]   <- sum((f$b - truth$b)^2) / sum(truth$b^2)
        mh <- exp(outer(f$a, rep(1, Tn)) + outer(f$b, f$k))
        eS[r, j] <- tryCatch(life_table(mh[, Tn])$e0, error = function(e) NA_real_)
      }
    }
    for (si in seq_along(SCN)) {
      s  <- SCN[[si]]
      mR <- mtrue / s                                   # 參考母體的真實死亡率
      set.seed(SEED + 11000 * ni + 97 * si + r)
      Dref <- matrix(rpois(length(Eref), Eref * mR), A, dimnames = dimnames(Eref))
      MRhat <- pmax(Dref, 0.5) / Eref                   # 參考母體的觀測死亡率
      for (nm in nmS) {
        j <- j + 1L
        f <- tryCatch(switch(nm,
            "Partial SMR"     = lc_psmr(D, E, MRhat),
            "Whittaker ratio" = lc_whit_ratio(D, E, MRhat),
            "零格以 SMR 填補" = lc_zero_smr(D, E, MRhat)),
          error = function(e) NULL)
        if (!is.null(f) && all(is.finite(f$a)) && all(is.finite(f$b))) {
          aM[r, , j] <- f$a - truth$a
          bS[r, j]   <- sum((f$b - truth$b)^2) / sum(truth$b^2)
          mh <- exp(outer(f$a, rep(1, Tn)) + outer(f$b, f$k))
          eS[r, j] <- tryCatch(life_table(mh[, Tn])$e0, error = function(e) NA_real_)
        }
      }
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }

  res[[ni]] <- data.frame(
    N = N, 估計量 = lab,
    alpha中位 = round(sapply(seq_len(K), function(j)
      median(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    alpha最大 = round(sapply(seq_len(K), function(j)
      max(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    SSE_beta  = round(apply(bS, 2, median, na.rm = TRUE), 4),
    e0偏誤    = round(apply(eS, 2, median, na.rm = TRUE) - e0_true, 3),
    e0標準差  = round(apply(eS, 2, sd, na.rm = TRUE), 3),
    e0四分位距 = round(apply(eS, 2, function(v) IQR(v, na.rm = TRUE)), 3),
    stringsAsFactors = FALSE)
  cat(sprintf("\n===== N = %.0e =====\n", N))
  print(res[[ni]][, -1], row.names = FALSE)
}
tab <- do.call(rbind, res)
dir.create("output/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(tab, "output/tables/tableM1_smr_compare.csv", row.names = FALSE)
cat("\n已輸出 tableM1_smr_compare.csv\n")
