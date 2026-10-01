###############################################################################
# 逐年齡的改善剖面（回應 0929 會議：alpha, beta 的改善情況可以更詳細，
# 由整體細到逐年齡，看是否有特定年齡或年齡段改善）
#
#   本檔把第伍節的整體指標拆到年齡維度。拆解的理由不只是呈現：閉式
#   b(mu; c0) 是期望死亡數的函數，而期望死亡數在年齡間相差數個數量級，
#   故閉式預測改善應\textbf{集中在死亡數少的年齡}，在死亡數多的年齡幾乎
#   沒有改善空間。逐年齡的結果因此是對推導的另一次檢驗，而不只是細節。
#
# 輸出：output/tables/tableM4_age_profile.csv（逐年齡）
#       output/tables/tableM4_age_band.csv（年齡段彙總）
#       output/figures/figT9_age_profile.png
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/fig_axis_utils.R")

SEED <- 20260926; REPS <- 100; NS <- c(1e4, 5e4); YEARS <- 2001:2024
dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0); ages <- rownames(D0)
truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage <- E0[, Tn] / sum(E0[, Tn]); BF <- make_logbias(0.5)

FIT <- list("標準 LC"      = function(D, E) lc_svd_fit(log(pmax(D, 0.5) / E)),
            "Firth"        = function(D, E) lc_poisson_firth(D, E, firth = TRUE),
            "中心化"       = function(D, E) lc_alpha_only(D, E, bfun = BF, weight = FALSE),
            "中心化＋加權" = function(D, E) lc_alpha_only(D, E, bfun = BF, weight = TRUE))
K <- length(FIT)
band <- cut(seq_len(A), breaks = c(0, 4, 8, 13, 18, A),
            labels = c("0-14 歲", "15-34 歲", "35-59 歲", "60-84 歲", "85 歲以上"))

prof <- list(); bnd <- list(); PR <- vector("list", length(NS))
for (ni in seq_along(NS)) {
  N <- NS[ni]; E <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu <- E * mtrue
  aM <- bM <- array(NA_real_, c(REPS, A, K))
  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000 * ni + r)                # 與第伍節同一組種子
    D <- matrix(rpois(length(E), mu), A, dimnames = dimnames(E))
    for (j in seq_len(K)) {
      f <- tryCatch(FIT[[j]](D, E), error = function(e) NULL)
      if (is.null(f) || !all(is.finite(f$a)) || !all(is.finite(f$b))) next
      aM[r, , j] <- f$a - truth$a
      bM[r, , j] <- f$b - truth$b
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }
  rmse <- function(X) sqrt(colMeans(X^2, na.rm = TRUE))
  med  <- function(X) apply(X, 2, median, na.rm = TRUE)
  aR <- sapply(seq_len(K), function(j) rmse(aM[, , j]))
  bR <- sapply(seq_len(K), function(j) rmse(bM[, , j]))
  aB <- sapply(seq_len(K), function(j) med(aM[, , j]))
  bB <- sapply(seq_len(K), function(j) med(bM[, , j]))
  mubar <- rowMeans(mu)
  PR[[ni]] <- list(aR = aR, bR = bR, mubar = mubar)

  prof[[ni]] <- data.frame(
    N = N, 年齡 = rep(ages, K), 年齡段 = rep(as.character(band), K),
    平均期望死亡數 = round(rep(mubar, K), 2),
    估計量 = rep(names(FIT), each = A),
    alpha偏誤 = round(as.vector(aB), 4), alpha均方根 = round(as.vector(aR), 4),
    beta偏誤  = round(as.vector(bB), 4), beta均方根  = round(as.vector(bR), 4),
    stringsAsFactors = FALSE)

  agg <- function(M) tapply(M, band, mean)
  bnd[[ni]] <- data.frame(
    N = N, 年齡段 = levels(band),
    平均期望死亡數 = round(as.vector(tapply(mubar, band, mean)), 2),
    alpha_標準LC = round(as.vector(agg(aR[, 1])), 4),
    alpha_Firth  = round(as.vector(agg(aR[, 2])), 4),
    alpha_中心化 = round(as.vector(agg(aR[, 3])), 4),
    alpha_中心化加權 = round(as.vector(agg(aR[, 4])), 4),
    alpha_改善倍數 = round(as.vector(agg(aR[, 1]) / agg(aR[, 4])), 2),
    beta_標準LC = round(as.vector(agg(bR[, 1])), 4),
    beta_Firth  = round(as.vector(agg(bR[, 2])), 4),
    beta_中心化 = round(as.vector(agg(bR[, 3])), 4),
    beta_中心化加權 = round(as.vector(agg(bR[, 4])), 4),
    beta_改善倍數 = round(as.vector(agg(bR[, 1]) / agg(bR[, 4])), 2),
    stringsAsFactors = FALSE)
  cat(sprintf("\n===== N = %.0e =====\n", N)); print(bnd[[ni]][, -1], row.names = FALSE)
}
write.csv(do.call(rbind, prof), "output/tables/tableM4_age_profile.csv", row.names = FALSE)
write.csv(do.call(rbind, bnd),  "output/tables/tableM4_age_band.csv",    row.names = FALSE)

COL <- c("#B03A2E", "#7D3C98", "#1F6FB2", "#1E8449"); PCH <- c(16, 17, 18, 15)
png_cjk("output/figures/figT9_age_profile.png", width = 2250, height = 1250, res = 160)
par(mfrow = c(2, 2), mar = c(5.2, 4.8, 3.2, 1.0), mgp = c(2.9, 0.8, 0))
for (ni in seq_along(NS)) for (which_p in c("alpha", "beta")) {
  M <- if (which_p == "alpha") PR[[ni]]$aR else PR[[ni]]$bR
  xx <- seq_len(A)
  plot(xx, M[, 1], type = "n", log = "y", xaxt = "n", xlab = "年齡組",
       ylab = if (which_p == "alpha") expression(alpha[x]~"的均方根誤差")
              else expression(beta[x]~"的均方根誤差"),
       ylim = range(M[M > 0], na.rm = TRUE),
       main = sprintf("%s：N = %s",
                      if (which_p == "alpha") "年齡水準" else "年齡敏感度",
                      formatC(NS[ni], format = "d", big.mark = ",")))
  grid(col = "grey88", lty = 1)
  for (j in seq_len(K))
    lines(xx, M[, j], col = COL[j], lwd = 2.4, type = "b", pch = PCH[j], cex = 0.8)
  axis(1, at = xx, labels = ages, las = 2, cex.axis = 0.78)
  if (ni == 1 && which_p == "alpha")
    legend("topright", bty = "n", cex = 0.92, legend = names(FIT),
           col = COL, lwd = 2.4, pch = PCH)
}
dev.off()
cat("\n已輸出 tableM4 與 figT9\n")
