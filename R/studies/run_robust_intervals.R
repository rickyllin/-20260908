###############################################################################
# 區間推論的第二來源：sigma 的估計，以及初始值的不確定性
#
#   報告 §捌之一 顯示現行作法的區間同時「太寬」且「涵蓋不足」：
#   寬度被灌水是因 kappa_hat 含估計雜訊使 sigma_hat 被高估，
#   涵蓋失敗則因漂移項被壓平使區間擺錯位置。
#
#   本腳本檢驗兩項針對性的處置：
#     R   以狀態空間 REML 分離變異成分（處理寬度）
#     RK  再加計 Kalman 濾波給出的初始值變異數 P_T（處理位置）
#
#   對照：O 不可約下限（參數已知）、P 現行作法（SVD 插入式）。
#
# 輸出：output/tables/tableL_sigma.csv     sigma 的估計
#       output/tables/tableM_interval.csv  區間寬度與涵蓋率
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
source("R/core/reml_kappa.R")

SEED <- 20260926; REPS <- 100; H <- 10; M <- 2000
PROBS <- c(0.10, 0.90)            # 名目 80%
NS <- c(5e4, 2e5, 1e6)
YEARS <- 2001:2024

dat  <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0)

truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn])
dk    <- diff(truth$k)
drift_true <- mean(dk); sigma_true <- sd(dk)
cat(sprintf("真值：drift = %.4f，sigma = %.4f\n", drift_true, sigma_true))

grid_n <- make_grid(A)$n
e0_of  <- function(m) life_table(m, grid_n)$e0

#' 由 (a, b, 初始值, 漂移, sigma) 模擬 h 期後的 e0 分布並取分位數
e0_band <- function(a, b, k_last, drift, sigma, jump_var = 0) {
  kh <- k_last + rowSums(matrix(rnorm(M * H, drift, sigma), M, H))
  if (jump_var > 0) kh <- kh + rnorm(M, 0, sqrt(jump_var))
  quantile(vapply(kh, function(kk) e0_of(exp(a + b * kk)), 0), PROBS, names = FALSE)
}

VARS <- c("O 不可約下限", "P 現行作法（SVD）",
          "R 卜瓦松 + 狀態空間 REML", "RK  R + Kalman 初始值變異數")

sig_rows <- list(); int_rows <- list()

for (N in NS) {
  E <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu_exp <- E * mtrue
  W <- C <- matrix(NA_real_, REPS, length(VARS))
  sg <- matrix(NA_real_, REPS, 3)   # naive / 狀態空間 ML / 狀態空間 REML

  for (r in seq_len(REPS)) {
    set.seed(SEED + 100 * which(NS == N) + r)
    ## 真實的未來路徑（所有變體共用，故為成對比較）
    k_fut  <- truth$k[Tn] + sum(rnorm(H, drift_true, sigma_true))
    target <- e0_of(exp(truth$a + truth$b * k_fut))

    D <- matrix(rpois(length(E), mu_exp), A, dimnames = dimnames(E))

    ## --- O：參數已知，只有未來創新項 ---
    bd <- e0_band(truth$a, truth$b, truth$k[Tn], drift_true, sigma_true)
    W[r, 1] <- diff(bd); C[r, 1] <- target >= bd[1] && target <= bd[2]

    ## --- P：SVD 插入式（現行作法）---
    fp <- lc_svd_fit(log(pmax(D, 0.5) / E))
    dp <- (fp$k[Tn] - fp$k[1]) / (Tn - 1)
    sp <- sd(diff(fp$k))
    bd <- e0_band(fp$a, fp$b, fp$k[Tn], dp, sp)
    W[r, 2] <- diff(bd); C[r, 2] <- target >= bd[1] && target <= bd[2]

    ## --- R / RK：卜瓦松（Firth）+ 狀態空間 ---
    ls <- tryCatch(lc_sigma(D, E, fitter = function(D, E)
                     lc_poisson_firth(D, E, firth = TRUE)),
                   error = function(e) NULL)
    if (!is.null(ls) && is.finite(ls$sigma)) {
      sg[r, ] <- c(ls$sigma_naive, NA, ls$sigma)
      ml <- tryCatch(sigma_state_space(ls$k, ls$tau^2, restricted = FALSE)$sigma,
                     error = function(e) NA_real_)
      sg[r, 2] <- ml
      bd <- e0_band(ls$a, ls$b, ls$k[Tn], ls$drift_ls, ls$sigma)
      W[r, 3] <- diff(bd); C[r, 3] <- target >= bd[1] && target <= bd[2]

      kj <- tryCatch(kalman_jumpoff(ls$k, ls$tau^2, ls$sigma, ls$drift_ls),
                     error = function(e) NULL)
      if (!is.null(kj)) {
        bd <- e0_band(ls$a, ls$b, kj$kappa_T, ls$drift_ls, ls$sigma,
                      jump_var = kj$var_T)
        W[r, 4] <- diff(bd); C[r, 4] <- target >= bd[1] && target <= bd[2]
      }
    }
    if (r %% 20 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }

  int_rows[[as.character(N)]] <- data.frame(
    N = N, 變體 = VARS,
    寬度 = round(apply(W, 2, median, na.rm = TRUE), 3),
    涵蓋率 = round(colMeans(C, na.rm = TRUE), 3),
    可行率 = round(colMeans(!is.na(W)), 2), stringsAsFactors = FALSE)
  sig_rows[[as.character(N)]] <- data.frame(
    N = N, 方法 = c("naive sd(diff(kappa_hat))", "狀態空間 ML", "狀態空間 REML"),
    中位數 = round(apply(sg, 2, median, na.rm = TRUE), 4),
    偏誤 = round(apply(sg, 2, median, na.rm = TRUE) - sigma_true, 4),
    IQR = round(apply(sg, 2, IQR, na.rm = TRUE), 3),
    可行率 = round(colMeans(!is.na(sg)), 2), stringsAsFactors = FALSE)

  cat(sprintf("\n===== N = %.0e =====\n", N))
  print(sig_rows[[as.character(N)]][, -1], row.names = FALSE)
  print(int_rows[[as.character(N)]][, -1], row.names = FALSE)
}

dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)
write.csv(do.call(rbind, sig_rows), "output/tables/tableL_sigma.csv", row.names = FALSE)
write.csv(do.call(rbind, int_rows), "output/tables/tableM_interval.csv", row.names = FALSE)
cat("\n已輸出 tableL / tableM\n")
