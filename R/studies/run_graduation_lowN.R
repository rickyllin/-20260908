###############################################################################
# 修勻能否把中心化校正的適用下界推低（分支研究：修勻 × 貝氏混合）
#
#   規模極限_待討論.md 記下兩件事：N <= 5000 時整列全零的年齡開始出現
#   （5000 時平均 1.9 個、2000 時 4.3 個），校正達成率掉到 72% 與 55%，
#   而 Firth 反而較穩，因為它不需要先估 mu_hat。該檔第五節問「下界能否
#   藉混合規則推低」。tableG2 顯示輸送版的 Y-prime 在 N=1e4 把 e0 偏誤由
#   -1.15 壓到 -0.00，故真正的候選答案是修勻而非切換到 Firth：
#   整列全零時，修勻由鄰近年齡借到的資訊正是 c0 所無法提供的。
#   本程式往下掃到 N=2000。
#
#   tableG1 顯示「先修勻、再施以中心化」的 alpha 最大偏誤隨平滑量單調爆開
#   （N=5e4：lambda=0.3 時 0.672、1 時 2.18、3 時 5.17、10 時 11.07）。
#   成因可直接由代數讀出，不必猜。
#
#   b(mu; c0) 是一個特定泛函的偏誤，即 E[log max(D, c0)] - log mu。
#   修勻之後的資料是 ltilde = S l，其中 S = (W + h Delta' Delta)^{-1} W，
#   故
#        E[ltilde] = S log mu + S b(mu)
#   而目標仍是 log mu。於是修勻後的偏誤是
#        (S - I) log mu   +   S b(mu)
#   第一項是修勻本身的偏誤，第二項是被修勻算子重新分配過的 b。
#   **修勻改變了估計量，b 因而不再是它的偏誤**；原樣扣 b 是誤設，
#   且誤設量隨 h 增大，正是實測所見的單調爆開。
#
#   本程式檢驗三個版本：
#     乙   原樣扣 b（tableG1 的作法，作為對照）
#     乙'  扣 S b(mu_hat)                       只輸送，不含修勻自身的偏誤
#     乙'' 扣 S b(mu_hat) + (S - I) log mu_hat  完整的兩項
#
# 輸出：output/tables/tableG3_graduation_lowN.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/partial_smr.R")

SEED <- 20261006; REPS <- 100; NS <- c(2e3, 5e3, 1e4); YEARS <- 2001:2024
LAM  <- c(1, 3, 10, 30); PEN_FROM <- 3L

dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0)
truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn]); e0_true <- life_table(mtrue[, Tn])$e0
BF <- make_logbias(0.5)

#' 逐年的修勻算子 S_j（A x A）。前兩個年齡組不受懲罰，故該區塊為單位矩陣。
smoother_list <- function(W, lam, z = 2, from = PEN_FROM) {
  A <- nrow(W); idx <- from:A
  Dm <- diag(length(idx)); for (i in seq_len(z)) Dm <- diff(Dm)
  lapply(seq_len(ncol(W)), function(j) {
    w <- W[idx, j]
    S <- diag(A)
    S[idx, idx] <- solve(diag(w) + lam * mean(w) * crossprod(Dm), diag(w))
    S
  })
}
apply_S <- function(SL, M) {
  out <- M
  for (j in seq_len(ncol(M))) out[, j] <- SL[[j]] %*% M[, j]
  out
}

#' 修勻後的中心化。mode 決定扣除項：
#'   "raw"       扣 b(mu_hat)
#'   "transport" 扣 S b(mu_hat)
#'   "full"      扣 S b(mu_hat) + (S - I) log mu_hat
lc_grad_center <- function(lm_, E, SL, bfun, mode, maxit = 20, tol = 1e-9) {
  Tn <- ncol(lm_); lt <- apply_S(SL, lm_)
  f <- lc_svd_fit(lt); a <- f$a; b <- f$b; k <- f$k
  for (it in seq_len(maxit)) {
    mu_hat <- E * exp(outer(a, rep(1, Tn)) + outer(b, k))
    Bm <- bfun(mu_hat)
    corr <- switch(mode,
      "raw"       = Bm,
      "transport" = apply_S(SL, Bm),
      "full"      = apply_S(SL, Bm) + apply_S(SL, log(mu_hat)) - log(mu_hat))
    bbar <- rowMeans(corr)
    W  <- mu_hat / mean(mu_hat)
    am <- rowSums(W * lt) / rowSums(W)
    Z  <- (lt - am) * sqrt(W)
    sv <- svd(Z); u1 <- sv$u[, 1]; v1 <- sv$v[, 1]
    if (sum(u1) < 0) { u1 <- -u1; v1 <- -v1 }
    b2 <- u1 / sum(u1)
    k2 <- sv$d[1] * v1 * sum(u1) / sqrt(pmax(colMeans(W), 1e-12))
    k2 <- k2 - mean(k2); sb <- sum(b2); b2 <- b2 / sb; k2 <- k2 * sb
    a2 <- am - bbar
    dif <- max(abs(a2 - a), abs(b2 - b), abs(k2 - k) / max(1, max(abs(k))))
    a <- a2; b <- b2; k <- k2
    if (dif < tol) break
  }
  list(a = a, b = b, k = k, iter = it)
}

MODES <- c("乙' 扣 S·b" = "transport")  #

lab <- c("標準 Lee-Carter", "Firth", "中心化＋加權",
         as.vector(t(outer(names(MODES), sprintf(" λ=%g", LAM), paste0))))
K <- length(lab)

res <- list()
for (ni in seq_along(NS)) {
  N <- NS[ni]
  E  <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu <- E * mtrue
  aM <- array(NA_real_, c(REPS, A, K)); bS <- eS <- matrix(NA_real_, REPS, K)

  for (r in seq_len(REPS)) {
    set.seed(SEED + 13000 * ni + r)                  # 與 tableG1 同一組樣本
    D   <- matrix(rpois(length(E), mu), A, dimnames = dimnames(E))
    lm_ <- log(pmax(D, 0.5) / E)
    f0  <- lc_svd_fit(lm_)
    Wmu <- E * exp(outer(f0$a, rep(1, Tn)) + outer(f0$b, f0$k))

    fits <- vector("list", K); fits[[1]] <- f0
    fits[[2]] <- tryCatch(lc_poisson_firth(D, E, firth = TRUE), error = function(e) NULL)
    fits[[3]] <- tryCatch(lc_alpha_only(D, E, bfun = BF, weight = TRUE),
                          error = function(e) NULL)
    j <- 3L
    for (md in MODES) {
      for (lam in LAM) { j <- j + 1L
        SL <- smoother_list(Wmu, lam)
        fits[[j]] <- tryCatch(lc_grad_center(lm_, E, SL, BF, md),
                              error = function(e) NULL) }
    }
    ## 重排為 label 的順序（label 是 mode 內層、lambda 外層交錯）
    for (j in seq_len(K)) {
      f <- fits[[j]]
      if (is.null(f) || !all(is.finite(f$a)) || !all(is.finite(f$b))) next
      aM[r, , j] <- f$a - truth$a
      bS[r, j]   <- sum((f$b - truth$b)^2) / sum(truth$b^2)
      mh <- exp(outer(f$a, rep(1, Tn)) + outer(f$b, f$k))
      eS[r, j] <- tryCatch(life_table(mh[, Tn])$e0, error = function(e) NA_real_)
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }

  lab2 <- c(lab[1:3], as.vector(outer(sprintf(" λ=%g", LAM), names(MODES),
                                      function(a, b) paste0(b, a))))
  res[[ni]] <- data.frame(
    N = N, 估計量 = lab2,
    alpha中位 = round(sapply(seq_len(K), function(j)
      median(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    alpha最大 = round(sapply(seq_len(K), function(j)
      max(abs(apply(aM[, , j], 2, median, na.rm = TRUE)))), 4),
    SSE_beta  = round(apply(bS, 2, median, na.rm = TRUE), 4),
    e0偏誤    = round(apply(eS, 2, median, na.rm = TRUE) - e0_true, 3),
    e0標準差  = round(apply(eS, 2, sd, na.rm = TRUE), 3),
    e0均方根  = round(sqrt(apply(eS, 2, function(v)
      mean((v - e0_true)^2, na.rm = TRUE))), 3),
    stringsAsFactors = FALSE)
  cat(sprintf("\n===== N = %.0e =====\n", N))
  print(res[[ni]][, -1], row.names = FALSE)
}
tab <- do.call(rbind, res)
write.csv(tab, "output/tables/tableG3_graduation_lowN.csv", row.names = FALSE)
cat("\n已輸出 tableG3_graduation_lowN.csv\n")
