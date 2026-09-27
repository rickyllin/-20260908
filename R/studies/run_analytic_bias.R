###############################################################################
# 解析偏誤校正的實測：由推導得來的估計量能否勝過現有作法
#
#   受檢驗的假設。R/core/analytic_bias_lc.R 推導出對數轉換加零格替代的
#   逐格偏誤閉式 b(mu; c0)，並已驗證它能預測標準 LC 的逐年齡 alpha 偏誤
#   （相關 1.000 / 1.000 / 0.988）。若該預測正確，扣除 b(mu_hat) 應
#   大幅移除 alpha 的偏誤。
#
#   可事先寫下的三項預測：
#     1. alpha 的偏誤應大幅下降（因為它被完全預測）
#     2. 漂移項的偏誤\textbf{不應}大幅下降——因為第伍節已證明漂移項的
#        偏誤來自等權重分解對異質變異的敏感，與對數無關（變體 E）
#     3. 因 mu_hat 本身有誤差，完全扣除（shrink = 1）可能過頭，
#        故一併掃描 shrink
#
#   對照：標準 LC、卜瓦松 MLE、Firth。
#   種子與 run_eda_to_correction.R / run_mquantile_scan.R 完全相同。
#
# 輸出：output/tables/tableT_analytic_bias.csv
#       output/tables/tableT_analytic_by_age.csv
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
source("R/core/analytic_bias_lc.R")

SEED <- 20260926; REPS <- 100
NS <- c(1e4, 5e4, 2e5)
YEARS <- 2001:2024

dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0)
truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn])
drift_true <- (truth$k[Tn] - truth$k[1]) / (Tn - 1)
cat(sprintf("真值：drift = %.4f\n", drift_true))

BF <- make_logbias(0.5)

SPEC <- list(
  list(id = "A  標準 LC",                f = function(D,E) lc_svd_fit(log(pmax(D,0.5)/E))),
  list(id = "B  卜瓦松 MLE",             f = function(D,E) lc_poisson_firth(D,E,firth=FALSE)),
  list(id = "C  Firth",                  f = function(D,E) lc_poisson_firth(D,E,firth=TRUE)),
  list(id = "D  解析校正 shrink=1.0",    f = function(D,E) lc_analytic(D,E,bfun=BF,shrink=1.0)),
  list(id = "E  解析校正 shrink=0.8",    f = function(D,E) lc_analytic(D,E,bfun=BF,shrink=0.8)),
  list(id = "F  解析校正 shrink=0.6",    f = function(D,E) lc_analytic(D,E,bfun=BF,shrink=0.6)),
  list(id = "G  解析校正 + 加權 SVD",    f = function(D,E) {
        g <- lc_analytic(D,E,bfun=BF,shrink=1.0)
        mu <- E*exp(outer(g$a,rep(1,ncol(D)))+outer(g$b,g$k))
        lmc <- log(pmax(D,0.5)/E) - BF(mu)
        W <- mu/mean(mu)
        am <- rowSums(W*lmc)/rowSums(W); Z <- (lmc-am)*sqrt(W)
        sv <- svd(Z); u1<-sv$u[,1]; v1<-sv$v[,1]
        if (sum(u1)<0){u1<--u1; v1<--v1}
        bb <- u1/sum(u1); kk <- sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
        kk <- kk-mean(kk); sb<-sum(bb); list(a=am,b=bb/sb,k=kk*sb,iter=g$iter)})
)
ESTS <- vapply(SPEC, `[[`, character(1), "id"); J <- length(SPEC)
ok_fit <- function(f) all(is.finite(f$a)) && all(is.finite(f$b)) &&
                      all(is.finite(f$k)) && max(abs(f$k)) < 1e3

res <- list(); res_age <- list()
for (N in NS) {
  E <- outer(N * wage, rep(1, Tn)); dimnames(E) <- dimnames(D0)
  mu_exp <- E * mtrue
  aM <- array(NA_real_, c(REPS, A, J))
  bS <- dS <- cS <- matrix(NA_real_, REPS, J); div <- integer(J)
  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000 * which(NS == N) + r)
    D <- matrix(rpois(length(E), mu_exp), A, dimnames = dimnames(E))
    for (j in seq_len(J)) {
      f <- tryCatch(SPEC[[j]]$f(D, E), error = function(e) NULL)
      if (is.null(f) || !ok_fit(f)) { div[j] <- div[j] + 1L; next }
      aM[r, , j] <- f$a - truth$a
      bS[r, j] <- sum((f$b - truth$b)^2) / sum(truth$b^2)
      cS[r, j] <- if (sd(f$b) < 1e-12) NA_real_ else cor(f$b, truth$b)
      dS[r, j] <- (f$k[Tn] - f$k[1]) / (Tn - 1) - drift_true
    }
    if (r %% 25 == 0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }
  res[[as.character(N)]] <- data.frame(
    N = N, 估計量 = ESTS,
    alpha中位 = round(sapply(seq_len(J), function(j)
      median(abs(apply(aM[,,j],2,median,na.rm=TRUE)))), 4),
    alpha最大 = round(sapply(seq_len(J), function(j)
      max(abs(apply(aM[,,j],2,median,na.rm=TRUE)))), 4),
    SSE_beta = round(apply(bS,2,median,na.rm=TRUE),4),
    cor_beta = round(apply(cS,2,median,na.rm=TRUE),3),
    漂移偏誤 = round(apply(dS,2,median,na.rm=TRUE),4),
    發散率 = round(div/REPS,3), stringsAsFactors = FALSE)
  res_age[[as.character(N)]] <- data.frame(
    N=N, 年齡=rownames(D0), 期望死亡=round(rowMeans(mu_exp),2),
    setNames(as.data.frame(lapply(seq_len(J), function(j)
      round(apply(aM[,,j],2,median,na.rm=TRUE),4))), ESTS),
    check.names=FALSE, stringsAsFactors=FALSE)
  cat(sprintf("\n===== N = %.0e =====\n", N))
  print(res[[as.character(N)]][,-1], row.names=FALSE)
}
tab <- do.call(rbind,res)
dir.create("output/tables", recursive=TRUE, showWarnings=FALSE)
write.csv(tab,"output/tables/tableT_analytic_bias.csv",row.names=FALSE)
write.csv(do.call(rbind,res_age),"output/tables/tableT_analytic_by_age.csv",row.names=FALSE)
cat("\n=== 幼年組的 alpha 偏誤 ===\n")
for (N in NS) {
  d <- res_age[[as.character(N)]]
  print(d[d$年齡 %in% c("1-4","5-9","10-14","45-49","85-89"), -1], row.names=FALSE); cat("\n")
}
cat("已輸出 tableT\n")
