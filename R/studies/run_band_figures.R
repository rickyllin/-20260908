###############################################################################
# 估計值的逐年齡分布圖：帶寬與真值的相對位置
#
#   目的。第肆節的分解指出 alpha 的偏誤幾乎全為確定性、beta 的偏誤九成以上
#   為隨機。這兩件事在圖上有截然不同的樣子，而圖比表更能說明：
#     alpha：帶\textbf{窄}但整條帶\textbf{偏離}真值（偏誤主導）
#     beta ：帶\textbf{寬}且大致罩住真值（變異數主導）
#   圖中另把閉式預測 bbar_x = mean_t b(mu_xt) 疊在 alpha 的帶上，
#   若推導正確，該線應落在帶的中線上。
#
# 輸出：output/figures/figT1_alpha_band.png
#       output/figures/figT2_beta_band.png
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/fig_axis_utils.R")

SEED <- 20260926; REPS <- 100; NS <- c(1e4, 5e4); YEARS <- 2001:2024
dat <- load_data(sex="Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]; A<-nrow(D0); Tn<-ncol(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage <- E0[,Tn]/sum(E0[,Tn]); ages <- rownames(D0)
BF <- make_logbias(0.5)

fit_w <- function(D,E) {                      # 中心化 + 資訊加權
  g <- lc_analytic(D,E,bfun=BF,shrink=1)
  mu <- E*exp(outer(g$a,rep(1,ncol(D)))+outer(g$b,g$k))
  lmc <- log(pmax(D,0.5)/E) - BF(mu); W <- mu/mean(mu)
  am <- rowSums(W*lmc)/rowSums(W); Z <- (lmc-am)*sqrt(W)
  sv <- svd(Z); u1<-sv$u[,1]; v1<-sv$v[,1]
  if (sum(u1)<0){u1<--u1; v1<--v1}
  bb<-u1/sum(u1); kk<-sv$d[1]*v1*sum(u1)/sqrt(pmax(colMeans(W),1e-12))
  kk<-kk-mean(kk); sb<-sum(bb); list(a=am,b=bb/sb,k=kk*sb)
}
SPEC <- list(
  list(id="標準 LC",       col="#B03A2E", f=function(D,E) lc_svd_fit(log(pmax(D,0.5)/E))),
  list(id="中心化",        col="#1F6FB2", f=function(D,E) lc_analytic(D,E,bfun=BF,shrink=1)),
  list(id="中心化＋加權",  col="#1E8449", f=fit_w),
  list(id="Firth",         col="#7D3C98", f=function(D,E) lc_poisson_firth(D,E,firth=TRUE))
)
J <- length(SPEC)

## ---- 模擬：保留逐年齡逐重複的估計值 ----
AA <- BB <- vector("list", length(NS)); PRED <- vector("list", length(NS))
for (i in seq_along(NS)) {
  N <- NS[i]
  E <- outer(N*wage, rep(1,Tn)); dimnames(E) <- dimnames(D0); mu <- E*mtrue
  PRED[[i]] <- rowMeans(BF(mu))                                  # 閉式預測 bbar_x
  aA <- array(NA_real_, c(REPS,A,J)); bA <- array(NA_real_, c(REPS,A,J))
  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000*i + r)
    D <- matrix(rpois(length(E), mu), A, dimnames=dimnames(E))
    for (j in seq_len(J)) {
      f <- tryCatch(SPEC[[j]]$f(D,E), error=function(e) NULL)
      if (is.null(f) || !all(is.finite(f$a)) || !all(is.finite(f$b))) next
      aA[r,,j] <- f$a - truth$a; bA[r,,j] <- f$b
    }
    if (r%%25==0) cat(sprintf("  N=%.0e rep %d/%d\n", N, r, REPS))
  }
  AA[[i]] <- aA; BB[[i]] <- bA
}

band <- function(M, q=c(0.05,0.95)) {      # M: REPS x A
  list(lo=apply(M,2,quantile,q[1],na.rm=TRUE),
       md=apply(M,2,median,na.rm=TRUE),
       hi=apply(M,2,quantile,q[2],na.rm=TRUE))
}
tp <- function(col, alpha = 0.18) {        # 半透明色
  v <- col2rgb(col) / 255
  grDevices::rgb(v[1], v[2], v[3], alpha)
}
xx <- seq_len(A)

## =================== 圖一：alpha 的偏誤帶 ===================
png_cjk("output/figures/figT1_alpha_band.png", width=2100, height=950, res = 211)
par(mfrow=c(1,2), mar=c(5.2,4.4,3.2,1.0), mgp=c(2.6,0.8,0))
for (i in seq_along(NS)) {
  use <- c(1,2)                                    # 標準 LC 與 中心化
  bd  <- lapply(use, function(j) band(AA[[i]][,,j]))
  ## ylim 以各帶的 2%/98% 分位為界並留邊，避免極端年齡把全圖壓扁
  allv <- unlist(lapply(bd, function(b) c(b$lo, b$hi)))
  yl <- range(c(quantile(allv, c(0.02, 0.98), na.rm = TRUE), 0, PRED[[i]]),
              na.rm = TRUE)
  yl <- yl + c(-1, 1) * 0.10 * diff(yl)
  plot(NA, xlim=c(1,A), ylim=yl, xaxt="n", xlab="年齡組",
       ylab=expression(hat(alpha)[x] - alpha[x]),
       main=sprintf("N = %s", format(NS[i], big.mark=",", scientific=FALSE)))
  axis(1, at=xx, labels=ages, las=2, cex.axis=0.62)
  abline(h=0, col="grey40", lwd=1.1)
  abline(v=xx, col="grey92", lty=3)
  for (m in seq_along(use)) {
    j <- use[m]; b <- bd[[m]]; cl <- SPEC[[j]]$col
    polygon(c(xx,rev(xx)), c(b$lo,rev(b$hi)), col=tp(cl), border=NA)
    ## 虛線標出區間的實際位置（5% 與 95% 分位）
    lines(xx, b$lo, col=cl, lwd=1.1, lty=2)
    lines(xx, b$hi, col=cl, lwd=1.1, lty=2)
    lines(xx, b$md, col=cl, lwd=2.2)
  }
  lines(xx, PRED[[i]], col="black", lwd=2.0, lty=2)
  points(xx, PRED[[i]], pch=4, cex=0.7, col="black")
  legend("topright", bty="n", cex=0.74,
    legend=c("標準 LC：中線","標準 LC：5% 與 95%","中心化：中線",
             "中心化：5% 與 95%", expression(paste("閉式預測  ", bar(b)[x]))),
    col=c(SPEC[[1]]$col, SPEC[[1]]$col, SPEC[[2]]$col, SPEC[[2]]$col, "black"),
    lwd=c(2.2,1.1,2.2,1.1,2.0), lty=c(1,2,1,2,2))
}
dev.off()

## =================== 圖二：beta 的估計帶（每格一個估計量）===================
## 三條帶疊在同一格會糊成一片，故改為 2 列（人口規模）x 3 行（估計量），
## 每格只畫一條帶與真值，使「帶寬的寬窄」可以橫向直接比較。
png_cjk("output/figures/figT2_beta_band.png", width = 2250, height = 1250, res = 205)
par(mfrow = c(2, 3), mar = c(5.0, 4.2, 3.0, 0.8), mgp = c(2.5, 0.75, 0))
use <- c(1, 2, 3)
for (i in seq_along(NS)) {
  ## 同一列（同一人口規模）共用 ylim，使三個估計量的帶寬可橫向比較；
  ## 以該列三條帶的 3%/97% 分位為界，避免單一極端年齡把全列壓扁。
  rowv <- unlist(lapply(use, function(j) {
    b <- band(BB[[i]][, , j]); c(b$lo, b$hi) }))
  ## 取 12%/88% 分位而非全距：帶的極端處會把中線與真值壓扁，
  ## 而帶寬的確切數值已另以文字標出，故此處優先保留中線的可讀性。
  ylrow <- range(c(quantile(rowv, c(0.12, 0.88), na.rm = TRUE), truth$b),
                 na.rm = TRUE)
  ylrow <- ylrow + c(-1, 1) * 0.10 * diff(ylrow)
  for (m in seq_along(use)) {
    j  <- use[m]; b <- band(BB[[i]][, , j]); cl <- SPEC[[j]]$col
    yl <- ylrow
    plot(NA, xlim = c(1, A), ylim = yl, xaxt = "n", xlab = "年齡組",
         ylab = expression(hat(beta)[x]),
         main = sprintf("%s ．N = %s", SPEC[[j]]$id,
                        format(NS[i], big.mark = ",", scientific = FALSE)))
    axis(1, at = xx, labels = ages, las = 2, cex.axis = 0.58)
    abline(h = 0, col = "grey40"); abline(v = xx, col = "grey93", lty = 3)
    polygon(c(xx, rev(xx)),
            c(pmax(b$lo, yl[1]), rev(pmin(b$hi, yl[2]))),
            col = tp(cl, 0.22), border = NA)
    ## 虛線標出區間的實際位置；超出繪圖範圍者截在邊界上
    lines(xx, pmax(pmin(b$lo, yl[2]), yl[1]), col = cl, lwd = 1.1, lty = 2)
    lines(xx, pmax(pmin(b$hi, yl[2]), yl[1]), col = cl, lwd = 1.1, lty = 2)
    lines(xx, b$md, col = cl, lwd = 2.2)
    lines(xx, truth$b, col = "black", lwd = 2.4)
    points(xx, truth$b, pch = 16, cex = 0.6)
    legend("bottomleft", bty = "n", cex = 0.66,
           legend = c("真值", "估計中線", "5% 與 95%（虛線）"),
           col = c("black", cl, cl), lwd = c(2.4, 2.2, 1.1),
           lty = c(1, 1, 2))
    text(A, yl[2] - 0.04 * diff(yl), sprintf("帶寬中位 %.3f", median(b$hi - b$lo)),
         adj = c(1, 1), cex = 0.82, col = cl, font = 2)
  }
}
dev.off()

## ---- 帶寬的數值摘要，供文中引用 ----
cat("\n=== 帶寬（5–95% 區間長度）的中位數，以及中線離真值的距離 ===\n")
for (i in seq_along(NS)) {
  cat(sprintf("\n--- N = %.0e ---\n", NS[i]))
  for (j in seq_len(J)) {
    ba <- band(AA[[i]][,,j]); bb <- band(BB[[i]][,,j])
    cat(sprintf("%-14s alpha: 帶寬中位 %.4f  |中線| 中位 %.4f  最大 %.4f | beta: 帶寬中位 %.4f  |中線-真值| 中位 %.4f\n",
      SPEC[[j]]$id,
      median(ba$hi-ba$lo), median(abs(ba$md)), max(abs(ba$md)),
      median(bb$hi-bb$lo), median(abs(bb$md - truth$b))))
  }
  cat(sprintf("  閉式預測與標準 LC 中線的最大差 = %.4f\n",
      max(abs(PRED[[i]] - band(AA[[i]][,,1])$md))))
}
cat("\n已輸出 figT1_alpha_band.png / figT2_beta_band.png\n")
