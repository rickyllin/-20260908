###############################################################################
# L1 懲罰（trend filtering）下兩項不變性的核對（進度報告_1012 第肆節）
#
#   報告第參節的命題對任一凸的 J 成立，故 J = h||.||_1 是其特例，兩項
#   不變性與範數的選擇無關——不變性來自差分算子的零空間。本程式以 ADMM
#   解 L1 問題後核對，並比較同一粗糙度下 L2 與 L1 解的型態。
#
#   z = 1 時斜率不受保護，其機制由證明直接讀出：Delta^1 t = 1，故
#   t'(Delta^1)' g = sum_i g_i 不必為零。z = 2 時 Delta^2 t = 0，保證無條件成立。
#
#   解法器先以 h -> 無窮大核對：z=1 應退化為加權常數、z=2 應退化為加權
#   最小平方直線。（初版的 ADMM 曾把分裂變數初始化為 Delta y，使第一個
#   迭代恰等於 y 而停止準則立即成立，解形同未施懲罰；此處初始化為 0，
#   並改以原始殘差與對偶殘差同時收斂為準。）
#
# 輸出：終端表格（報告表 4、表 5）
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/analytic_bias_lc.R")
c0 <- 0.5
Dmat <- function(T,z){M<-diag(T); for(i in seq_len(z)) M<-diff(M); M}
#' trend filtering：min 1/2 (r-y)'W(r-y) + h||D r||_1，ADMM（修正初始化與停止準則）
tf <- function(y, w, h, z = 2, rho = NULL, maxit = 200000, tol = 1e-10) {
  T <- length(y); D <- Dmat(T,z); m <- nrow(D)
  if (is.null(rho)) rho <- max(mean(w), 1e-8)
  A <- solve(diag(w) + rho*crossprod(D))
  r <- y; v <- rep(0,m); u <- rep(0,m)            # v 為分裂變數，初始為 0
  for (k in seq_len(maxit)) {
    r <- as.vector(A %*% (w*y + rho*t(D)%*%(v-u)))
    Dr <- as.vector(D %*% r)
    vo <- v
    v <- sign(Dr+u)*pmax(abs(Dr+u)-h/rho, 0)
    u <- u + Dr - v
    if (max(abs(Dr-v)) < tol && max(abs(v-vo)) < tol) break
  }
  list(r=r, Dr=as.vector(D%*%r), kinks=which(abs(as.vector(D%*%r))>1e-5), iter=k)
}
cat("### 零、解法器核對（h 很大時 z=2 應退化為直線、z=1 應退化為常數）\n")
set.seed(7); T<-24; tt<-seq_len(T); y<-2*tt+rnorm(T,,3); w<-rep(1,T)
for (z in c(1,2)) { f<-tf(y,w,h=1e5,z=z)
  cat(sprintf("  z=%d h=1e5：折點 %d 個，與最小平方%s的最大差 %.2e\n", z, length(f$kinks),
    ifelse(z==1,"常數","直線"),
    max(abs(f$r - if(z==1) mean(y) else fitted(lm(y~tt)))))) }

cat("\n### 一、L1 的兩個不變性（KKT：1'D'=0 給均值；z=2 另有 t'D'=0 給斜率）\n")
set.seed(11); w<-runif(T,0.2,5)
for (nm in c("雜訊無趨勢","強單調上升")) {
  set.seed(11); y <- if(nm=="雜訊無趨勢") rnorm(T,,2)+0.3*tt else 2*tt+rnorm(T,,0.3)
  for (z in c(1,2)) {
    f<-tf(y,w,h=4,z=z); s0<-coef(lm(y~tt,weights=w))[2]; s1<-coef(lm(f$r~tt,weights=w))[2]
    cat(sprintf("  %s z=%d：加權均值差 %.1e｜斜率 %.6f→%.6f 差 %.2e｜Σsign(Dr)=%d｜折點%d\n",
      nm,z,abs(sum(w*f$r)-sum(w*y))/sum(w),s0,s1,abs(s1-s0),
      sum(sign(f$Dr)),length(f$kinks)))
  }
}
cat("\n### 二、臺灣女性 kappa：L1 與 L2 在同一粗糙度下的解的型態\n")
dat<-load_data(sex="Female"); keep<-which(dat$years %in% 2001:2024)
D0<-dat$D[,keep,drop=FALSE]; E0<-dat$E[,keep,drop=FALSE]
k<-lc_svd_fit(log(pmax(D0,c0)/E0))$k; yrs<-dat$years[keep]; Tk<-length(k)
rough<-function(v) sum(diff(diff(v))^2); P<-crossprod(Dmat(Tk,2))
l2<-function(h) as.vector(solve(diag(Tk)+h*P,k))
tgt<-rough(k)*0.5
hh<-10^seq(-4,4,length.out=800)
h2<-hh[which.min(abs(sapply(hh,function(h) rough(l2(h)))-tgt))]
hl<-hh[which.min(abs(sapply(hh,function(h) rough(tf(k,rep(1,Tk),h,2)$r))-tgt))]
cat(sprintf("  kappa 原始粗糙度 %.3f，兩法皆調到 %.3f\n", rough(k), tgt))
cat(sprintf("  L2 h=%.4f → %.3f，非零二階差分 %d/%d（全部非零，解處處彎曲）\n",
            h2, rough(l2(h2)), sum(abs(diff(diff(l2(h2))))>1e-8), Tk-2))
fl<-tf(k,rep(1,Tk),hl,2)
cat(sprintf("  L1 h=%.4f → %.3f，折點 %d/%d 個，位於 %s\n", hl, rough(fl$r),
            length(fl$kinks), Tk-2, paste(yrs[fl$kinks+1],collapse=", ")))
cat(sprintf("  原始 kappa 二階差分最大的五年：%s\n",
  paste(sprintf("%d(%+.2f)", yrs[order(-abs(diff(diff(k))))[1:5]+1],
        diff(diff(k))[order(-abs(diff(diff(k))))[1:5]]),collapse=", ")))
