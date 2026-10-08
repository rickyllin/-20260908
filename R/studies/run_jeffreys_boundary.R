###############################################################################
# 核對：Jeffreys 懲罰概似在 T >= X + 2 時無上界
# （1008策略文件_可行性評估.md 第一節；對應該文件二的定理 3.4）
#
#   該定理斷言沿 beta -> 無窮、kappa -> 0（eta 仍收斂）的退化路徑，
#       l_J(theta(s)) = l(eta_infty) + (T - X - 1) log s + O(1)，
#   故 T >= X + 2 時對任何資料全域最大值不存在。台灣五歲組為 X=22、T=24。
#
#   本程式自建該路徑 beta_x(s) = s * bt_x + 1/X（故 sum beta = 1 恆成立）、
#   kappa_t(s) = kt_t / s，並以 log det I 對 log s 的斜率核對。
#   切空間基底取 alpha 的 X 維、V_beta = {sum=0} 的 X-1 維、
#   V_kappa = {sum=0} 的 T-1 維，共 q = X + (X-1) + (T-1)。
#
#   核對結果：六組維度全部吻合（(22,24) 得 1.941 對理論 2），
#   且該路徑上對數概似收斂而 (1/2) log det I 持續上升，故合計無上界。
#
#   與本研究的關係：1012 V1 第柒節之二所證的是修正分數在全零列的 alpha_x
#   方向會變號、故有有限的「局部」根；此定理說的是「全域」上界。
#   兩者不矛盾，但 V1 的措辭須改（見評估第四節）。
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
## 沿退化路徑 beta(s) = s*bt + 1/X、kappa(s) = kt/s（alpha 固定）計算 log det I
## 切空間基底：alpha 的 X 維、V_beta = {sum=0} 的 X-1 維、V_kappa = {sum=0} 的 T-1 維
onb0 <- function(n) {            # {sum=0} 的正交單位基底
  M <- diag(n) - 1/n
  q <- qr(M); Q <- qr.Q(q)[, 1:(n-1), drop=FALSE]
  qr.Q(qr(Q))[, 1:(n-1), drop=FALSE]
}
logdetI <- function(a, b, k, E) {
  X <- length(a); T <- length(k); n <- X*T
  Bb <- onb0(X); Bk <- onb0(T)
  q <- X + (X-1) + (T-1)
  mu <- as.vector(E * exp(outer(a, rep(1,T)) + outer(b, k)))   # 欄優先 (x,t)
  J <- matrix(0, n, q)
  xi <- rep(seq_len(X), times = T); ti <- rep(seq_len(T), each = X)
  for (j in seq_len(X))        J[, j]            <- as.numeric(xi == j)
  for (j in seq_len(X-1))      J[, X+j]          <- Bb[xi, j] * k[ti]
  for (j in seq_len(T-1))      J[, X+(X-1)+j]    <- b[xi] * Bk[ti, j]
  Jw <- J * sqrt(mu)
  sv <- svd(Jw)$d
  list(ld = 2*sum(log(sv)), rank = sum(sv > max(sv)*1e-12), q = q)
}
run <- function(X, T, seed=1) {
  set.seed(seed)
  a <- seq(-9, -1, length.out = X)
  bt <- rnorm(X); bt <- bt - mean(bt); bt <- bt/sum(abs(bt))   # sum = 0
  kt <- rnorm(T); kt <- kt - mean(kt)
  E <- matrix(1e4/X, X, T)
  ss <- 10^seq(2, 4, length.out = 7)
  ld <- sapply(ss, function(s) logdetI(a, s*bt + 1/X, kt/s, E)$ld)
  sl <- coef(lm(ld ~ log(ss)))[2]
  r <- logdetI(a, ss[1]*bt+1/X, kt/ss[1], E)
  sprintf("(%2d,%2d)  數值斜率 %8.3f   理論 2(T-X-1) = %4d   秩 %d/%d",
          X, T, sl, 2*(T-X-1), r$rank, r$q)
}
cat("### 定理 3.4 的獨立核對：log det I 沿退化路徑對 log s 的斜率\n\n")
for (p in list(c(4,6), c(4,5), c(4,4), c(5,4), c(22,24), c(22,30))) cat(" ", run(p[1],p[2]), "\n")
cat("\n### 該路徑上概似是否真的收斂（故 log det 的發散不被概似抵銷）\n")
set.seed(7); X<-22; T<-24
a <- seq(-9,-1,length.out=X); bt<-rnorm(X); bt<-bt-mean(bt); bt<-bt/sum(abs(bt))
kt<-rnorm(T); kt<-kt-mean(kt); E<-matrix(1e4/X,X,T)
mu0 <- E*exp(outer(a,rep(1,T))+outer(bt,kt)); set.seed(9)
D <- matrix(rpois(length(E), mu0), X)
ll <- function(s){ b<-s*bt+1/X; k<-kt/s
  mu<-E*exp(outer(a,rep(1,T))+outer(b,k)); sum(D*log(mu)-mu) }
for (s in c(1,10,100,1000,10000)) cat(sprintf("  s=%7.0f  loglik %14.4f  (1/2)logdetI %10.4f  合計 %14.4f\n",
  s, ll(s), 0.5*logdetI(a, s*bt+1/X, kt/s, E)$ld,
  ll(s)+0.5*logdetI(a, s*bt+1/X, kt/s, E)$ld))
