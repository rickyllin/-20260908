###############################################################################
# 不變性能否推廣到收縮先驗（進度報告_1012 第肆節之四，命題 5）
#
#   命題 2 要求懲罰為凸，而該假設可以去掉：論證其實只需要「先驗只透過
#   Delta^z r 依賴 r」。令 P 為對零空間 span{1,t} 的 W-正交投影，則先驗
#   在 N 的方向上為平坦，而高斯工作概似在 W-正交分解下可因式分解
#   （交叉項為零），故 Pr 的後驗為 N(Py, sigma^2 (P'WP)^{-1}) 且與其餘
#   獨立，於是後驗平均精確保留資料的零空間成分。
#
#   這覆蓋 L2、L1、NEG、horseshoe 與動態收縮——五者都只透過 Delta^z r
#   依賴 r。意涵是本文可採用 Faulkner-Minin (2018)、Kowal et al. (2019)
#   一路的任何先驗而不失去 alpha_hat 與漂移的保證。
#
#   本程式以 horseshoe 的 Gibbs 抽樣核對後驗平均（非凸，故不能用 KKT）。
#   資料刻意含一個轉折，以確保後驗真的在平滑而非恆等映射。
#
# 輸出：終端四行（報告第肆節之四所引的數字）
###############################################################################
# 命題二能否推廣到任何「只透過 Delta^z r 依賴 r」的先驗？
#   論證：令 P 為對 N = span{1,t} 的 W-正交投影。先驗 pi(r) = f(Delta^z r)
#   只依賴 (I-P)r，而高斯工作概似在 W-正交分解下可因式分解，
#   故 Pr 的後驗為 N(Py, (P'WP)^{-1})，與其餘獨立 => E[Pr|y] = Py 精確成立。
#   以 horseshoe（尺度混合，非凸）的 Gibbs 抽樣核對後驗平均。
set.seed(20261006)
T <- 24; tt <- seq_len(T); z <- 2
D <- diag(T); for (i in seq_len(z)) D <- diff(D); m <- nrow(D)
w <- runif(T, 0.2, 5); W <- diag(w)
y <- 2 + 0.3*tt + c(rep(0,12), rep(2.5,12)) + rnorm(T,,0.6)   # 含一個轉折
NIT <- 60000; BURN <- 10000
r <- y; lam <- rep(1,m); tau <- 1; sig2 <- 1
acc <- rep(0,T); n <- 0
for (it in seq_len(NIT)) {
  Pr <- crossprod(D, D/ (lam^2*tau^2)) / sig2                  # 先驗精確度
  A  <- W/sig2 + Pr
  Rc <- chol(A); mu <- backsolve(Rc, forwardsolve(t(Rc), (W %*% y)/sig2))
  r  <- as.vector(mu + backsolve(Rc, rnorm(T)))
  Dr <- as.vector(D %*% r)
  lam <- sqrt(1/rgamma(m, 1, 1 + Dr^2/(2*tau^2*sig2)))         # 近似半柯西
  tau <- sqrt(1/rgamma(1, (m+1)/2, 1 + sum(Dr^2/lam^2)/(2*sig2)))
  sig2 <- 1/rgamma(1, (T+m)/2, (sum(w*(y-r)^2) + sum(Dr^2/(lam^2*tau^2)))/2)
  if (it > BURN) { acc <- acc + r; n <- n + 1 }
}
rh <- acc/n
cm <- function(v) sum(w*v)/sum(w)
sl <- function(v) { tb <- cm(tt); sum(w*(tt-tb)*v)/sum(w*(tt-tb)^2) }
cat(sprintf("horseshoe 後驗平均（%d 次保留）\n", n))
cat(sprintf("  加權平均：資料 %.8f｜後驗 %.8f｜差 %.2e\n", cm(y), cm(rh), abs(cm(y)-cm(rh))))
cat(sprintf("  加權斜率：資料 %.8f｜後驗 %.8f｜差 %.2e\n", sl(y), sl(rh), abs(sl(y)-sl(rh))))
cat(sprintf("  （蒙地卡羅標準誤約 %.2e，故上列之差應落在同一量級）\n",
    sd(y-rh)/sqrt(n)))
cat(sprintf("  後驗確實平滑了資料：||r-y||=%.4f，Delta^2 的最大值落在第 %d 點\n",
    sqrt(sum((rh-y)^2)), which.max(abs(as.vector(D %*% rh)))))
