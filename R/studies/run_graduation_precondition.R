###############################################################################
# 修勻的事前先決條件診斷（進度報告_1012 第陸節）
#
#   命題（報告第陸節之一）。平滑後的均方誤差為
#       MSE(h) = (c'(S_h - I)theta + c'S_h delta)^2 + c'S_h V S_h' c，
#   於 h = 0 微分，若 c'delta = 0 則第一項消失而第二項為負，故存在 h > 0
#   使均方誤差嚴格改善；若 c'delta != 0 則第一項為一階，改善沒有保證。
#   亦即修勻的改善只在它所施加的估計量對目標泛函已經無偏時才有保證，
#   而中心化校正所做的正是把 delta 移除——先校正、後修勻。
#
#   命題（報告第陸節之二）。(S_h - I)theta = 0 對一切 h 成立，若且唯若
#   theta 落在差分懲罰的零空間，即格子指標的 z-1 次以下多項式。z = 2 時
#   零空間為對數死亡率在格子上成一次式：沿年齡即 Gompertz 法則、
#   沿時間即固定漂移。先決條件因而是「懲罰的零空間須包含所相信的模型」，
#   而本文兩個方向的結構假設恰好都是 z = 2 的零空間。
#
#   本程式算式 (ratio)：逐年齡的 MSE 比。每一項都有閉式，所需輸入只有
#   曝露數與一組粗略死亡率，與現行的事前可行性檢定相同。
#   V_t 為對角，蓋因各格死亡數獨立，故無須估計共變。
#
#   注意：不可用一階近似 S = I - h W^{-1} D'D。資訊權重跨四個數量級，
#   單一 h 使幼年組落在幾乎完全平滑的區域而非擾動區域，一階展開在那裡
#   給出 10^3 量級的荒謬值。此處直接用 S 的定義計算。
#
# 輸出：終端表格（報告表 7、表 8）
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/analytic_variance.R")
c0 <- 0.5; PEN <- 3L
dat <- load_data(sex="Female"); keep <- which(dat$years %in% 2001:2024)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]
A <- nrow(D0); Tn <- ncol(D0); ages <- rownames(D0)
truth <- lc_poisson_firth(D0,E0,firth=FALSE)
mtrue <- exp(outer(truth$a,rep(1,Tn))+outer(truth$b,truth$k))
wage <- E0[,Tn]/sum(E0[,Tn]); lmt <- log(mtrue)
idx <- PEN:A; Dm <- diag(length(idx)); for(i in 1:2) Dm <- diff(Dm); P <- crossprod(Dm)
Smat <- function(w, lam) { S <- diag(A)
  S[idx,idx] <- solve(diag(w[idx]) + lam*mean(w[idx])*P, diag(w[idx])); S }

cat("### 事前診斷：修勻對 alpha_hat 的 MSE 是改善還是惡化\n")
cat("    均為閉式，只需曝露數與一組粗略死亡率，與現行事前可行性檢定同樣的輸入。\n")
for (N in c(1e4, 5e4)) {
  E <- outer(N*wage, rep(1,Tn)); MU <- E*mtrue
  vv <- matrix(logvar_exact(as.vector(MU), c0), A)
  var_raw <- rowSums(vv)/Tn^2
  lam <- if (N==1e4) 3 else 1
  bias <- numeric(A); var_sm <- numeric(A)
  for (t in seq_len(Tn)) {
    S <- Smat(MU[,t], lam)
    bias <- bias + ((S - diag(A)) %*% lmt[,t])/Tn
    var_sm <- var_sm + diag(S %*% diag(vv[,t]) %*% t(S))/Tn^2
  }
  bias <- as.vector(bias)
  mse_raw <- var_raw; mse_sm <- bias^2 + var_sm
  cat(sprintf("\n  N=%.0e，lam=%g（實測較佳值）\n", N, lam))
  cat(sprintf("  %-7s %9s %9s %9s %9s %8s %s\n","年齡","平滑偏誤","sd 原","sd 修勻",
              "√MSE 比","判讀","|偏誤|/sd原"))
  for (x in seq_len(A)) {
    r <- sqrt(mse_sm[x]/mse_raw[x])
    cat(sprintf("  %-7s %+9.4f %9.4f %9.4f %9.3f %8s %6.2f\n", ages[x], bias[x],
      sqrt(var_raw[x]), sqrt(var_sm[x]), r,
      ifelse(x < PEN, "區外", ifelse(r < 0.9, "改善", ifelse(r <= 1.05, "持平", "惡化"))),
      abs(bias[x])/sqrt(var_raw[x])))
  }
  ok <- idx[sqrt(mse_sm[idx]/mse_raw[idx]) < 0.9]
  cat(sprintf("  → 改善的年齡：%s\n", paste(ages[ok], collapse=", ")))
  cat(sprintf("  → 全體 √(ΣMSE) 比：%.4f\n", sqrt(sum(mse_sm)/sum(mse_raw))))
}
