###############################################################################
# 不取對數時的估計：轉換的偏誤、秩一適足性與擬分數的對應
# （進度報告_理論初探_擬概似與變異安定轉換）
#
#   回應老師的問題：若不經過對數轉換，其他方法在小樣本下是否仍穩健。
#   把各法寫成同一族估計方程
#       sum w psi(h(D) - h(mu)) d h(mu)/d theta = 0，
#   問題即化為 (h, w, psi) 三個選擇。本程式算三組量。
#
#   一、各轉換的偏誤 E[h(D)] - h(mu)，折算為對數尺度（除以 h'(mu) mu）後比較。
#       恆等轉換一列全為零不是數值巧合，而是 E[D] - mu = 0（報告命題 1）。
#       對數轉換於 mu=0.02 達 +3.24 而平方根族停在約 -0.31。
#       Freeman-Tukey 的大值來自其慣用中心化 2 sqrt(mu+1/2) 是大 mu 的近似：
#       mu->0 時 E[h(D)] -> 1 而該式趨於 1.414，相差一個不隨 mu 消失的常數。
#
#   二、秩一結構在哪個尺度上成立。對數尺度最好（第一奇異值佔 83.9%、殘差
#       0.0657），平方根 75.7%/0.0840，率本身 72.1%/0.1356，m^(1/3) 居中。
#       殘差一律換算為「配適值取對數後與 log m 之差」的均方根，故可比較。
#
#   三、命題 2 的核對：對數尺度加權最小平方對應的變異函數為 mu^2/w，
#       故 w=mu 給 Poisson、w=1 給固定變異係數。於真值處比較分數向量
#       （而非參數估計值——各法識別條件不同，alpha 相差位移與尺度），
#       cor(U_{w=mu}, U_P) 由 0.588 升至 0.99998，而 cor(U_{w=1}, U_P)
#       停在 0.12 附近。
#
# 輸出：終端三組表格（報告表 2、表 4、表 3）
###############################################################################
## ---- 一、各轉換的偏誤 ----------------------------------------------------
c0 <- 0.5
Eh <- function(mu, h) sapply(mu, function(m){
  K <- max(80, ceiling(m + 12*sqrt(m) + 12)); d <- 0:K
  p <- dpois(d, m); sum(p*h(d))/sum(p) })
H <- list(
  "對數＋零格替代 log max(D,c0)" = list(h=function(d) log(pmax(d,c0)),   g=function(m) log(m)),
  "Anscombe 2sqrt(D+3/8)"       = list(h=function(d) 2*sqrt(d+3/8),     g=function(m) 2*sqrt(m+3/8)),
  "Freeman-Tukey sqrt(D)+sqrt(D+1)"= list(h=function(d) sqrt(d)+sqrt(d+1), g=function(m) 2*sqrt(m+1/2)),
  "Bartlett sqrt(D+3/8)"        = list(h=function(d) sqrt(d+3/8),       g=function(m) sqrt(m+3/8)),
  "恆等 D（不轉換）"             = list(h=function(d) d,                 g=function(m) m))
cat("### 一、各轉換的偏誤 E[h(D)] - h(mu)，以相對尺度比較\n")
cat("    相對尺度 = 偏誤 / h'(mu)，即折算回 log mu 的單位後才可互相比較\n\n")
mus <- c(0.02,0.05,0.1,0.25,0.5,1,2,5,20,100)
dh <- list(function(m) 1/m, function(m) 1/sqrt(m+3/8),
           function(m) 1/sqrt(m+1/2), function(m) 0.5/sqrt(m+3/8), function(m) 1)
cat(sprintf("%-32s %s\n","mu =", paste(sprintf("%8.2f",mus), collapse="")))
for (i in seq_along(H)) {
  b <- Eh(mus, H[[i]]$h) - H[[i]]$g(mus)
  rel <- b / dh[[i]](mus)
  cat(sprintf("%-32s %s\n", names(H)[i], paste(sprintf("%8.4f", rel), collapse="")))
}
cat("\n  （恆等轉換一列全為零，這不是數值巧合：E[D] - mu = 0 恆成立。）\n")
cat("\n### 二、各轉換偏誤的上界（折算回 log mu 單位後掃描 mu）\n")
mg <- exp(seq(log(0.01), log(500), length.out=3000))
for (i in seq_along(H)) {
  rel <- (Eh(mg, H[[i]]$h) - H[[i]]$g(mg)) / dh[[i]](mg)
  j <- which.max(abs(rel))
  cat(sprintf("  %-32s 最大 |偏誤| %8.4f 於 mu=%.3f｜mu->0 時 %s\n",
      names(H)[i], abs(rel[j]), mg[j],
      sprintf("%+.4f", rel[1])))
}

## ---- 二、秩一適足性 ------------------------------------------------------
dat <- load_data(sex="Female"); keep <- which(dat$years %in% 2001:2024)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]
Tn <- ncol(D0); A <- nrow(D0)
m <- D0/E0                                    # 全國觀測死亡率，抽樣誤差可忽略
cat("### 秩一結構在哪個尺度上成立：全國女性 2001-2024\n")
cat("    以列中心化後的奇異值衡量：第一個奇異值佔平方和的比例愈高，秩一愈好\n\n")
scales <- list(
  "log m（本文的模型尺度）" = log(m),
  "2 sqrt(m)（Anscombe 尺度）" = 2*sqrt(m),
  "sqrt(m)" = sqrt(m),
  "m（率本身）" = m,
  "m^(1/3)" = m^(1/3))
cat(sprintf("%-28s %10s %10s %12s\n","尺度","第一比例","前二比例","殘差 RMS"))
for (nm in names(scales)) {
  X <- scales[[nm]]; Xc <- X - rowMeans(X)
  sv <- svd(Xc)$d; p1 <- sv[1]^2/sum(sv^2); p2 <- sum(sv[1:2]^2)/sum(sv^2)
  f <- svd(Xc); r1 <- f$u[,1,drop=FALSE] %*% (f$d[1]*t(f$v[,1,drop=FALSE]))
  ## 殘差換算為對數尺度的 RMS，使各尺度可比
  fit <- switch(nm,
    "log m（本文的模型尺度）" = exp(rowMeans(X) + r1),
    "2 sqrt(m)（Anscombe 尺度）" = ((rowMeans(X) + r1)/2)^2,
    "sqrt(m)" = (rowMeans(X) + r1)^2,
    "m（率本身）" = rowMeans(X) + r1,
    "m^(1/3)" = (rowMeans(X) + r1)^3)
  rms <- sqrt(mean((log(pmax(fit,1e-12)) - log(m))^2))
  cat(sprintf("%-28s %10.5f %10.5f %12.5f\n", nm, p1, p2, rms))
}
cat("\n  殘差 RMS 為「配適值取對數後與 log m 之差」的均方根，故各尺度在同一單位下比較。\n")
cat("\n### 參照：Poisson 最大概似與 Firth 的既有結果（解析偏誤校正_推導與實測.md）\n")
cat("    人數一萬：標準 LC 的 alpha 最大偏誤 2.3309｜Poisson MLE 13.0317（完全分離而崩壞）\n")
cat("                Firth 0.4267｜本文中心化+加權 0.7399\n")
cat("    亦即「不轉換」使估計方程無偏，卻讓解本身在小樣本不存在或不穩定。\n")

## ---- 三、擬分數的對應 ----------------------------------------------------
