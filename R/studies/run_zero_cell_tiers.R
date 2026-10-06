###############################################################################
# 零格的分布與三個層級（進度報告_1012 第陸節之五）
#
#   一個零格所能借到的資訊有三個來源：秩一結構、鄰近年齡（修勻）、
#   參考母體。零格替代值 c0 是退化情形——它不向任何地方借，只填一個常數。
#   適用範圍由兩個量決定，而兩者都算得出來。
#
#     列資訊量 I_x = sum_t mu_{x,t}。Poisson 對 log mu 的 Fisher 資訊為 mu，
#       故 I_x 即該年齡全期對 alpha_x 的資訊總量。
#     局部曲度 |Delta^2 log mu_x|，即零空間條件（報告命題 4）。
#
#   甲級（I_x 充足）校正即足；乙級（I_x 不足而曲度小）可沿年齡修勻；
#   丙級（I_x 不足且曲度大）須借用外部資訊。
#
#   丙級恰為 0 至 20-24 的六組，與 run_graduation_precondition.R 由均方誤差
#   比判為惡化的年齡完全相同——兩者所用的量不同，結論一致，故該分界不是
#   任一條規則的產物。這同時解釋殘餘的 alpha 最大偏誤為何固定落在 1-4 歲
#   且不隨平滑量變動。
#
#   附帶結果：整列全零的期望年齡數 sum_x exp(-I_x) 亦為閉式，人數兩千、
#   五千、一萬時為 4.26、2.09、0.90，而規模極限_待討論.md 以模擬量到的
#   對應值為 4.3、1.9、1.0。該欄原以模擬取得，此處指出它有閉式。
#
# 輸出：終端表格（報告表 9）
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R"); source("R/core/analytic_variance.R")
c0 <- 0.5
dat <- load_data(sex="Female"); keep <- which(dat$years %in% 2001:2024)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]
A <- nrow(D0); Tn <- ncol(D0); ages <- rownames(D0)
tr <- lc_poisson_firth(D0,E0,firth=FALSE)
mt <- exp(outer(tr$a,rep(1,Tn))+outer(tr$b,tr$k)); wage <- E0[,Tn]/sum(E0[,Tn])
th <- rowMeans(log(mt)); d2 <- c(NA,NA,diff(diff(th)))

cat("### 零格的分布、列資訊量與曲度：三個層級的判準\n")
cat("    列資訊量 I_x = sum_t mu_{x,t}（該年齡各年期望死亡數之和，即 log mu 的 Fisher 資訊）\n")
cat("    整列全零機率 = exp(-I_x)；期望零格數 = sum_t exp(-mu_{x,t})（共 24 格）\n\n")
for (N in c(1e4, 5e3, 2e3)) {
  E <- outer(N*wage, rep(1,Tn)); MU <- E*mt
  I <- rowSums(MU); pz <- exp(-I); nz <- rowSums(exp(-MU))
  tier <- ifelse(I >= 20, "甲 校正即足",
           ifelse(!is.na(d2) & abs(d2) < 0.3, "乙 可沿年齡修勻", "丙 須參考母體"))
  cat(sprintf("  N=%.0e\n  %-7s %9s %9s %9s %9s  %s\n", N,
      "年齡","I_x","零格數","P(全零)","|Δ²logm|","層級"))
  for (x in seq_len(A)) cat(sprintf("  %-7s %9.2f %9.1f %9.2e %9s  %s\n", ages[x],
     I[x], nz[x], pz[x], ifelse(is.na(d2[x]),"—",sprintf("%.4f",abs(d2[x]))), tier[x]))
  cat(sprintf("  → 期望零格總數 %.1f / %d 格（%.1f%%）；整列全零的期望年齡數 %.2f\n",
      sum(nz), A*Tn, 100*sum(nz)/(A*Tn), sum(pz)))
  cat(sprintf("  → 丙級年齡：%s\n\n", paste(ages[tier=="丙 須參考母體"], collapse=", ")))
}
