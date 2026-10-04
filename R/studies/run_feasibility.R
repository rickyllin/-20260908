###############################################################################
# 事前可行性檢定的分級表（臺灣女性 2024 年齡結構）
#
#   程式：R/core/lc_feasibility.R。本檔只做兩件事：
#     1. 產生不同人口規模下的分級表（表 F1）；
#     2. 產生逐年齡的預測剖面（表 F2），供判讀哪些年齡拖累分級。
#   參考死亡率取全國配適值，與 1005 報告的真值相同，故兩者可對照。
#
# 輸出：output/tables/tableF1_feasibility_grade.csv
#       output/tables/tableF2_feasibility_byage.csv
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
source("R/core/analytic_bias_lc.R")
source("R/core/lc_feasibility.R")

YEARS <- 2001:2024; C0 <- 0.5
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0)
truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mref  <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn])

NS <- c(5e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e5, 1e6)
F1 <- NULL; F2 <- NULL
for (N in NS) {
  E <- N * outer(wage, rep(1, Tn)); rownames(E) <- rownames(D0)
  r <- lc_feasibility(E, mref, C0, ages = rownames(D0))
  F1 <- rbind(F1, cbind(N = N, r$summary[, -1]))
  F2 <- rbind(F2, cbind(N = N, r$byage))
}
write.csv(F1, "output/tables/tableF1_feasibility_grade.csv", row.names = FALSE)
write.csv(F2, "output/tables/tableF2_feasibility_byage.csv", row.names = FALSE)

cat("== 表 F1：分級表（c0 =", C0, "，mu* =", round(F1$mustar[1], 4), "）==\n")
print(data.frame(
  N              = format(F1$N, big.mark = ",", scientific = FALSE),
  `格子低於mu星` = sprintf("%.1f%%", 100 * F1$share_below_mustar),
  `預測alpha最大` = sprintf("%+.3f", F1$alpha_max_pred),
  `最差年齡`     = F1$worst_age,
  `正偏誤最大`   = sprintf("%+.3f", F1$b_pos_max),
  `負偏誤最小`   = sprintf("%+.3f", F1$b_neg_min),
  `主因`         = F1$driver,
  `低於ln2`      = sprintf("%.1f%%", 100 * F1$share_below_ln2),
  `分級`         = F1$grade,
  `e0的alpha通道`= sprintf("%+.3f", F1$e0_alpha_channel),
  check.names = FALSE), row.names = FALSE)

cat("\n== 表 F2：N = 50,000 的逐年齡剖面（前 8 列）==\n")
print(within(subset(F2, N == 5e4)[1:8, -1], {
  mu_median <- round(mu_median, 3); share_low <- round(share_low, 3)
  b_bar <- round(b_bar, 3); m_ratio <- round(m_ratio, 2)}), row.names = FALSE)

## ---- 圖 F1：逐年齡的預測偏誤剖面與分級門檻 ----------------------------
source("R/core/fig_axis_utils.R")
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
show <- c(1e4, 5e4, 2e5, 1e6)
cols <- c("#9e2a2b", "#e07a00", "#1f4e79", "#2d6a4f")
png_cjk("output/figures/figF1_feasibility.png", width = 2250, height = 900, res = 190)
par(mfrow = c(1, 2), mar = c(5.6, 4.4, 2.8, 0.8), mgp = c(2.7, 0.8, 0))

# (a) 逐年齡的預測偏誤
yl <- range(subset(F2, N %in% show)$b_bar)
plot(NA, xlim = c(1, A), ylim = yl + c(-0.15, 0.25), xaxt = "n",
     xlab = "", ylab = expression(bar(b)[x]~"（預測的"~alpha[x]~"偏誤）"),
     main = "(a) 逐年齡的預測偏誤剖面")
axis(1, at = seq(1, A, by = 3), labels = rownames(D0)[seq(1, A, by = 3)],
     las = 2, cex.axis = 0.85)
abline(h = 0, lty = 3, col = "grey50")
for (j in seq_along(show)) {
  d <- subset(F2, N == show[j])
  lines(seq_len(A), d$b_bar, lwd = 2.2, col = cols[j], type = "o", pch = 16, cex = 0.7)
}
legend("topright", legend = format(show, big.mark = ",", scientific = FALSE),
       col = cols, lwd = 2.2, pch = 16, bty = "n", cex = 0.9, title = "單一性別人數")

# (b) 預測的 alpha 最大偏誤隨人口規模，附分級門檻
plot(F1$N, 100 * F1$share_below_mustar, log = "x", type = "o", pch = 16,
     lwd = 2.4, col = "#1f4e79", xlab = "單一性別人數 N", ylim = c(0, 60),
     ylab = expression("期望死亡數 <"~mu^"*"~"的格子（%）"),
     main = "(b) 分級判準與門檻")
abline(h = c(2, 15, 35), lty = 2, col = "grey45")
text(5.3e3, 4.5, "A 可直接編表",    pos = 4, cex = 0.9)
text(1e6, 8.0,  "B 建議校正",       pos = 2, cex = 0.9)
text(1e6, 24,   "C 必須校正",       pos = 2, cex = 0.9)
text(1e6, 46,   "D 不建議單獨編表", pos = 2, cex = 0.9)
axis(4, at = c(2, 15, 35), labels = FALSE, tcl = -0.3)
dev.off()
cat("\n圖已輸出：output/figures/figF1_feasibility.png\n")
