###############################################################################
# 二階解析量的實測：對數轉換的變異數、偏誤的導數，與校正的變異數代價
#
#   受檢驗的三項陳述（皆由 R/core/analytic_variance.R 推導）：
#     V1  v(mu; c0) = Var[log max(D, c0)] 可預測標準 LC 逐年齡 alpha 的
#         抽樣變異數，Var(alpha_hat_x) = T^{-2} sum_t v(mu_{x,t})。
#     V2  b'(mu; c0) = E[log(D+1)] - E[log max(D,c0)] - 1/mu 為精確式，
#         且 mu b'(mu) 的極限為 -1（mu->0）與 0（mu->inf）。
#     V3  校正的變異數膨脹因子為 (1 - c_x)^2，c_x = mean_t[mu b'(mu)]。
#
#   設定與 run_analytic_bias.R 完全相同（女性 2001-2024、全國配適為真值、
#   N = 1e4/5e4/2e5、重複 100 次、種子 20260926 + 1000 i + r），
#   故各表可與 1005 版的表逐格對照。
#
# 輸出：output/tables/tableV1_logvar_props.csv
#       output/tables/tableV2_alpha_sd_byage.csv
#       output/tables/tableV3_inflation_byband.csv
#       output/figures/figV1_logvar.png
###############################################################################

source("R/core/lc_poisson_lasso.R")
source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R")
source("R/core/analytic_bias_lc.R")
source("R/core/analytic_variance.R")
source("R/core/fig_axis_utils.R")

SEED <- 20260926; REPS <- 100; C0 <- 0.5
NS <- c(1e4, 5e4, 2e5); YEARS <- 2001:2024
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)

dat <- load_data(sex = "Female"); keep <- which(dat$years %in% YEARS)
D0 <- dat$D[, keep, drop = FALSE]; E0 <- dat$E[, keep, drop = FALSE]
A <- nrow(D0); Tn <- ncol(D0)
truth <- lc_poisson_firth(D0, E0, firth = FALSE)
mtrue <- exp(outer(truth$a, rep(1, Tn)) + outer(truth$b, truth$k))
wage  <- E0[, Tn] / sum(E0[, Tn])
BF <- make_logbias(C0)

## ---- 表 V1：三個解析量在代表性 mu 上的值，附模擬核對 -------------------
grid <- c(0.02, 0.05, 0.1, 0.3, 0.5, 0.6931, 0.9234, 1, 2, 5, 10, 20, 50, 100, 300)
set.seed(20261002)
v_sim <- vapply(grid, function(m) var(log(pmax(rpois(2e6, m), C0))), numeric(1))
b_num <- vapply(grid, function(m) {
  h <- m * 1e-5
  (logbias_exact(m + h, C0) - logbias_exact(m - h, C0)) / (2 * h)
}, numeric(1))
V1 <- data.frame(
  mu          = grid,
  b           = logbias_exact(grid, C0),
  v_closed    = logvar_exact(grid, C0),
  v_sim       = v_sim,
  v_ratio     = logvar_exact(grid, C0) / v_sim,
  bprime      = logbias_deriv(grid, C0),
  bprime_num  = b_num,
  bprime_ratio= logbias_deriv(grid, C0) / b_num,
  mu_bprime   = logbias_elast(grid, C0),
  inflate     = (1 - logbias_elast(grid, C0))^2
)
write.csv(V1, "output/tables/tableV1_logvar_props.csv", row.names = FALSE)
cat("== 表 V1 ==\n"); print(round(V1, 4))

## ---- 表 V2：逐年齡的預測與實測 -----------------------------------------
V2 <- NULL
for (ni in seq_along(NS)) {
  N <- NS[ni]
  E  <- N * outer(wage, rep(1, Tn)); MU <- E * mtrue
  pred <- alpha_var_predict(MU, C0)
  raw <- corr <- matrix(NA_real_, A, REPS)
  for (r in seq_len(REPS)) {
    set.seed(SEED + 1000 * ni + r)
    D  <- matrix(rpois(A * Tn, MU), A, Tn)
    f0 <- lc_svd_fit(log(pmax(D, C0) / E)); raw[, r] <- f0$a
    muh <- E * exp(outer(f0$a, rep(1, Tn)) + outer(f0$b, f0$k))
    corr[, r] <- f0$a - rowMeans(matrix(BF(as.vector(muh)), A, Tn))
  }
  sd_raw <- apply(raw, 1, sd); sd_cor <- apply(corr, 1, sd)
  V2 <- rbind(V2, data.frame(
    N = N, age_group = seq_len(A), age = rownames(D0),
    mu_median = apply(MU, 1, median),
    sd_closed = pred$sd_raw, sd_sim = sd_raw,
    sd_ratio  = sd_raw / pred$sd_raw,
    c_x = pred$c_x,
    inflate_closed = pred$inflate, inflate_sim = (sd_cor / sd_raw)^2))
}
write.csv(V2, "output/tables/tableV2_alpha_sd_byage.csv", row.names = FALSE)
cat("\n== 表 V2 摘要（各規模的中位數）==\n")
print(aggregate(cbind(sd_closed, sd_sim, sd_ratio, inflate_closed, inflate_sim)
                ~ N, V2, median))

## ---- 表 V3：依期望死亡數分段彙總，顯示一階展開的有效範圍 ---------------
V2$band <- cut(V2$mu_median, breaks = c(0, 0.5, 1, 2, 10, Inf),
               labels = c("mu<0.5", "0.5-1", "1-2", "2-10", "mu>10"),
               right = FALSE)
V3 <- do.call(rbind, lapply(split(V2, list(V2$N, V2$band), drop = TRUE), function(d)
  data.frame(N = d$N[1], band = as.character(d$band[1]), n_age = nrow(d),
             sd_ratio      = median(d$sd_ratio),
             inflate_closed= median(d$inflate_closed),
             inflate_sim   = median(d$inflate_sim),
             inflate_err   = median(d$inflate_closed - d$inflate_sim))))
V3 <- V3[order(V3$N, V3$band), ]
write.csv(V3, "output/tables/tableV3_inflation_byband.csv", row.names = FALSE)
cat("\n== 表 V3 ==\n"); print(V3, row.names = FALSE)

## ---- 圖 V1：v(mu)、mu b'(mu) 與膨脹因子 --------------------------------
g  <- exp(seq(log(0.01), log(500), length.out = 400))
vg <- logvar_exact(g, C0); eg <- logbias_elast(g, C0); bg <- logbias_exact(g, C0)
png_cjk("output/figures/figV1_logvar.png", width = 2250, height = 820, res = 185)
par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.6, 0.8), mgp = c(2.6, 0.8, 0))
plot(g, vg, log = "x", type = "l", lwd = 2.4, col = "#1f4e79",
     xlab = expression(mu~"（單年期望死亡數）"), ylab = expression(v(mu*";"*c[0])),
     main = "(a) 對數尺度的變異數")
abline(v = g[which.max(vg)], lty = 3, col = "grey40")
text(0.012, max(vg) * 0.90,
     sprintf("極大於 mu = %.2f\n v = %.3f", g[which.max(vg)], max(vg)),
     pos = 4, cex = 0.95)
plot(g, eg, log = "x", type = "l", lwd = 2.4, col = "#9e2a2b",
     xlab = expression(mu), ylab = expression(mu*b*minute*(mu)),
     main = "(b) 偏誤導數的彈性", ylim = c(-1.05, 0.25))
abline(h = c(-1, 0), lty = 3, col = "grey40")
text(0.02, -0.93, "極限 -1", pos = 4, cex = 0.95)
text(120, 0.07, "極限 0", pos = 2, cex = 0.95)
plot(g, (1 - eg)^2, log = "xy", type = "l", lwd = 2.4, col = "#2d6a4f",
     xlab = expression(mu), ylab = expression((1-c[x])^2),
     main = "(c) 校正的變異數膨脹因子")
abline(h = c(1, 4), lty = 3, col = "grey40")
text(0.02, 3.4, "極限 4", pos = 4, cex = 0.95)
text(0.012, 1.03, "極限 1", pos = 4, cex = 0.95)
text(4.43, 0.775, sprintf("極小 %.3f 於 mu = %.2f", 0.7640, 4.4288),
     pos = 1, offset = 0.9, cex = 0.95)
dev.off()
cat("\n圖已輸出：output/figures/figV1_logvar.png\n")
