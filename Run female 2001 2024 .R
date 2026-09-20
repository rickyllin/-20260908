###############################################################################
# Female, 2001-2024 的分析
# 先 source 主檔與修正檔，再 source 本檔。
###############################################################################

source("rabbi_mazzuco_replication.R")
source("Patch adjust kappa.R")

## ---- 在 load_data() 外面加年份篩選 ---------------------------------------
subset_years <- function(dat, y0, y1) {
  keep <- which(dat$years >= y0 & dat$years <= y1)
  stopifnot(length(keep) > 5)
  list(D = dat$D[, keep, drop = FALSE], E = dat$E[, keep, drop = FALSE],
       ages = dat$ages, years = dat$years[keep])
}

dat <- subset_years(load_data(sex = "Female"), 2001, 2024)
D <- dat$D; E <- dat$E
cat(sprintf("Female %d-%d：%d 個年齡組 x %d 年；零死亡格 %d；最小死亡數 %.0f\n",
            min(dat$years), max(dat$years), nrow(D), ncol(D),
            sum(D == 0), min(D)))

## ---- (a) e0 與 e0dagger ---------------------------------------------------
mobs <- D / E
e0 <- apply(mobs, 2, function(m) life_table(m)$e0)
ed <- apply(mobs, 2, function(m) life_table(m)$edag)
print(round(data.frame(Year = dat$years, e0 = e0, e0_dagger = ed), 2),
      row.names = FALSE)

## ---- (b) 修勻 -------------------------------------------------------------
logm <- log(mobs)
sm   <- smooth_l1_2d(logm, lxx = 5, lxt = 1, ltt = 5, method = "rq")
mid  <- make_grid(nrow(D))$mid
ra <- function(M) mean(abs(D2_uneven(mid) %*% M))
rt <- function(M) mean(abs(M %*% t(D2_uneven(seq_len(ncol(M)) - 1))))
cat(sprintf("\n修勻 MAE=%.5f；年齡粗糙度 %.4f->%.4f；時間粗糙度 %.4f->%.4f\n",
            mean(abs(sm - logm)), ra(logm), ra(sm), rt(logm), rt(sm)))

## ---- (c) 四種變體 ---------------------------------------------------------
fits <- lapply(c("LC","LM","LCP","LCedag"), function(v) fit_variant(D, E, v))
names(fits) <- c("LC","LM","LCP","LCedag")
cat("\nkappa 性質：\n")
for (v in names(fits)) {
  k <- fits[[v]]$k; n <- length(k)
  cat(sprintf("  %-7s drift=%.4f  離線性殘差 sd=%.3f  sd(diff)=%.3f  NA=%d\n",
              v, (k[n]-k[1])/(n-1), sd(residuals(lm(k ~ seq_len(n)))),
              sd(diff(k)), sum(is.na(k))))
}

## ---- (d) 樣本外：只有 24 年，h=10 會只剩 14 年配適，建議 h=5 -------------
for (h in c(5, 10)) {
  cat(sprintf("\n=== 留最後 %d 年（配適 %d 年）===\n", h, ncol(D) - h))
  print(oos_evaluate(D, E, h = h), digits = 4, row.names = FALSE)
}

## ---- 圖 -------------------------------------------------------------------
op <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
matplot(dat$years, t(logm), type = "l", col = rgb(0,0,0,.2), lty = 1,
        xlab = "Year", ylab = "log m", main = "Observed (Female)")
matplot(dat$years, t(sm), type = "l", col = rgb(0,0,1,.35), lty = 1,
        xlab = "Year", ylab = "log m", main = "L1 smoothed")
cols <- c(LC="black", LM="blue", LCP="grey50", LCedag="red")
plot(dat$years, fits$LC$k, type = "n",
     ylim = safe_range(sapply(fits, `[[`, "k")),
     xlab = "Year", ylab = expression(kappa[t]), main = "kappa_t")
for (v in names(fits)) lines(dat$years, fits[[v]]$k, col = cols[v], lwd = 2)
legend("topright", names(cols), col = cols, lwd = 2, bty = "n", cex = .8)
plot(dat$years, e0, type = "b", pch = 16, xlab = "Year", ylab = expression(e[0]),
     main = "Female e0")
par(op)

