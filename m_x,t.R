###############################################################################
# Lee-Carter model on Taiwan 5-year age-group data (1970-2024)
#   Death   : 五齡Death.xlsx      (Year, Age, Female, Male, Total)
#   Exposure: 5-age_exposure.xls  (Year, Age, Female, Male, Total)
#
# 注意：exposure 檔的 Age 欄被 Excel 誤判成日期（1-4 / 5-9 / 10-14 變成日期），
#       本程式以「位置」重新指派年齡標籤，不依賴檔案內的 Age 文字。
###############################################################################

library(readxl)

## ---------------------------------------------------------------- 0. 設定 --
setwd("/Users/linzili/政大課堂/碩一/碩士論文/碩士論文＿20260908/data")
path_death <- "五齡Death.xlsx"
path_expo  <- "5-age exposure.xls"

sex        <- "Total"   # "Total" / "Male" / "Female"
open_age   <- "100+"    # 併尾組：110+ 在 1970-2009 暴露數為 0，必須合併
#   可選 "95+", "100+", "none"

## ------------------------------------------------- 1. 讀檔並修正年齡標籤 --
age_labels <- c("0", "1-4", "5-9", "10-14",
                paste0(seq(15, 105, by = 5), "-", seq(19, 109, by = 5)),
                "110+")                                   # 共 24 組

read_mort <- function(path, labs) {
  df <- as.data.frame(read_excel(path, col_types = c("numeric", "text",
                                                     "numeric", "numeric",
                                                     "numeric")))
  names(df)[1:2] <- c("Year", "Age")
  df <- df[order(df$Year), ]                # 每年 24 列，順序固定
  stopifnot(all(table(df$Year) == length(labs)))
  df$Age <- rep(labs, times = length(unique(df$Year)))
  df
}

death <- read_mort(path_death, age_labels)
expo  <- read_mort(path_expo,  age_labels)
death <- death[which(death$Year>2000),]
expo <- expo[which(expo$Year>2000),]
## ------------------------------------------------------ 2. 轉成矩陣 (x,t) --
to_matrix <- function(df, value, labs) {
  m <- tapply(df[[value]], list(factor(df$Age, levels = labs), df$Year), sum)
  m[labs, , drop = FALSE]
}

Dx <- to_matrix(death, sex, age_labels)
Ex <- to_matrix(expo,  sex, age_labels)
stopifnot(identical(dimnames(Dx), dimnames(Ex)))

## --------------------------------------------------------- 3. 合併尾年齡 --
collapse_top <- function(M, cut) {
  if (cut == "none") return(M)
  top <- switch(cut,
                "95+"  = c("95-99", "100-104", "105-109", "110+"),
                "100+" = c("100-104", "105-109", "110+"))
  keep <- setdiff(rownames(M), top)
  rbind(M[keep, , drop = FALSE],
        matrix(colSums(M[top, , drop = FALSE]), nrow = 1,
               dimnames = list(cut, colnames(M))))
}

Dx <- collapse_top(Dx, open_age)
Ex <- collapse_top(Ex, open_age)

if (any(Ex <= 0)) stop("仍有暴露數為 0 的儲存格，請再往下併年齡組。")
cat(sprintf("維度: %d 個年齡組 x %d 年；最小暴露數 = %.1f\n",
            nrow(Ex), ncol(Ex), min(Ex)))

ages  <- rownames(Ex)
years <- as.numeric(colnames(Ex))

## ------------------------------------------ 4. Lee-Carter：SVD 一階段估計 --
# log m_xt = a_x + b_x k_t + e_xt,  限制 sum(b) = 1, sum(k) = 0
mxt   <- Dx / Ex
logm  <- log(mxt)
if (any(!is.finite(logm)))
  warning("有 death = 0 的格子，log(m) 為 -Inf；考慮再併組或改用 Poisson ML。")

ax <- rowMeans(logm)
Z  <- logm - ax                     # 中心化後每列和為 0 => sum(k)=0 自動成立

sv <- svd(Z)
u1 <- sv$u[, 1]; v1 <- sv$v[, 1]; d1 <- sv$d[1]
if (sum(u1) < 0) { u1 <- -u1; v1 <- -v1 }   # 固定符號

bx <- u1 / sum(u1)
kt <- d1 * v1 * sum(u1)
names(ax) <- names(bx) <- ages
names(kt) <- years

cat(sprintf("第一奇異值解釋變異比例 = %.4f\n", sv$d[1]^2 / sum(sv$d^2)))
cat(sprintf("sum(b) = %.6f, sum(k) = %.2e\n", sum(bx), sum(kt)))

## --------------------------------- 5. 第二階段：調整 k_t 使總死亡數吻合 --
# 解 sum_x E_xt * exp(a_x + b_x k_t) = D_.t
adjust_kt <- function(ax, bx, Ex, Dx, k0) {
  sapply(seq_along(k0), function(j) {
    Dt <- sum(Dx[, j])
    f  <- function(k) sum(Ex[, j] * exp(ax + bx * k)) - Dt
    uniroot(f, interval = c(k0[j] - 50, k0[j] + 50), extendInt = "yes",
            tol = 1e-10)$root
  })
}

kt_adj <- adjust_kt(ax, bx, Ex, Dx, kt)
names(kt_adj) <- years

## --------------------------------------------------------- 6. 配適度診斷 --
fit_logm <- outer(ax, rep(1, length(kt_adj))) + outer(bx, kt_adj)
Dhat     <- Ex * exp(fit_logm)

dev_res  <- sign(Dx - Dhat) * sqrt(2 * (ifelse(Dx > 0, Dx * log(Dx / Dhat), 0)
                                        - (Dx - Dhat)))
cat(sprintf("Poisson deviance = %.1f (df = %d)\n",
            sum(dev_res^2, na.rm = TRUE),
            length(Dx) - (length(ax) + length(bx) - 1 + length(kt_adj) - 1)))
cat(sprintf("log(m) 的 R^2 = %.4f\n",
            1 - sum((logm - fit_logm)^2) / sum((logm - ax)^2)))

## ------------------------------------------------------------- 7. 畫圖 ----
op <- par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
plot(seq_along(ax), ax, type = "b", pch = 16, xaxt = "n",
     xlab = "Age group", ylab = expression(alpha[x]), main = "alpha_x")
axis(1, at = seq_along(ax), labels = ages, las = 2, cex.axis = .6)

plot(seq_along(bx), bx, type = "b", pch = 16, xaxt = "n",
     xlab = "Age group", ylab = expression(beta[x]), main = "beta_x")
axis(1, at = seq_along(bx), labels = ages, las = 2, cex.axis = .6)
abline(h = 0, lty = 3)

plot(years, kt, type = "l", lty = 2, col = "grey40",
     xlab = "Year", ylab = expression(kappa[t]), main = "kappa_t")
lines(years, kt_adj, lwd = 2)
legend("topright", c("SVD", "adjusted"), lty = c(2, 1),
       col = c("grey40", "black"), lwd = c(1, 2), bty = "n", cex = .8)

matplot(years, t(dev_res), type = "l", col = rgb(0, 0, 0, .3), lty = 1,
        xlab = "Year", ylab = "Deviance residual", main = "Residuals by age")
abline(h = 0, col = 2)
par(op)

## ------------------------------------------ 8. 預測：random walk with drift --
h     <- 20
n     <- length(kt_adj)
drift <- (kt_adj[n] - kt_adj[1]) / (n - 1)
sigma <- sd(diff(kt_adj) - drift)

kt_fc <- kt_adj[n] + drift * (1:h)
se    <- sigma * sqrt(1:h)
kt_lo <- kt_fc - 1.96 * se
kt_hi <- kt_fc + 1.96 * se
yr_fc <- max(years) + (1:h)

cat(sprintf("drift = %.4f / year, sigma = %.4f\n", drift, sigma))

mxt_fc <- exp(outer(ax, rep(1, h)) + outer(bx, kt_fc))
dimnames(mxt_fc) <- list(ages, yr_fc)

plot(c(years, yr_fc), c(kt_adj, kt_fc), type = "n",
     xlab = "Year", ylab = expression(kappa[t]), main = "kappa_t forecast")
polygon(c(yr_fc, rev(yr_fc)), c(kt_lo, rev(kt_hi)),
        col = rgb(0, 0, 1, .15), border = NA)
lines(years, kt_adj, lwd = 2); lines(yr_fc, kt_fc, col = 4, lwd = 2)

## ------------------------------------------------------------- 9. 輸出 ----
params <- data.frame(Age = ages, alpha = ax, beta = bx, row.names = NULL)
kappa  <- data.frame(Year = years, kappa_svd = kt, kappa_adj = kt_adj,
                     row.names = NULL)

write.csv(params, "lc_params_ax_bx.csv", row.names = FALSE)
write.csv(kappa,  "lc_kappa_t.csv",      row.names = FALSE)
write.csv(data.frame(Year = yr_fc, kappa = kt_fc, lo = kt_lo, hi = kt_hi),
          "lc_kappa_forecast.csv", row.names = FALSE)

print(head(params, 24))

###############################################################################
# 附：用 StMoMo 交叉驗證（Poisson ML 而非 SVD/最小平方）
#
# library(StMoMo); library(demography)
# dat <- structure(list(year = years, age = 0:(nrow(Dx)-1),
#                       rate = list(total = Dx/Ex), pop = list(total = Ex),
#                       type = "mortality", label = "TWN"), class = "demogdata")
# fit <- fit(lc(link = "log"), data = StMoMoData(dat, series = "total"))
# 兩者的 a_x 幾乎相同，b_x 在暴露數小的高齡組會有明顯差異——
# 這正是你論文要量化的地方。
###############################################################################

