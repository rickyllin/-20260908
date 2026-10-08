###############################################################################
# 跨國的適用性評估（一）：可行性分級的門檻是否隨國家而變
#
#   老師的問題是本文的方法是否也適用於其他地區的資料。理論上的答案可以
#   先寫出來：b(mu;c0) 只依賴 mu 與 c0，是卜瓦松分布的性質而非任何國家的
#   性質，故閉式本身逐國成立、不需重新推導。會隨國家而變的是
#   「人口規模 N 對應到哪一組 mu」這個映射，因為
#       mu_{x,t} = N * w_x * m_{x,t}，
#   其中 w_x 為年齡結構、m_{x,t} 為死亡率水準，兩者皆逐國不同。
#
#   據此，分級規則（低於 mu* 的格數佔比，門檻 2%/15%/35%）本身是國家無關的，
#   但 N 到分級的映射不是。本程式逐國算出該映射，並找一個可預測門檻的
#   單一摘要量，使跨國使用不必逐國重跑。
#
#   資料：HMD 的澳洲、加拿大、日本、紐西蘭、美國，加臺灣，女性 1998-2021。
#
# 輸出：output/tables/tableX1_cross_country_grade.csv
###############################################################################
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R");  source("R/core/analytic_bias_lc.R")
source("R/core/lc_feasibility.R");   source("R/core/hmd_load.R")
c0 <- 0.5; MS <- mustar(c0)
cat(sprintf("變號點 mu* = %.4f（c0 = %.1f）\n\n", MS, c0))

## 逐國取真值：以全國資料的 Poisson 配適為真，取末年的死亡率與年齡結構
info <- list()
for (nm in ALL_COUNTRIES) {
  x <- get_country(nm); A <- nrow(x$D); Tn <- ncol(x$D)
  f <- lc_poisson_firth(round(x$D), x$E, firth = FALSE)
  mt <- exp(outer(f$a, rep(1, Tn)) + outer(f$b, f$k))
  w  <- x$E[, Tn] / sum(x$E[, Tn])          # 年齡結構（末年）
  info[[nm]] <- list(mt = mt, w = w, A = A, Tn = Tn,
                     cdr = sum(x$D) / sum(x$E),            # 粗死亡率
                     e0 = life_table(mt[, Tn])$e0)
}

## 給定 N，算低於 mu* 的格數佔比
share <- function(nm, N) {
  z <- info[[nm]]
  mu <- outer(N * z$w, rep(1, z$Tn)) * z$mt
  mean(mu < MS)
}
grade_of <- function(p) if (p < 0.02) "A" else if (p < 0.15) "B" else
                        if (p < 0.35) "C" else "D"
## 解出各級的 N 門檻
thresh <- function(nm, p) {
  f <- function(lN) share(nm, exp(lN)) - p
  if (f(log(1e2)) < 0) return(NA_real_)
  if (f(log(5e7)) > 0) return(Inf)
  exp(uniroot(f, c(log(1e2), log(5e7)), tol = 1e-8)$root)
}

cat(sprintf("%-10s %8s %8s %12s %12s %12s\n", "國家", "e0", "粗死亡率",
            "A 級門檻", "B 級門檻", "C 級門檻"))
rows <- list()
for (nm in ALL_COUNTRIES) {
  tA <- thresh(nm, 0.02); tB <- thresh(nm, 0.15); tC <- thresh(nm, 0.35)
  rows[[nm]] <- data.frame(國家 = nm, e0 = round(info[[nm]]$e0, 2),
    粗死亡率 = signif(info[[nm]]$cdr, 4),
    A門檻 = round(tA), B門檻 = round(tB), C門檻 = round(tC))
  cat(sprintf("%-10s %8.2f %8.5f %12.0f %12.0f %12.0f\n",
      nm, info[[nm]]$e0, info[[nm]]$cdr, tA, tB, tC))
}
cat("\n（門檻為單一性別人數；低於 C 級門檻者為 D 級。）\n")

cat("\n### 同一人數下各國的分級\n")
NS <- c(5e3, 1e4, 2e4, 5e4, 1e5, 2e5)
cat(sprintf("%-10s %s\n", "國家",
            paste(sprintf("%10s", formatC(NS, format="d", big.mark=",")), collapse="")))
for (nm in ALL_COUNTRIES) {
  g <- sapply(NS, function(N) sprintf("%s(%.0f%%)", grade_of(share(nm,N)), 100*share(nm,N)))
  cat(sprintf("%-10s %s\n", nm, paste(sprintf("%10s", g), collapse="")))
}
tab <- do.call(rbind, rows)
write.csv(tab, "output/tables/tableX1_cross_country_grade.csv", row.names = FALSE)
cat("\n已輸出 tableX1_cross_country_grade.csv\n")
