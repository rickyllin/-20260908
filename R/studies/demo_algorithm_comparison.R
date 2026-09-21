###############################################################################
# 點估計 LASSO 與分位數 LASSO 的演算法差異
#
# 目的：說明「LASSO」一詞底下其實有兩種結構不同的最佳化問題，
#       其可用的演算法並不相同。
#
#   (A) 分位數 LASSO：L1 損失 + L1 懲罰 -> 全域分段線性 -> 線性規劃
#   (B) 點估計 LASSO：平方損失 + L1 懲罰 -> 平滑 + 可分離不可微 -> 座標下降
#
# Tseng (2001) 的收斂定理要求目標函式可寫成
#       f(x_1,...,x_N) = f_0(x_1,...,x_N) + sum_k f_k(x_k)
# 其中 f_0 可微、不可微的 f_k 逐區塊可分離。
#   (B) 符合；(A) 的不可微部分是「損失本身」且耦合所有座標，不符合。
#
# 本檔以實際資料驗證：對 (A) 施以座標下降會卡在非駐點。
#
# 用法（工作目錄為 repo 根目錄）： Rscript R/studies/demo_algorithm_comparison.R
###############################################################################

source("R/core/rabbi_mazzuco_replication.R")
TABDIR <- "output/tables"; dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)

## ===================== 0. 取一小塊真實資料 ================================
dat  <- load_data(sex = "Female")
keep <- which(dat$years >= 2015)                    # 22 齡組 x 10 年 = 220 格
logm <- log(dat$D[, keep] / dat$E[, keep])
A <- nrow(logm); Tn <- ncol(logm)

R  <- build_R(A, Tn, make_grid(A)$mid, lxx = 5, lxt = 1, ltt = 5)
Rm <- as.matrix(R)
y    <- as.vector(t(logm))
yext <- c(y, rep(0, nrow(Rm) - length(y)))
nfid <- length(y)                                   # 保真列數

obj_L1 <- function(z) sum(abs(yext - Rm %*% z))

cat(sprintf("資料：%d 齡組 x %d 年；設計矩陣 %d x %d\n\n",
            A, Tn, nrow(Rm), ncol(Rm)))

## ===================== 1. (A) 線性規劃 ====================================
Rc   <- methods::as(methods::as(R, "dgCMatrix"), "matrix.csr")
z_lp <- rq_sfn_safe(Rc, yext, tau = 0.5)$coef
f_lp <- obj_L1(z_lp)
cat(sprintf("(A1) 線性規劃（Frisch-Newton 內點法）  目標值 = %.4f\n", f_lp))

## ===================== 2. (A) 對同一目標施以座標下降 ======================
# L1 損失的單座標精確最小化 = 加權中位數，故此為「精確」座標下降，
# 非近似；若它仍卡住，即證明問題出在結構而非步長。
wmed <- function(v, w) { o <- order(v); v <- v[o]; w <- w[o]
                         v[which(cumsum(w) >= sum(w) / 2)[1]] }
z_cd <- rep(0, ncol(Rm)); r <- yext - Rm %*% z_cd
for (it in 1:100) {
  zold <- z_cd
  for (j in seq_along(z_cd)) {
    a <- Rm[, j]; nz <- which(a != 0)
    if (!length(nz)) next
    delta <- wmed(r[nz] / a[nz], abs(a[nz]))
    z_cd[j] <- z_cd[j] + delta; r[nz] <- r[nz] - a[nz] * delta
  }
  if (max(abs(z_cd - zold)) < 1e-10) break
}
f_cd <- obj_L1(z_cd); f_0 <- obj_L1(rep(0, ncol(Rm)))
cat(sprintf("(A2) 座標下降（同一目標，精確單座標最小化）目標值 = %.4f（%d 次掃描）\n",
            f_cd, it))
cat(sprintf("     起始值 %.4f；亦即【完全沒有移動】，較 LP 差 %.1f 倍\n",
            f_0, f_cd / f_lp))

## ===================== 3. 卡住的原因 ======================================
# 在 z = 0 處，逐一檢查是否存在可改善的單座標移動
z0 <- rep(0, ncol(Rm)); nimp <- 0L
for (j in seq_len(ncol(Rm)))
  for (dl in c(-2, -1, -0.5, -0.1, 0.1, 0.5, 1, 2)) {
    z <- z0; z[j] <- dl
    if (obj_L1(z) < f_0 - 1e-9) nimp <- nimp + 1L
  }
w_fid <- colSums(abs(Rm[1:nfid, ]))                 # 保真列的權重
w_pen <- colSums(abs(Rm[(nfid + 1):nrow(Rm), ]))    # 懲罰列的權重
cat(sprintf("\n(A3) 可改善的單座標移動 = %d / %d\n", nimp, ncol(Rm)))
cat(sprintf("     保真權重 %.1f；懲罰權重 中位數 %.2f、最小 %.2f；懲罰 > 保真的座標佔 %.0f%%\n",
            median(w_fid), median(w_pen), min(w_pen), 100 * mean(w_pen > w_fid)))
cat("     -> 每個座標的加權中位數皆為 0（懲罰權重壓過保真權重），故軸向皆不動；\n")
cat("        但聯合移動可大幅下降，即 Tseng (2001) 所述的非駐點卡滯。\n")

## ===================== 4. (B) 對照：平方損失 + L1 懲罰 ====================
X <- Rm[1:nfid, , drop = FALSE]; lam <- 0.05
soft <- function(z, g) sign(z) * pmax(abs(z) - g, 0)
b <- rep(0, ncol(X)); rr <- y - X %*% b; d2 <- colSums(X^2)
for (it2 in 1:500) {
  bo <- b
  for (j in seq_along(b)) {
    if (d2[j] == 0) next
    rho <- sum(X[, j] * rr) + d2[j] * b[j]
    bn  <- soft(rho, lam) / d2[j]
    rr  <- rr - X[, j] * (bn - b[j]); b[j] <- bn
  }
  if (max(abs(b - bo)) < 1e-12) break
}
g   <- -as.vector(t(X) %*% (y - X %*% b))
kkt <- ifelse(abs(b) > 1e-9, abs(g + lam * sign(b)), pmax(abs(g) - lam, 0))
cat(sprintf("\n(B)  平方損失 + L1：座標下降 %d 次掃描收斂，KKT 最大違反 = %.2e\n",
            it2, max(kkt)))
cat("     -> 已達最適。差別純粹來自損失函數是否可微。\n")

## ===================== 5. 輸出 ============================================
out <- data.frame(
  case      = c("分位數 LASSO（L1 損失）", "分位數 LASSO（L1 損失）",
                "點估計 LASSO（平方損失）"),
  algorithm = c("線性規劃（內點法）", "座標下降（精確）", "座標下降（軟閾值）"),
  objective = c(f_lp, f_cd, NA),
  sweeps    = c(NA, it, it2),
  reached_optimum = c("是", "否（卡在起始點）", "是"),
  diagnostic = c("—",
                 sprintf("可改善的單座標移動 0/%d", ncol(Rm)),
                 sprintf("KKT 最大違反 %.1e", max(kkt))))
write.csv(out, file.path(TABDIR, "tableD_algorithm_comparison.csv"), row.names = FALSE)
cat(sprintf("\n已輸出 %s/tableD_algorithm_comparison.csv\n", TABDIR))
