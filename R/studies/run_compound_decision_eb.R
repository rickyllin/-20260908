###############################################################################
# Poisson 複合決策的非參數經驗貝氏（研究筆記/複合決策與經驗貝氏_評估.md）
#
#   Favaro and Fortini (2024, arXiv:2411.07651) 以 Smith-Makov / Newton 的
#   預測遞迴估計混合分布 G，再以 Robbins 形式 (y+1)p(y+1)/p(y) 給出 Poisson
#   均值的 Bayes 解。遞迴為
#       G_{n+1} = (1-a_{n+1}) G_n + a_{n+1} Pois(Y|theta)G_n / \int(...)，
#   a_n = (a+n)^{-gamma}，gamma in (1/2, 1]。
#
#   本程式檢驗三件事（對應該筆記第三節）：
#     3.1 可交換性。把 528 格併成同一個 G 是誤設——mu 跨 475 倍，5-9 歲
#         真值 0.05 會被拉到 0.38（高估 7.6 倍）。正確單位是固定年齡、
#         跨 368 個鄉鎮市區。
#     3.2 曝露不等。鄉鎮市區人口差數個數量級，G 須放在率 lambda 上並把
#         E_k 放進概似；照該文把 G 放在計數均值上效果大打折扣
#         （5-9 歲 3.16 對 0.115，oracle 0.063）。遞迴照樣跑得動，但該文
#         定理 2.1 與 2.4 假設同分布，offset 的情形超出其理論。
#     3.3 對數泛函。E[log theta|y] 在理論上最適（oracle 確認），但以估出的
#         G 計算時會繼承網格下界——下界由 1e-2 降到 1e-6，均方誤差惡化 88
#         倍，而 log E[theta|y] 只動 1.8 倍。網格下界在此扮演的角色與 c0
#         在零格替代中完全相同。
#
#   順序依賴：預測遞迴依觀測次序而定，故各量皆取二十次隨機排列的平均。
#
# 輸出：終端三組表格
###############################################################################
## ---- 3.1 可交換性 --------------------------------------------------------
source("R/core/lc_poisson_lasso.R"); source("R/core/rabbi_mazzuco_replication.R")
source("R/core/fh_firth_kalman.R"); source("R/core/analytic_bias_lc.R")
c0 <- 0.5
## Newton / Smith-Makov 預測遞迴（Favaro and Fortini 2024, 式 (2)）
pr <- function(y, grid, a = 1, gam = 0.67, perms = 20, seed = 1) {
  m <- length(grid); G <- matrix(0, perms, m)
  for (p in seq_len(perms)) {
    set.seed(seed + p); ord <- sample(length(y)); g <- rep(1/m, m)
    for (i in seq_along(ord)) {
      lik <- dpois(y[ord[i]], grid); num <- lik * g
      al <- (a + i)^(-gam)
      g <- (1 - al) * g + al * num / sum(num)
    }
    G[p, ] <- g
  }
  colMeans(G)
}
post <- function(y, grid, g) {            # 回傳 E[theta|y] 與 E[log theta|y]
  t(sapply(y, function(yy) {
    w <- dpois(yy, grid) * g; w <- w / sum(w)
    c(sum(grid * w), sum(log(grid) * w))
  }))
}

cat("### 一、可交換性：把 528 格視為同一個 G 的抽樣是否合理\n")
dat <- load_data(sex="Female"); keep <- which(dat$years %in% 2001:2024)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]
A <- nrow(D0); Tn <- ncol(D0); ages <- rownames(D0)
tr <- lc_poisson_firth(D0, E0, firth=FALSE)
mt <- exp(outer(tr$a,rep(1,Tn))+outer(tr$b,tr$k)); wage <- E0[,Tn]/sum(E0[,Tn])
N <- 1e4; E <- outer(N*wage, rep(1,Tn)); MU <- E*mt
cat(sprintf("  mu 的全距：%.3f（%s）至 %.1f（%s），跨 %.0f 倍\n",
    min(MU), ages[which.min(apply(MU,1,mean))], max(MU),
    ages[which.max(apply(MU,1,mean))], max(MU)/min(MU)))
set.seed(99); D <- matrix(rpois(length(E), MU), A, dimnames=dimnames(E))
grid <- exp(seq(log(0.01), log(600), length.out = 300))
g_all <- pr(as.vector(D), grid)
pp <- post(c(0,0,1,2,50,300), grid, g_all)
cat("  以全部 528 格估一個 G，再看幾個觀測值的後驗：\n")
for (i in seq_len(6)) cat(sprintf("    y=%3d → E[θ|y]=%8.3f, E[logθ|y]=%+7.3f, logE[θ|y]=%+7.3f（差 %.3f）\n",
    c(0,0,1,2,50,300)[i], pp[i,1], pp[i,2], log(pp[i,1]), log(pp[i,1])-pp[i,2]))
tgt <- c(which(ages=="5-9"), which(ages=="85-89"))
for (x in tgt) {
  o <- post(D[x,], grid, g_all)
  cat(sprintf("  %s：真 mu 均 %.2f，E[θ|y] 均 %.2f（%s）\n", ages[x],
      mean(MU[x,]), mean(o[,1]),
      ifelse(mean(o[,1])>mean(MU[x,])*1.5,"被高齡格拉高",
      ifelse(mean(o[,1])<mean(MU[x,])*0.67,"被低齡格拉低","尚可"))))
}

cat("\n### 二、正確的可交換單位：固定年齡、跨地區\n")
K <- 368                                   # 臺灣鄉鎮市區數
for (xn in c("5-9","1-4","85-89")) {
  x <- which(ages == xn)
  mu_nat <- mean(mt[x,])                    # 全國該年齡死亡率
  set.seed(7)
  area_e <- exp(rnorm(K, 0, 0.25))          # 地區效應，對數常態
  Ek <- N * wage[x]                         # 各區該年齡曝露
  muk <- Ek * mu_nat * area_e
  y <- rpois(K, muk)
  gk <- pr(y, grid <- exp(seq(log(1e-3), log(max(10*max(muk),1)), length.out=300)))
  o  <- post(y, grid, gk)
  lt <- log(muk)
  m1 <- mean((log(pmax(y,c0)/1) - lt)^2)    # 零格替代（以 E=1 尺度比較對數偏移）
  m2 <- mean((log(o[,1]) - lt)^2)           # log E[theta|y]
  m3 <- mean((o[,2] - lt)^2)                # E[log theta|y]  ← 正確的泛函
  cat(sprintf("  %-6s mu 均 %.2f，零格 %d/%d｜MSE(log)：零格替代 %.4f｜logE[θ|y] %.4f｜E[logθ|y] %.4f\n",
      xn, mean(muk), sum(y==0), K, m1, m2, m3))
}

## ---- 3.3 網格下界與兩個泛函 ----------------------------------------------

## ---- 3.2 曝露不等 --------------------------------------------------------
pr_mu <- function(y, grid, a=1, gam=0.67, perms=20, seed=1) {
  m<-length(grid); G<-matrix(0,perms,m)
  for (p in seq_len(perms)) { set.seed(seed+p); o<-sample(length(y)); g<-rep(1/m,m)
    for (i in seq_along(o)) { num<-dpois(y[o[i]],grid)*g
      al<-(a+i)^(-gam); g<-(1-al)*g+al*num/sum(num) }
    G[p,]<-g }; colMeans(G) }
pr_rate <- function(y, Ek, grid, a=1, gam=0.67, perms=20, seed=1) {
  m<-length(grid); G<-matrix(0,perms,m)
  for (p in seq_len(perms)) { set.seed(seed+p); o<-sample(length(y)); g<-rep(1/m,m)
    for (i in seq_along(o)) { j<-o[i]; num<-dpois(y[j], Ek[j]*grid)*g
      al<-(a+i)^(-gam); g<-(1-al)*g+al*num/sum(num) }
    G[p,]<-g }; colMeans(G) }

dat<-load_data(sex="Female"); keep<-which(dat$years %in% 2001:2024)
D0<-dat$D[,keep,drop=FALSE]; E0<-dat$E[,keep,drop=FALSE]; Tn<-ncol(D0); ages<-rownames(D0)
tr<-lc_poisson_firth(D0,E0,firth=FALSE)
mt<-exp(outer(tr$a,rep(1,Tn))+outer(tr$b,tr$k)); wage<-E0[,Tn]/sum(E0[,Tn])
K<-368
cat("### 曝露不等時，G 放在計數均值上是否仍可行\n")
cat("    鄉鎮市區人口由 1,000 至 500,000（對數均勻），各區死亡率另有對數常態擾動\n\n")
set.seed(11); Npop <- exp(seq(log(1e3), log(5e5), length.out=K))[sample(K)]
cat(sprintf("  %-7s %7s %6s %11s %11s %11s %11s\n","年齡","率均","零格",
            "零格替代","G on mu","G on lambda","oracle"))
for (xn in c("1-4","5-9","15-19","40-44","85-89")) {
  x<-which(ages==xn); mu_nat<-mean(mt[x,])
  Ek <- Npop*wage[x]
  set.seed(7); lam_k <- mu_nat*exp(rnorm(K,0,0.25)); muk <- Ek*lam_k
  y <- rpois(K, muk); lt <- log(lam_k)                 # 目標為對數死亡「率」
  gr_m <- exp(seq(log(1e-3), log(10*max(muk)), length.out=400))
  gm <- pr_mu(y, gr_m)
  em <- sapply(y, function(yy){w<-dpois(yy,gr_m)*gm; w<-w/sum(w); sum(gr_m*w)})
  gr_l <- exp(seq(log(mu_nat/300), log(mu_nat*300), length.out=400))
  gl <- pr_rate(y, Ek, gr_l)
  el <- sapply(seq_len(K), function(k){w<-dpois(y[k],Ek[k]*gr_l)*gl; w<-w/sum(w); sum(gr_l*w)})
  ora <- sapply(seq_len(K), function(k){
    w<-dpois(y[k],Ek[k]*gr_l)*dlnorm(gr_l/mu_nat,0,0.25); w<-w/sum(w); sum(gr_l*w)})
  cat(sprintf("  %-7s %7.1e %3d/%d %11.4f %11.4f %11.4f %11.4f\n", xn, mean(lam_k),
      sum(y==0), K,
      mean((log(pmax(y,c0)/Ek)-lt)^2),
      mean((log(em/Ek)-lt)^2),
      mean((log(el)-lt)^2),
      mean((log(ora)-lt)^2)))
}
cat("\n  三欄皆為 log 死亡率的均方誤差。G on mu 忽略曝露差異，G on lambda 把 E_k 放進概似。\n")
