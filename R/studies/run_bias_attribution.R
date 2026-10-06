###############################################################################
# 偏誤落在哪個參數、哪些年齡、哪些年度（研究筆記/偏誤的歸屬_alpha_beta_kappa.md）
#
#   陷阱：alpha, beta, kappa 不加約束不可分別識別，故「beta 的偏誤」隨正規化
#   而變。此處一律在可識別的量上計算，並把三者換算到同一尺度，即對數死亡率
#   曲面的分解
#       log mhat - log m = (ahat-a) + (bhat-b)k + b(khat-k) + (bhat-b)(khat-k)。
#   四項並非正交，故佔比可超過一百（相消時），因而並列 RMS 與佔比。
#
#   三段計算：
#     一、曲面偏誤的四項分解，以及 e0 的歸因（一次只放一個參數的估計值）。
#         關鍵發現：標準 LC 的曲面偏誤由 alpha 主導（98%/85%）而其 e0 偏誤
#         由 kappa 主導（五萬時 -1.476 對合計 -1.481）——不矛盾，因為 alpha
#         的偏誤集中在對 e0 影響小的幼年組。又五萬時校正版的 e0 偏誤 +0.073
#         是 alpha(-0.326)、beta(+1.492)、kappa(-0.603) 三項相消的結果，
#         該相消沒有理由在其他泛函上成立。
#     二、kappa 的全距壓縮須拆成「逐次衰減」與「形狀不一致」兩件事。逐次全距
#         衰減到 53-57%，而中位路徑只有 16-18%，後者小得多是因為各次形狀彼此
#         不一致（與真值相關僅 0.34-0.39）。只引中位路徑會把衰減說得過重。
#     三、校正對 kappa 的改善原先未被記錄：五萬時 cor(khat,k) 由 0.393 升到
#         0.941、逐次全距由 5.602 升到 7.620。報告所記「校正後時間路徑更不
#         平滑、擺幅更大」在該規模上其實是衰減更少。
#
# 輸出：終端三組表格
###############################################################################
## ---- 一、曲面分解與逐參數的最大偏誤 -------------------------------------
SEED <- 20261006; REPS <- 100; c0 <- 0.5
dat <- load_data(sex="Female"); keep <- which(dat$years %in% 2001:2024)
D0 <- dat$D[,keep,drop=FALSE]; E0 <- dat$E[,keep,drop=FALSE]
A <- nrow(D0); Tn <- ncol(D0); ages <- rownames(D0); yrs <- dat$years[keep]
tr <- lc_poisson_firth(D0,E0,firth=FALSE)
a0 <- tr$a; b0 <- tr$b; k0 <- tr$k
mt <- exp(outer(a0,rep(1,Tn))+outer(b0,k0)); wage <- E0[,Tn]/sum(E0[,Tn])
BF <- make_logbias(c0)
for (N in c(1e4, 5e4)) {
  E <- outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu <- E*mt
  AA <- array(NA_real_,c(REPS,A,2)); BB <- array(NA_real_,c(REPS,A,2))
  KK <- array(NA_real_,c(REPS,Tn,2))
  for (r in seq_len(REPS)) {
    set.seed(SEED+17000+r)
    D <- matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    f1 <- lc_svd_fit(log(pmax(D,c0)/E))
    f2 <- lc_alpha_only(D,E,bfun=BF,weight=TRUE)
    for (j in 1:2) { f <- list(f1,f2)[[j]]
      AA[r,,j]<-f$a-a0; BB[r,,j]<-f$b-b0; KK[r,,j]<-f$k-k0 }
  }
  cat(sprintf("\n================ N = %.0e ================\n", N))
  for (j in 1:2) {
    nm <- c("標準 Lee-Carter","中心化＋加權")[j]
    da <- apply(AA[,,j],2,median); db <- apply(BB[,,j],2,median); dk <- apply(KK[,,j],2,median)
    ## 把三者的貢獻換算到同一個尺度：對數死亡率曲面
    Ta <- outer(da, rep(1,Tn))                 # alpha 項
    Tb <- outer(db, k0)                        # (beta_hat-beta) * kappa
    Tk <- outer(b0, dk)                        # beta * (kappa_hat-kappa)
    Ti <- outer(db, dk)                        # 交互項
    tot <- Ta+Tb+Tk+Ti
    rms <- function(M) sqrt(mean(M^2))
    cat(sprintf("\n-- %s --\n  曲面偏誤的 RMS 分解（對數死亡率單位）\n", nm))
    cat(sprintf("    alpha 項 %.4f｜beta 項 %.4f｜kappa 項 %.4f｜交互 %.4f｜合計 %.4f\n",
        rms(Ta), rms(Tb), rms(Tk), rms(Ti), rms(tot)))
    cat(sprintf("    佔合計平方和：alpha %.0f%%｜beta %.0f%%｜kappa %.0f%%｜交互 %.0f%%\n",
        100*sum(Ta^2)/sum(tot^2),100*sum(Tb^2)/sum(tot^2),
        100*sum(Tk^2)/sum(tot^2),100*sum(Ti^2)/sum(tot^2)))
    cat(sprintf("  kappa：全距 %.3f 對真值 %.3f（%.0f%%）｜漂移偏誤 %+.4f｜dk 與 -kappa 的相關 %+.3f\n",
        diff(range(dk+k0)), diff(range(k0)), 100*diff(range(dk+k0))/diff(range(k0)),
        (dk[Tn]-dk[1])/(Tn-1), cor(dk, -k0)))
    o <- order(-abs(da))[1:4]
    cat(sprintf("  alpha 偏誤最大的四個年齡：%s\n",
        paste(sprintf("%s(%+.3f)",ages[o],da[o]),collapse=", ")))
    o <- order(-abs(db))[1:4]
    cat(sprintf("  beta  偏誤最大的四個年齡：%s\n",
        paste(sprintf("%s(%+.4f)",ages[o],db[o]),collapse=", ")))
    o <- order(-abs(dk))[1:4]
    cat(sprintf("  kappa 偏誤最大的四個年度：%s\n",
        paste(sprintf("%d(%+.3f)",yrs[o],dk[o]),collapse=", ")))
  }
}

## ---- 二、kappa 的衰減與形狀 ---------------------------------------------
SEED<-20261006; REPS<-100; c0<-0.5
dat<-load_data(sex="Female"); keep<-which(dat$years %in% 2001:2024)
D0<-dat$D[,keep,drop=FALSE]; E0<-dat$E[,keep,drop=FALSE]
A<-nrow(D0); Tn<-ncol(D0); yrs<-dat$years[keep]
tr<-lc_poisson_firth(D0,E0,firth=FALSE); a0<-tr$a;b0<-tr$b;k0<-tr$k
mt<-exp(outer(a0,rep(1,Tn))+outer(b0,k0)); wage<-E0[,Tn]/sum(E0[,Tn]); BF<-make_logbias(c0)
cat("### kappa 的全距壓縮：逐次重複的全距，對中位路徑的全距\n")
cat("    若逐次全距接近真值而中位路徑被壓扁，則是形狀不一致而非衰減。\n\n")
cat(sprintf("  %-10s %-16s %12s %12s %12s\n","N","估計量","逐次全距中位","中位路徑全距","真值全距"))
for (N in c(1e4,5e4)) {
  E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mt
  KK<-array(NA_real_,c(REPS,Tn,2))
  for (r in seq_len(REPS)) { set.seed(SEED+17000+r)
    D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    KK[r,,1]<-lc_svd_fit(log(pmax(D,c0)/E))$k
    KK[r,,2]<-lc_alpha_only(D,E,bfun=BF,weight=TRUE)$k }
  for (j in 1:2) {
    per <- apply(KK[,,j],1,function(v) diff(range(v)))
    med <- apply(KK[,,j],2,median)
    cat(sprintf("  %-10.0e %-16s %12.3f %12.3f %12.3f\n", N,
        c("標準 LC","中心化＋加權")[j], median(per), diff(range(med)), diff(range(k0))))
  }
  ## 形狀一致性：各次 kappa 與真值 kappa 的相關
  for (j in 1:2) {
    cs <- apply(KK[,,j],1,function(v) cor(v,k0))
    cat(sprintf("      %s：cor(kappa_hat, kappa) 的中位 %.3f（四分位 %.3f, %.3f）\n",
        c("標準 LC","中心化＋加權")[j], median(cs), quantile(cs,.25), quantile(cs,.75)))
  }
}

## ---- 三、e0 的歸因 -------------------------------------------------------
SEED<-20261006; REPS<-100; c0<-0.5
dat<-load_data(sex="Female"); keep<-which(dat$years %in% 2001:2024)
D0<-dat$D[,keep,drop=FALSE]; E0<-dat$E[,keep,drop=FALSE]
A<-nrow(D0); Tn<-ncol(D0)
tr<-lc_poisson_firth(D0,E0,firth=FALSE); a0<-tr$a;b0<-tr$b;k0<-tr$k
mt<-exp(outer(a0,rep(1,Tn))+outer(b0,k0)); wage<-E0[,Tn]/sum(E0[,Tn]); BF<-make_logbias(c0)
e0t <- life_table(mt[,Tn])$e0
cat(sprintf("### e0 的歸因（真值 %.3f 歲）：一次只放入一個參數的偏誤\n\n", e0t))
cat(sprintf("  %-10s %-16s %9s %9s %9s %9s\n","N","估計量","只放 alpha","只放 beta","只放 kappa","全放"))
for (N in c(1e4,5e4)) {
  E<-outer(N*wage,rep(1,Tn)); dimnames(E)<-dimnames(D0); mu<-E*mt
  AA<-array(NA_real_,c(REPS,A,2)); BB<-array(NA_real_,c(REPS,A,2)); KK<-array(NA_real_,c(REPS,Tn,2))
  for (r in seq_len(REPS)) { set.seed(SEED+17000+r)
    D<-matrix(rpois(length(E),mu),A,dimnames=dimnames(E))
    f1<-lc_svd_fit(log(pmax(D,c0)/E)); f2<-lc_alpha_only(D,E,bfun=BF,weight=TRUE)
    for (j in 1:2){f<-list(f1,f2)[[j]]; AA[r,,j]<-f$a; BB[r,,j]<-f$b; KK[r,,j]<-f$k} }
  for (j in 1:2) {
    ah<-apply(AA[,,j],2,median); bh<-apply(BB[,,j],2,median); kh<-apply(KK[,,j],2,median)
    e <- function(a,b,k) life_table(exp(a + b*k[Tn]))$e0 - e0t
    cat(sprintf("  %-10.0e %-16s %+9.3f %+9.3f %+9.3f %+9.3f\n", N,
        c("標準 LC","中心化＋加權")[j],
        e(ah,b0,k0), e(a0,bh,k0), e(a0,b0,kh), e(ah,bh,kh)))
  }
}
cat("\n  各欄為只把該參數換成估計值（其餘用真值）所得的 e0 偏誤，單位為歲。\n")
cat("  三者相加不等於「全放」，蓋因 e0 對參數非線性且三項有交互。\n")
