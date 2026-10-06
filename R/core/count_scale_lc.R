###############################################################################
# 計數尺度上的 Lee-Carter 估計：Anscombe 非線性最小平方與有界影響擬概似
#
#   兩者都不對資料取對數，故依報告命題 1，其估計方程於真值處無偏。
#   實作共用一個解法器：對 alpha、beta、kappa 三個區塊輪流做 Fisher 評分，
#   每個區塊的估計方程皆為 sum u(theta) * d theta / d psi = 0 的形式。
#
#   Anscombe（anscombe1948）。h(D) = 2 sqrt(D + 3/8)，其變異在大 mu 時約為 1。
#   平均函數有兩種取法：
#     近似 g(mu) = 2 sqrt(mu + 3/8)，即慣用的變異安定式；
#     精確 g(mu) = E[h(D)]，以截斷級數算出。
#   權重亦有兩種：
#     單位權重（教科書作法，假設變異已安定）；
#     精確變異 v_A(mu) = Var[h(D)]。
#   須注意變異安定在小 mu 處失效：mu = 0.05 時 Var[h(D)] 約 0.057 而非 1，
#   蓋因該處 D = 0 幾乎必然、h(D) 幾乎為常數。以單位權重配適因而在
#   最需要資訊的年齡低估了資訊量。
#
#   有界影響擬概似（cantoni2001）。以 Huber 函數截斷 Pearson 殘差
#   r = (D - mu)/sqrt(mu)，並扣除 a(mu) = E[psi_c(r)] 以保 Fisher 一致性；
#   該修正項同樣以截斷級數精確算出，不以模擬近似。
###############################################################################

#' 截斷精確期望 E[f(D)]，D ~ Poisson(mu)
pois_E <- function(mu, f) sapply(mu, function(m) {
  K <- max(60, ceiling(m + 12 * sqrt(m) + 12)); d <- 0:K
  p <- dpois(d, m); sum(p * f(d)) / sum(p)
})

h_ans  <- function(d) 2 * sqrt(d + 3/8)
g_ans_approx <- function(mu) 2 * sqrt(mu + 3/8)
g_ans_exact  <- function(mu) pois_E(mu, h_ans)
v_ans_exact  <- function(mu) pois_E(mu, function(d) h_ans(d)^2) - g_ans_exact(mu)^2

#' Huber 的 Fisher 一致性修正 a(mu) = E[psi_c((D-mu)/sqrt(mu))]
huber_a <- function(mu, c = 1.345) sapply(mu, function(m) {
  K <- max(60, ceiling(m + 12 * sqrt(m) + 12)); d <- 0:K
  p <- dpois(d, m); r <- (d - m) / sqrt(m)
  sum(p * pmax(-c, pmin(c, r))) / sum(p)
})

#' 以插值加速：在對數網格上預先算好再內插
make_interp <- function(fun, lo = 1e-4, hi = 1e5, n = 800) {
  g <- exp(seq(log(lo), log(hi), length.out = n)); y <- fun(g)
  function(mu) approx(g, y, xout = pmin(pmax(mu, lo), hi), rule = 2)$y
}

#' 通用的計數尺度 Lee-Carter 解法器
#'
#' @param ufun 函數 (D, mu) -> 逐格的估計方程貢獻 u，其和對 theta 的導數
#'   以數值差分取得。估計方程為 sum u * d theta / d psi = 0。
#' @param mu_floor mu 的下界。整列全零時任何無偏的計數尺度估計方程都會把
#'   alpha 推向負無窮（報告第貳節），故須設下界並回報其是否生效。
lc_count_fit <- function(D, E, ufun, maxit = 200, tol = 1e-9,
                         mu_floor = 1e-6, a_lo = -30, a_hi = 5) {
  A <- nrow(D); Tn <- ncol(D)
  f0 <- lc_svd_fit(log(pmax(D, 0.5) / E))
  a <- f0$a; b <- f0$b; k <- f0$k
  mu_of <- function(a, b, k) pmax(E * exp(outer(a, rep(1, Tn)) + outer(b, k)), mu_floor)
  step <- function(par, grad_dir, blk) {
    ## 單一 Fisher 評分步：以數值差分取 d(sum u)/d par
    eps <- 1e-5
    U  <- function(p) {
      if (blk == "a") mu <- mu_of(p, b, k)
      else if (blk == "b") mu <- mu_of(a, p, k)
      else mu <- mu_of(a, b, p)
      u <- ufun(D, mu)
      if (blk == "a") rowSums(u)
      else if (blk == "b") rowSums(u * outer(rep(1, A), k))
      else colSums(u * outer(b, rep(1, Tn)))
    }
    u0 <- U(par); u1 <- U(par + eps)
    dd <- (u1 - u0) / eps
    dd[abs(dd) < 1e-12] <- sign(dd[abs(dd) < 1e-12] + 1e-30) * 1e-12
    par - u0 / dd
  }
  hit <- FALSE
  for (it in seq_len(maxit)) {
    a_new <- pmin(pmax(step(a, NULL, "a"), a_lo), a_hi)
    if (any(a_new <= a_lo + 1e-8)) hit <- TRUE
    a <- 0.5 * a + 0.5 * a_new
    b <- 0.5 * b + 0.5 * step(b, NULL, "b")
    k <- 0.5 * k + 0.5 * step(k, NULL, "k")
    k <- k - mean(k); sb <- sum(b)
    if (abs(sb) < 1e-8) sb <- 1e-8
    b <- b / sb; k <- k * sb
    if (it > 3) {
      dif <- max(abs(a - a_prev), abs(b - b_prev), abs(k - k_prev) / max(1, max(abs(k))))
      if (is.finite(dif) && dif < tol) break
    }
    a_prev <- a; b_prev <- b; k_prev <- k
  }
  names(a) <- rownames(D); names(b) <- rownames(D); names(k) <- colnames(D)
  list(a = a, b = b, k = k, iter = it, floor_hit = hit)
}

## --- 各方法的 u(D, mu) -------------------------------------------------------
u_anscombe <- function(exact_mean = TRUE, unit_weight = TRUE) {
  gE <- if (exact_mean) make_interp(g_ans_exact) else g_ans_approx
  vE <- make_interp(v_ans_exact)
  function(D, mu) {
    g <- gE(mu); dg <- mu / sqrt(mu + 3/8)        # d g / d theta
    w <- if (unit_weight) 1 else 1 / pmax(vE(mu), 1e-8)
    (h_ans(D) - g) * dg * w
  }
}
u_poisson <- function() function(D, mu) D - mu
u_huber <- function(c = 1.345) {
  aE <- make_interp(function(m) huber_a(m, c))
  function(D, mu) {
    r <- (D - mu) / sqrt(mu)
    (pmax(-c, pmin(c, r)) - aE(mu)) * sqrt(mu)
  }
}
