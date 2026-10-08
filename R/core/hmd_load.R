###############################################################################
# HMD 五國資料的載入，對齊到本文的 22 個年齡組與共同年窗
#
#   HMD 的 Deaths_5x1 與 Exposures_5x1 為五歲組、逐年，共 24 組
#   （0、1-4、5-9、…、105-109、110+）。本文用 22 組（…、95-99、100+），
#   故把 100-104、105-109、110+ 併為 100+。
#
#   兩項須注意的事。其一，HMD 的死亡數因 Lexis 三角分割而為非整數
#   （例如澳洲 1921 年 0 歲為 3842.31），故不可直接視為卜瓦松計數；
#   本文的用法是以各國全國資料配適得到真值的死亡率，再據以模擬整數計數，
#   如此 HMD 的非整數性不影響推論。其二，各國涵蓋的年度不同
#   （澳 1921-2021、加 1921-2023、日 1947-2024、紐 1948-2021、美 1933-2024），
#   故取 1998-2021 這 24 年為共同年窗，使 T 與臺灣的設定一致。
###############################################################################
HMD_DIR <- c(Australia="data/Australia", Canada="data/Canada",
             Japan="data/Japan", NZ="data/NZ", US="data/US")
HMD_YEARS <- 1998:2021

#' 讀一個 HMD 的 5x1 檔，回傳 年齡 x 年度 的矩陣（指定性別）
hmd_read <- function(path, sex = "Female", years = HMD_YEARS) {
  L <- readLines(path, warn = FALSE)
  L <- L[-(1:2)]                                  # 前兩列為標題與欄名
  hdr <- strsplit(trimws(L[1]), "\\s+")[[1]]
  col <- match(sex, hdr)
  stopifnot(!is.na(col))
  d <- do.call(rbind, strsplit(trimws(L[-1]), "\\s+"))
  d <- d[nchar(d[,1]) > 0, , drop = FALSE]
  yr <- as.integer(d[,1]); ag <- d[,2]
  v  <- suppressWarnings(as.numeric(d[,col]))     # "." 表遺漏
  keep <- yr %in% years
  yr <- yr[keep]; ag <- ag[keep]; v <- v[keep]
  ## 併高齡：100-104、105-109、110+ -> 100+
  ag[ag %in% c("100-104","105-109","110+")] <- "100+"
  AG <- c("0","1-4", paste0(seq(5,95,5),"-",seq(9,99,5)), "100+")
  M <- matrix(0, length(AG), length(years), dimnames = list(AG, as.character(years)))
  idx <- cbind(match(ag, AG), match(as.character(yr), colnames(M)))
  ok <- !is.na(idx[,1]) & !is.na(idx[,2]) & !is.na(v)
  for (i in which(ok)) M[idx[i,1], idx[i,2]] <- M[idx[i,1], idx[i,2]] + v[i]
  M
}

#' 載入一國的死亡數與曝露數
hmd_country <- function(name, sex = "Female", years = HMD_YEARS) {
  dir <- HMD_DIR[[name]]
  list(D = hmd_read(file.path(dir, "Deaths_5x1.txt"), sex, years),
       E = hmd_read(file.path(dir, "Exposures_5x1.txt"), sex, years),
       name = name)
}

#' 臺灣（沿用本文既有的載入器），對齊到同一年窗
tw_country <- function(sex = "Female", years = HMD_YEARS) {
  dat <- load_data(sex = sex); k <- which(dat$years %in% years)
  list(D = dat$D[, k, drop = FALSE], E = dat$E[, k, drop = FALSE], name = "Taiwan")
}

ALL_COUNTRIES <- c("Taiwan", names(HMD_DIR))
get_country <- function(nm, sex = "Female", years = HMD_YEARS)
  if (nm == "Taiwan") tw_country(sex, years) else hmd_country(nm, sex, years)
