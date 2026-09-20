# 小人口 Lee-Carter 模型之估計與區間推論

碩士論文的分析程式碼與報告。指導教授：余清祥。

以臺灣五齡組死亡資料（1970–2024）研究 **Lee-Carter（LC）模型在小樣本下的
估計問題**，主軸為三個依序遞進的研究問題：

1. **復刻**　Rabbi & Mazzuco (2021) 的「修勻＋壽命差異調整」在國家級資料上
   成立的優勢，改換資料與配適期後是否仍成立？
2. **診斷**　小樣本下 LC 失效的機制是**變異數增大**還是**有方向的偏誤**？
3. **改善**　把懲罰項直接加在年齡參數 $\beta_x$ 上是否可行？
   LASSO／Ridge／Elastic Net 如何取捨？

---

## 目錄結構

```
R/
  core/      核心工具（被其他腳本 source）
    lc_poisson_lasso.R          load_data()、卜瓦松 LC、L1 懲罰
    rabbi_mazzuco_replication.R 差分算子、2D L1 修勻、生命表、LC 變體
    patch_adjust_kappa.R        修正 adjust_kappa（e0† 對 κ 為單峰）
    penalized_lc.R              Elastic Net 族：LASSO/Ridge/EN/MCP/鬆弛/自適應
    fig_axis_utils.R            座標軸格式（1e-03 → 10^-3）
  studies/   實驗腳本（各自可獨立 Rscript 執行）
    論文復刻.R                   復刻 R&M 全部圖表
    run_female_2001_2024.R      女性 2001–2024 單次分析
    mc_exposure.R               暴露數 × 估計方法的蒙地卡羅
    lambda_selection.R          λ 的四種選法 vs oracle
    lambda3_and_resampling.R    三方向 λ 格點、重抽樣
    run_unified_mc.R            統一模擬規模（reps=100、固定種子）
    run_penalized_study.R       懲罰式 LC：收斂、偏誤—變異數、變異數速率
    run_appendix_TN.R           附錄的逐年齡比較、T×N 二因子、EN 變體
    make_figures.R              由 output/tables/ 重繪出版用圖
  legacy/
    lee_carter_basic.R          早期的獨立腳本，保留作對照
data/        原始資料
output/
  figures/   所有 .png
  tables/    所有 .csv
報告/         LaTeX 報告（詳細版／簡版／簡報），見該目錄的 README
```

---

## 快速開始

```r
install.packages(c("readxl", "quantreg", "SparseM", "Matrix"))
```

所有腳本都以 **repo 根目錄**為工作目錄（用 `碩士論文＿20260908.Rproj` 開啟
RStudio 即可），彼此以 `source("R/core/xxx.R")` 相對引用。

```sh
Rscript R/studies/論文復刻.R          # 復刻全部圖表
Rscript R/studies/run_unified_mc.R    # 統一規模的蒙地卡羅（reps=100）
Rscript R/studies/run_penalized_study.R
Rscript R/studies/run_appendix_TN.R
Rscript R/studies/make_figures.R      # 重繪出版用圖
```

各腳本結尾的示範區塊包在 `if (sys.nframe() == 0)` 裡：
**`Rscript` 執行會跑示範，被 `source()` 則只載入函式**，可安全互相引用。

---

## 資料

| 檔案 | 用途 |
|---|---|
| `data/五齡Death.xlsx` | 死亡數（Year, Age, Female, Male, Total） |
| `data/5-age exposure.xls` | 暴露人數，同結構 |
| 其餘 `.txt` / `歷年人口資料.xlsx` | 單齡組原始資料，目前腳本未使用，保留供單齡版分析 |

讀檔由 `load_data()` 負責，回傳 `list(D, E, ages, years)`：

- **年齡標籤以「位置」重新指派**——exposure 檔的 `1-4`/`5-9`/`10-14` 被 Excel
  誤判成日期。
- 尾組併至 `100+`（`110+` 在 1970–2009 暴露數為 0），併組後 **22 個年齡組**。
- 路徑相對於 repo 根目錄，`load_data(sex = "Female")` 即可。

---

## 主要發現

**問題一：復刻的優勢有界限。** LASSO 修勻的配適誤差確實較小（MAE 0.0413 對
樣條的 0.0538／0.0759），$e_0^\dagger$ 調整亦使 $\kappa_t$ 更線性。但擴充為
「性別 × 起始年」12 格面板後，優勢僅在 5 格成立，其餘 7 格由 Lee-Miller 奪冠，
長配適期**無一例外**；且其預測區間是六種方法中**最寬**者。

**問題二：失效的機制是偏誤，不是變異數。** $\alpha_x$ 的偏誤在 $N=10^4$ 時達
$+0.439$、方向一致；漂移項在所有 $N$ 下一律偏向零。因此文獻共識
「LC 區間過窄」在小樣本**不成立**——區間其實過寬，涵蓋失敗來自位置偏誤。
以卜瓦松最大概似取代 SVD 可消除漂移項衰減，配合拔靴法去除 $\hat\sigma$ 的
估計雜訊，$N\ge10^6$ 時區間寬度幾乎追平理論下限。

**另一項發現：配適期長度比人口規模更關鍵。** $T=10$ 年時即使 $N=10^6$，
$\hat\beta$ 與真值的相關僅 0.120；$T=55$ 年時 $N$ 僅兩萬也達 0.746。
$T\le20$ 年時任何人口規模都無法達到 0.8 的門檻。

**問題三：懲罰可行，但收縮目標不能取零。** 在識別條件 $\sum_x\beta_x=1$ 之下，
$\sum_x|\beta_x| \equiv 1$ 是常數，**朝零收縮的懲罰項完全失效**（實測 $\lambda$
加到 $10^6$ 估計值不動）。改取均勻向量 $\mathbf{1}/A$ 為目標後，$N=5\times10^4$
時均方誤差降低 94%。惟純看均方誤差 **Ridge 優於 LASSO**（稀疏性並非正確的
先驗），LASSO 的價值在可解釋性。

各腳本結尾都有「先跑過的結果」註解，記錄數值與解讀，是讀這個 repo 最快的路徑。

---

## 與 Rabbi & Mazzuco (2021) 原文的差異

1. 原文用單齡組 0–100+，本 repo 是五齡組 22 組；年齡方向的二階差分改用
   不等距分割差商。
2. $e_0^\dagger$ 改用簡略生命表（Chiang）的離散近似。
3. 原文用 `smoothAPC`；本 repo 直接把式 (4)(5) 的堆疊設計矩陣丟給
   `quantreg::rq.fit.sfn`。
4. 原文比較的一維／二維樣條，本 repo 以**同一組差分算子的 L2 版本**代替，
   是乾淨的受控比較，但數值不能與原文逐項對照。
5. Fig 4–7、10 的面板維度是「性別 × 配適起始年」，不是原文的 20 個國家。
6. 留出期 2020–2024 涵蓋 COVID-19，正式版將另跑 2001–2019 對照。
7. HU / HUR / HUW 未復刻（需 `demography` 的函數型多主成分分解）。

---

## 注意事項

- **模擬規模已統一**為 `reps = 100`、`seed = 20260921`（`run_unified_mc.R`）。
  重跑後 λ 選法的相對損失由早期的「約七倍」修正為 3.49 倍、三方向格點的
  建議組合由 (1,2,1) 改為 (1,5,5)；結論方向不變。
- `output/` 下的檔案皆可由 `R/studies/` 重新產生。
- `參考資料/`（論文 PDF 與個人筆記）與 `範例論文/` 不納入版控。
