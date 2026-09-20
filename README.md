# 台灣死亡率 Lee-Carter 模型：Rabbi & Mazzuco (2021) 復刻與延伸

碩士論文的分析程式碼。以台灣五齡組死亡資料（1970–2024）為基礎，做三件事：

1. **復刻** Rabbi & Mazzuco (2021, *European Journal of Population* 37:97–120) 的
   `LC_e0†` 方法——先用 2D L1（quantile LASSO）修勻死亡率曲面，再跑標準 Lee-Carter，
   最後調整 κ_t 使配適的壽命差異 e0† 吻合觀察值。
2. **延伸**：把懲罰從「曲面」搬到「參數」——帶 L1 懲罰的 Poisson Lee-Carter，
   β_x 朝參考曲線收縮。
3. **模擬檢驗**：蒙地卡羅研究暴露數大小對估計的影響、λ 的選擇準則比較、
   以及重抽樣的參數／預測不確定性。

---

## 快速開始

```r
# 環境：R 4.2+，需要 readxl, quantreg, SparseM, Matrix
install.packages(c("readxl", "quantreg", "SparseM", "Matrix"))
```

所有腳本都用**相對路徑** `source()` 彼此，**工作目錄必須是本 repo 根目錄**
（用 `碩士論文＿20260908.Rproj` 開啟 RStudio 即可）。

```sh
Rscript "論文復刻.R"          # 主要成果：全部圖表輸出到 rm_figs/
Rscript "Run female 2001 2024 .R"   # Female 2001-2024 的單次分析
Rscript "Mc exposure.R"       # 蒙地卡羅（較慢，數分鐘）
```

每個腳本結尾的示範區塊都包在 `if (sys.nframe() == 0)` 裡：
**用 `Rscript` 執行會跑示範；被 `source()` 進來則只載入函式**，可以安全地互相引用。

---

## 資料

| 檔案 | 用途 |
|---|---|
| `data/五齡Death.xlsx` | 死亡數（Year, Age, Female, Male, Total），24 個年齡組 × 55 年 |
| `data/5-age exposure.xls` | 暴露人數，同結構 |
| `data/Death.txt`、`Population size.txt`、`Exposure To Risk.txt`、`曝險人數資料.txt`、`歷年人口資料.xlsx` | 單齡組原始資料，**目前的腳本都沒有讀**，保留供未來單齡版分析用 |

讀檔由 `lc_poisson_lasso.R` 的 `load_data()` 負責，回傳 `list(D, E, ages, years)`：

- **年齡標籤以「位置」重新指派**，不依賴檔案內的 Age 文字——exposure 檔的
  `1-4` / `5-9` / `10-14` 被 Excel 誤判成日期。
- 預設把 `100-104`、`105-109`、`110+` 併成開放組 `100+`（`110+` 在 1970–2009
  暴露數為 0），併組後為 **22 個年齡組**。

`load_data()` 的預設路徑是**相對於 repo 根目錄**的 `data/`，所以工作目錄設對就能直接用：

```r
load_data(sex = "Female")                    # 走預設的 data/
load_data(sex = "Female", path_death = "…")  # 要換資料來源才需明寫
```

---

## 檔案結構與相依關係

```
lc_poisson_lasso.R              ← 最底層：load_data(), lc_poisson(), lc_path()
  └── rabbi_mazzuco_replication.R   ← 核心工具箱
        ├── Patch adjust kappa.R        （必須在主檔之後 source，會覆寫 adjust_kappa）
        ├── 論文復刻.R                   ← 主要成果
        ├── Run female 2001 2024 .R
        └── Mc exposure.R
              ├── Lambda selection.R
              └── lambda3_and_resampling.R
```

### 各檔案內容

**`lc_poisson_lasso.R`** — 讀檔 + 帶 L1 懲罰的 Poisson Lee-Carter。
`D_xt ~ Poisson(E_xt · exp(α_x + β_x κ_t))`，對 `|β_x − β0_x|` 加懲罰，
以 proximal Newton + 軟閾值 + 二分法滿足 `Σβ = 1`。
另有 `lc_lambda_max()`（glmnet 式的 λ_max）與 `lc_path()`（λ 路徑 + AIC/BIC）。

> 檔尾有一段重要的恆等式說明：在 `Σβ_x = 1` 的識別條件下，
> `β0 = 0` 時 L1 懲罰項**恆為常數**（只罰 β 的負部），估計值一個位元都不會動。
> 要讓 L1 真的產生選擇效果，收縮目標必須是非零的參考曲線。

**`rabbi_mazzuco_replication.R`** — 核心工具箱，其他檔案都從這裡取用：

| 函式 | 功能 |
|---|---|
| `make_grid()`, `D1_uneven()`, `D2_uneven()` | 五齡組的不等距差商算子（以組中點為節點） |
| `smooth_l1_2d()` | 2D L1 修勻（原文式 (4)(5)），堆疊成 τ=0.5 的分位數迴歸 |
| `life_table()` | 簡略生命表（Chiang），回傳 `e0` 與 `edag` |
| `lc_svd_fit()` | 標準 LC 的 SVD 一階段估計 |
| `adjust_kappa()` | 逐年解 κ_t 使某目標量吻合觀察值 |
| `fit_variant()` | LC / LM / LCP / LCedag 四種變體 |
| `oos_evaluate()` | 留出後 h 年的樣本外評估（含 jump-off 修正） |

**`Patch adjust kappa.R`** — 修正 `adjust_kappa`：e0† 對 κ 是**單峰而非單調**，
原版的寬區間 + `extendInt` 會兩端同號而失敗。改成網格掃描 + 取最接近原始 κ 的根。
也定義了 `safe_range()`（濾掉非有限值的 `range`）與診斷用的 `plot_target_curve()`。
**一定要在 `rabbi_mazzuco_replication.R` 之後 source。**

**`論文復刻.R`** — 主要成果。復刻原文所有圖表，輸出到 `rm_figs/`：

| 函式 | 對應原文 | 輸出 |
|---|---|---|
| `fig1_smoothing()` | §3.1 Fig 1 | 三種修勻法（L1 2D / L2 2D / L2 1D）比較 |
| `fig2_params()` | §3.2 Fig 2 | α_x 與 β_x |
| `fig3_kappa()` | §3.2 Fig 3 / Table 1 | κ_t 與 random walk drift |
| `fig8_e0_edag()` | §3.2 Fig 8 | e0 與 e0† 的關係 |
| `fig4to7_accuracy()` | §3.3 Fig 4–7 | 樣本外準確度面板（**性別 × 配適起始年**，非原文的 20 國） |
| `fig9_e0_forecast()` | §3.4 Table 2 / Fig 9 | 2050 年 e0 預測與區間 |
| `fig10_11_coverage()` | §3.4 Fig 10–11 | 區間預測涵蓋率與誤判年數 |

比較的六種方法：`LC`（對齊總死亡數）、`LCP`（Poisson ML）、`LM`（Lee-Miller，對齊 e0）、
`BMS`（簡化版，搜尋最佳配適起始年）、`LCedag`（本文復刻，修勻後對齊 e0†）、
`MLCedag`（修勻後的期望死亡數跑 Poisson，再對齊 e0†）。
**HU / HUR / HUW 未復刻**（需 `demography` 套件的函數型多主成分分解）。

**`Mc exposure.R`** — 蒙地卡羅：以 Female 2001–2024 的 LC 配適為真值，
`E = N · w`、`D ~ Poisson(E·m)`，比較純 LC 與 LASSO 修勻後的 LC，
評估 α/β/κ 的偏誤、變異數、MSE、drift 與崩潰率 `P(cor(β̂, β) < 0.8)`。

**`Lambda selection.R`** — 四種 λ 選法與 oracle 對照：
(A) 原文式（對觀察值的 MAE/MSE）、(B) 留出格交叉驗證、
(C) 留出年份的下游預測誤差、(D) 參數式 bootstrap 插入式選法。

**`lambda3_and_resampling.R`** — 兩部分：
(1) 把 λ 拆成三個方向 `(lxx, lxt, ltt)` 分開做格點搜尋；
(2) 重抽樣模組——參數不確定性（`poisson` 參數式 / `pearson` 半參數式）、
預測不確定性（κ 創新項的 `iid` / `block` 重抽）、以及
`resample_calibration()` 檢查重抽樣 SE 與蒙地卡羅真實 SD 的比值。

**`m_x,t.R`** — 早期的獨立腳本：Lee-Carter 從讀檔到預測的完整流程。
與主線不共用函式，輸出寫進 `data/`
（`lc_params_ax_bx.csv`、`lc_kappa_t.csv`、`lc_kappa_forecast.csv`）。
保留作為對照，新的分析請走 `rabbi_mazzuco_replication.R` 這條線。

---

## 輸出檔案

| 檔案 | 產生者 |
|---|---|
| `rm_figs/*.png`、`rm_figs/table*.csv`、`accuracy_panel.csv` | `論文復刻.R` |
| `mc_results.csv`、`mc_plots.png` | `Mc exposure.R` |
| `lambda_selection.csv` | `Lambda selection.R` |
| `lambda_grid3.csv` | `lambda3_and_resampling.R`（執行後才產生，repo 內尚無） |
| `lc_poisson_lasso_params.csv`、`lc_poisson_lasso_path.csv` | `lc_poisson_lasso.R` |
| `rabbi_mazzuco_oos.csv` | `rabbi_mazzuco_replication.R` |
| `data/lc_params_ax_bx.csv`、`lc_kappa_t.csv`、`lc_kappa_forecast.csv` | `m_x,t.R` |

---

## 已知結果與注意事項

每個腳本的**檔尾都有一段「先跑過的結果」註解**，記錄數值與解讀，是讀這個 repo 最快的路徑。
幾個關鍵發現：

- **e0† 調整讓 κ 變得比原始 LC 更不線性、更抖**（離線性殘差 sd：LC 0.849 → LCedag 1.420），
  與原文在瑞典資料上的結論**相反**。兩個候選解釋（五齡組使 e0† 的離散近似誤差偏大／
  台灣的死亡轉型較快且不規則）可用單齡組資料重跑來分離。
- **原文式的 λ 選法（對觀察值的 MAE/MSE）必然選到 λ = 0**——修勻一定會增加對觀察值的誤差。
  這是準則定義決定的，不是方法比較。
- **`lxt`（年齡×時間交叉項）是修勻的主力**。在 LC 底下
  `∂²log m/∂x∂t = β'_x · κ'_t`，所以交叉項本質上就是對 β_x 粗糙度的懲罰，
  是三個方向裡唯一直接瞄準 β 的。純年齡方向的修勻 `(20,0,0)` 比完全不修勻還糟。
- **小人口下 drift 被系統性壓平**（低估死亡改善速度）。這是有方向的偏誤，
  比變異數變大嚴重；LASSO 修勻在 `N ≥ 5e4` 時能把 drift 拉回真值附近且四分位距減半。
- **N ≤ 2e4 時零格佔四成**，主導一切的是零格處理（`pmax(D, 0.5)`）而非估計方法。
- 負面結果（已試過，不要再走）：把零格排除在保真項之外的加權版本**更差**——
  丟掉零格等於丟掉「此處死亡率很低」的訊息。

### 與原文的差異（報告時務必寫明）

1. 原文用單齡組 0–100+，本 repo 是五齡組 22 組；年齡方向的二階差分改用不等距分割差商。
2. e0† 改用簡略生命表（Chiang）的離散近似，精度低於單齡版本。
3. 原文用 `smoothAPC` 套件；本 repo 直接把式 (4)(5) 的堆疊設計矩陣丟給
   `quantreg::rq.fit.sfn`，結果相同但可控。
4. 原文比較的一維／二維樣條，本 repo 以**同一組差分算子的 L2 版本**代替，
   是乾淨的受控比較，但數值不能與原文逐項對照。
5. Fig 4–7、10 的面板維度是「性別 × 配適起始年」，不是原文的 20 個國家。
6. 留出期 2020–2024 涵蓋 COVID（女性 e0 在 2020 達 84.30，2022 降至 83.14，
   2024 才回到 83.94），所有方法的 ME(e0) 因此一致為負。
   要乾淨的比較，另跑 `Y1 = 2019` 的版本當對照。

### 其他坑

- 檔名含空格與中文，`source()` 時記得加引號。
- 腳本內的圖表標題用中文，在某些 graphics device 上會變成亂碼；
  `Patch adjust kappa.R` 的部分修正改用 ASCII 標題。
- `Mc exposure.R`、`Lambda selection.R`、`lambda3_and_resampling.R` 的示範區塊
  都用了**偏小的 `reps`** 以便試跑流程；正式結果請把 `reps` 調到 30 以上。
