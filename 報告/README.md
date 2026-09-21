# 進度報告（2026-09-21）

給指導教授的書面報告與簡報。內容以 `參考資料/文獻回顧_問題與發現.md` 為底，
數據取自 `output/tables/`，圖取自 `output/figures/`。

| 檔案 | 說明 |
|---|---|
| `進度報告.tex` → `.pdf` | **詳細版（會議版）**，33 頁＋附錄（`ctexart`） |
| `進度報告_v2.tex` → `.pdf` | **第二版**，36 頁。在會議版之上增補最佳化演算法的層次與可用性分析（§柒之九） |
| `進度報告_簡版.tex` → `.pdf` | **簡版**，5 頁，摘錄問題、關鍵數據與結論 |
| `進度報告_簡報.tex` → `.pdf` | 簡報，32 張（`ctexbeamer`，16:9） |
| `建議文獻清單.md` | 依指導方向整理的文獻，標明已有與待下載 |

簡版供快速瀏覽，各節標註了詳細版的對應節號；詳細版含完整的實驗設計、
判讀邏輯與書目。第二版另含 §柒之九「點估計 LASSO 與分位數 LASSO 的
演算法差異」，該節的實測由 `R/studies/demo_algorithm_comparison.R` 產生。

## 編譯

需要 **XeLaTeX**（中文用 `ctex` + `fontset=macnew`，即 macOS 的 PingFang／Heiti）。
本機 TinyTeX 裝在 `~/Library/TinyTeX` 未加入 PATH，Makefile 已代為處理：

```sh
make            # 四份都編
make full       # 會議版（詳細）
make v2         # 第二版
make brief      # 簡版
make slides     # 簡報
make clean      # 清掉中間檔
```

> **版本關係**：`進度報告.tex` 是 2026-09-21 會議當天的版本，內容凍結；
> 後續新發現寫入 `進度報告_v2.tex`，兩者在會議前的內容相同。

首次編譯若缺套件（`ctex`、`xecjk`、`beamer`、`caption`、`pdflscape` 等），
因本機 TinyTeX 是 TeX Live 2025 而 CTAN 預設已是 2026，需指向凍結的 2025 檔案庫：

```sh
export PATH="$HOME/Library/TinyTeX/bin/universal-darwin:$PATH"
REPO="https://ftp.math.utah.edu/pub/tex/historic/systems/texlive/2025/tlnet-final"
tlmgr update --self --repository "$REPO"
tlmgr option repository "$REPO"
tlmgr install ctex xecjk beamer caption ...
```

## 體例

行文與格式比照 `範例論文/` 中余清祥老師的著作（2025 JPS、NAAJ 2019 等）：

- 章節用「壹貳參肆伍」層級，子節用「一、二、三」。
- 專有名詞首次出現標註原文，如「修勻（Graduation）」、「偏誤（Bias）」；
  縮寫（LC、SVD、GLM、RWD）一律在緒論首次出現時給全稱。
- 符號在第貳節統一定義：純量細體小寫、向量粗體小寫、矩陣粗體大寫；
  附下標的 $\beta_x$ 是分量，不附下標的 $\boldsymbol{\beta}$ 是整個向量。
- 引用格式：中文為「王信忠等（2012）」、英文為「Wang et al.（2018）」。
- 參考文獻分「一、中文部分」與「二、英文部分」。
- 圖表敘述依固定脈絡：呈現什麼 → 回應哪個問題 → 指標高低代表什麼 →
  觀察到什麼 → 與假說吻合與否 → 反映什麼現象 → 接續如何設計。

## 圖片來源

`\graphicspath` 指向 `../output/figures/`，直接取用 repo 既有輸出，不另存副本。
圖由 `R/studies/make_figures.R` 統一產生，該檔會套用
`R/core/fig_axis_utils.R` 的座標軸格式（把 `1e-03` 改為 $10^{-3}$）。

交互參照一律使用 `\label`/`\ref`，不要寫死表圖編號（唯一例外是正文提到
Rabbi & Mazzuco 原文的「原文圖 8」）。

## 待處理

- 書目已完整（卷期頁碼全部核對過）。**尚缺 PDF 的 19 筆**列於
  `建議文獻清單.md` 第二節，依重要性排序。
