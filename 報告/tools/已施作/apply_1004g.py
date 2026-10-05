# -*- coding: utf-8 -*-
"""表圖標題改為簡潔形式、長說明移入表／圖註；並把宣告結論的標題改為主題式標題。
（就地改寫 V2，只跑一次）"""
import io, re
P = "報告/進度報告_1005_V2.tex"
s = io.open(P, encoding="utf-8").read()
LOG = []

def caption_of(label):
    """回傳 (起, 迄, 內容) 指向該 label 所屬環境的 \caption{...}。"""
    i = s.index(u"\\label{%s}" % label)
    b = max(s.rfind(u"\\begin{table}", 0, i), s.rfind(u"\\begin{figure}", 0, i))
    c = s.index(u"\\caption{", b)
    j = c + len(u"\\caption{"); d = 1
    while d:
        if s[j] == '{': d += 1
        elif s[j] == '}': d -= 1
        j += 1
    return c, j, s[c + len(u"\\caption{"):j - 1]

def fix(label, newcap, note=None):
    """縮短 caption，並把長說明放到表／圖註。"""
    global s
    c, j, old = caption_of(label)
    s = s[:c] + u"\\caption{%s}" % newcap + s[j:]
    if note:
        i = s.index(u"\\label{%s}" % label)
        b = max(s.rfind(u"\\begin{table}", 0, i), s.rfind(u"\\begin{figure}", 0, i))
        istab = s.rfind(u"\\begin{table}", 0, i) > s.rfind(u"\\begin{figure}", 0, i)
        if istab:
            k = s.rindex(u"\\end{tabular}", 0, i) + len(u"\\end{tabular}")
        else:
            k = s.index(u"\\caption{", b)
            k = s.index(u"}\n", k) + 1
        s = s[:k] + u"\n\n{\\footnotesize 註：%s}\n" % note + s[k:]
    LOG.append(u"  %-14s %d -> %d 字" % (label, len(old), len(newcap)))

# ---------- 圖 ----------
fix(u"fig:osc", u"小人口下 Lee-Carter 估計的震盪與偏誤（臺灣女性，以全國配適為真值）",
    u"(a) 單次實現的年齡別死亡率，空心點為零死亡格；"
    u"(b) $\\hat\\beta_x$ 的 20 次重複；"
    u"(c) 零歲平均餘命的 100 次重複，虛線為真值。")
fix(u"fig:v", u"四個二階解析量對期望死亡數的變化（$c_0=0.5$，橫軸取對數）",
    u"(a) $v(\\mu;c_0)$，於 $\\mu=2.04$ 取極大 $0.481$；"
    u"(b) $A(\\mu)=1+\\mu b'(\\mu)$，於 $\\mu=4.43$ 取極大 $1.126$；"
    u"(c) $1/A(\\mu)^2$，縱軸亦取對數，$\\mu\\to0$ 時發散、於 $\\mu=4.43$ 取極小 $0.789$；"
    u"(d) 相對於 Poisson 最大概似的效率，於 $\\mu=4.47$ 取極小 $0.878$。")
fix(u"fig:aband", u"$\\hat\\alpha_x$ 的偏誤：5--95\\% 帶與閉式預測",
    u"黑色虛線為 Corollary~\\ref{cor:alpha} 的預測值 $\\bar b_x$，"
    u"與標準 Lee-Carter 的中線最大差 $0.043$ 與 $0.019$。")
fix(u"fig:feas", u"事前可行性檢定的輸出",
    u"(a) 四個人口規模下逐年齡的預測偏誤 $\\bar b_x$；"
    u"(b) 分級判準（期望死亡數低於 $\\mu^\\ast$ 的格子比例）隨人口規模的變化與四級門檻。")
fix(u"fig:bband", u"$\\hat\\beta_x$ 的 5--95\\% 帶（上列人數一萬、下列人數五萬）",
    u"每格一個估計量，黑線為真值。")
fix(u"fig:cmp", u"三項指標隨人口規模的變化（縱軸取對數）",
    u"實線與實心點為 Poisson 模擬，虛線與空心點為實際觀測死亡數的 Poisson 抽薄。")
fix(u"fig:e0", u"配適末年零歲平均餘命的 100 次重複（人數五萬）",
    u"虛線為真值 $84.07$ 歲，各估計量上方標示其中位偏誤。")
fix(u"fig:ageprof", u"逐年齡的均方根誤差（縱軸取對數）",
    u"上列為人數一萬、下列為人數五萬；左欄為 $\\alpha_x$、右欄為 $\\beta_x$。"
    u"右欄中標準 Lee-Carter（紅）的線完全被中心化（藍）覆蓋，"
    u"因為僅校正 $\\alpha_x$ 不改變中心化矩陣。")
fix(u"fig:e0osc", u"零歲平均餘命的逐年序列（單次實現，亂數種子固定）",
    u"深紅為觀測序列、黑為真值、淺紅虛線為標準 Lee-Carter。")
fix(u"fig:pareto", u"$\\alpha$ 與 $\\beta$ 的取捨圖",
    u"橫軸為 SSE($\\beta$)、縱軸為 $\\alpha$ 最大偏誤，皆取對數，愈靠左下愈好；"
    u"箭頭由標準 Lee-Carter 出發。")
fix(u"fig:jointband",
    u"標準 Lee-Carter、Firth 與僅校正 $\\alpha$ 加權的逐年齡分布（人數五萬，5--95\\% 帶）",
    u"(a) 逐年齡偏誤；(b) $\\hat\\beta_x$ 的帶。")

# ---------- 表 ----------
fix(u"tab:wcases", u"三種損失的變異數最適權重",
    u"$\\beta_z=\\Ex[\\min(Z^2,z^2)]$，$Z\\sim N(0,1)$；$f$ 為 $u$ 的密度。")
fix(u"tab:popstr", u"人口年齡結構對偏誤的影響（人口總數固定於五萬，重複 100 次）",
    u"三種結構的真實死亡率相同。")
fix(u"tab:v1", u"四個二階解析量的值與數值核對（$c_0=0.5$）",
    u"$v$ 以 $2\\times10^{6}$ 次 Poisson 抽樣的樣本變異數為對照；"
    u"$b'$ 與步長 $10^{-5}\\mu$ 的中央差分比值在全部 15 個 $\\mu$ 上皆為 $1.0000$，故不另列。")
fix(u"tab:ageband", u"逐年齡段的均方根誤差與各段的平均期望死亡數（人數五萬，重複 100 次）")
fix(u"tab:e0osc", u"零歲平均餘命逐年序列的震盪幅度（100 次重複的中位數，單位：歲）",
    u"粗糙度取二階差分以濾去線性趨勢而只留不規則成分，"
    u"全距則用以判讀趨勢的幅度是否被壓縮。")
fix(u"tab:smr", u"與參考母體修勻法的比較（人數五萬，參考母體兩百萬，重複 100 次）",
    u"三個不使用參考母體的估計量與三個使用者並列，"
    u"後者分別在型態相同與型態不同的情境下運作。")

# ---------- 標題：宣告結論者改為主題式 ----------
TITLES = [
 (u"\\section{標準 Lee-Carter 的 Fisher 不一致與偏誤閉式}",
  u"\\section{取對數與零格替代所造成的偏誤}", u"參 節名"),
 (u"\\subsection{以外部資訊換取穩定：多母體、參考母體與修勻}",
  u"\\subsection{多母體、參考母體與修勻}", u"貳之二"),
 (u"\\subsection{以改動估計機制換取穩定：對數轉換、零格替代與有限樣本偏誤的修正}",
  u"\\subsection{對數轉換、零格替代與有限樣本偏誤的修正}", u"貳之三"),
 (u"\\subsection{可供借用的一般理論：$M$-估計、混合表示、損失與懲罰的分類，與信賴度標準}",
  u"\\subsection{$M$-估計、混合表示、損失與懲罰的分類，與信賴度標準}", u"貳之四"),
 (u"\\subsection{相對於 Poisson 最大概似的效率，與四個量的綜述}",
  u"\\subsection{相對於 Poisson 最大概似的效率}", u"肆之五"),
 (u"\\subsection{二階解析量的數值核對、逐年齡的變異數與校正代價的有效範圍}",
  u"\\subsection{二階解析量的數值核對}", u"伍之四"),
 (u"\\subsection{逐年齡的改善剖面}", u"\\subsection{逐年齡的誤差剖面}", u"捌之二"),
 (u"\\subsection{兩參數的聯合改善，與兩者在平均餘命上的分工}",
  u"\\subsection{校正與加權的施加位置，及其在平均餘命上的分工}", u"捌之六"),
 (u"\\subsection{本文修正的限制、偏誤的方向，與四項須誠實面對的界限}",
  u"\\subsection{修正的限制、偏誤的方向與適用的界限}", u"玖之一"),
]
for a, b, lab in TITLES:
    assert s.count(a) == 1, u"標題替換失敗：%s" % lab
    s = s.replace(a, b); LOG.append(u"  標題 %s" % lab)

io.open(P, "w", encoding="utf-8").write(s)
print(u"\n".join(LOG)); print(u"\n已改寫 %s" % P)
