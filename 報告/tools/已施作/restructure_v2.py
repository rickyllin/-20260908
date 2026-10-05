# -*- coding: utf-8 -*-
"""把 進度報告_1005_V2 的章節結構對齊 V3（就地改寫，只跑一次）。

用法：於專案根目錄執行 `python3 報告/tools/restructure_v2.py`。
本檔所做的是結構與行文慣例的對齊：合併相鄰且論述連續的小節、
把理論推導改為「一般形式 → 三種分布 → 應用」、緒論的貢獻移入結論、
全文不使用 LC 簡稱，以及合併後相對指涉的修正。
執行過一次之後 V2 即為新結構，再跑會因 assert 失敗而中止，此為預期行為。
"""
import io, re, sys

src = "報告/進度報告_1005_V2.tex"
dst = src
s = io.open(src, encoding="utf-8").read()
LOG = []

def sub(a, b, label=""):
    global s
    n = s.count(a)
    assert n == 1, u"替換失敗（%d 次）：%s" % (n, (label or a)[:46])
    s = s.replace(a, b)

def cut(start_marker, end_marker, label):
    """刪去 [start_marker, end_marker) 之間的整段。"""
    global s
    i = s.index(start_marker)
    j = s.index(end_marker, i)
    LOG.append(u"  刪段 %-30s %5d 字" % (label, j - i))
    s = s[:i] + s[j:]

def merge(host, *gone):
    """把 gone 各小節的標題刪去，其 label 移到 host 小節的標題之後。"""
    global s
    m = re.search(r"\\subsection\{[^\n]*?\}\\label\{%s\}" % re.escape(host), s)
    assert m, u"找不到存留小節 %s" % host
    extra = ""
    for g in gone:
        mg = re.search(r"\n\\subsection\{[^\n]*?\}\\label\{%s\}\n" % re.escape(g), s)
        assert mg, u"找不到待併小節 %s" % g
        s = s[:mg.start()] + "\n" + s[mg.end():]
        extra += "\\label{%s}" % g
    m = re.search(r"\\subsection\{[^\n]*?\}\\label\{%s\}" % re.escape(host), s)
    s = s[:m.end()] + extra + s[m.end():]
    LOG.append(u"  併入 %-22s <- %s" % (host, u"、".join(gone)))

def move_block(start_marker, end_marker, before_marker, label):
    """把 [start_marker, end_marker) 的整段移到 before_marker 之前。"""
    global s
    i = s.index(start_marker)
    j = s.index(end_marker, i)
    block = s[i:j]
    s = s[:i] + s[j:]
    k = s.index(before_marker)
    s = s[:k] + block + s[k:]
    LOG.append(u"  移段 %-26s %5d 字" % (label, len(block)))


def drop_exhibit(label):
    """刪去含該 label 的 table 或 figure 環境。"""
    global s
    i = s.index("\\label{%s}" % label)
    for env in ("table", "figure"):
        b = s.rfind("\\begin{%s}" % env, 0, i)
        e = s.find("\\end{%s}" % env, i)
        if b != -1 and e != -1 and s.find("\\end{%s}" % env, b) == e:
            e += len("\\end{%s}" % env)
            while e < len(s) and s[e] == "\n":
                e += 1
            LOG.append(u"  刪 %-6s %-14s %5d 字" % (env, label, e - b))
            s = s[:b] + s[e:]
            return
    raise AssertionError(u"找不到 %s 所在的環境" % label)

# ===================== 併標題 =====================
LOG.append(u"\n併標題：")
# 肆：綜述與漸近說明併入相對效率一節
sub(u"""\\subsection{相對於 Poisson 最大概似的效率}\\label{subsec:releff}""",
    u"""\\subsection{相對於 Poisson 最大概似的效率，與四個量的綜述}\\label{subsec:releff}""",
    "肆之四 標題")
sub(u"""把前四個小節收束起來，標準 LC 在對數尺度上的二階行為可以由 $\\mu$ 的四個
函數完全描述：變異數 $v(\\mu;c_0)$、靈敏度 $A(\\mu)$、
校正的變異數代價 $1/A(\\mu)^2$，以及相對效率 $\\mathrm{RE}(\\mu)$。
四者連同第\\ref{subsec:bprops}的偏誤 $b(\\mu;c_0)$，構成對該估計量的完整描述。
圖~\\ref{fig:v} 把後四者並列。""",
u"""標準 LC 在對數尺度上的二階行為至此可由 $\\mu$ 的三個函數描述：
變異數 $v(\\mu;c_0)$、靈敏度 $A(\\mu)$ 與相對效率 $\\mathrm{RE}(\\mu)$，
連同校正的變異數代價 $1/A(\\mu)^2$ 與第\\ref{subsec:bprops}的偏誤
$b(\\mu;c_0)$，構成對該估計量的完整描述；圖~\\ref{fig:v} 把後四者並列。""",
"肆之五 開場")
merge("subsec:releff", "subsec:secsum", "subsec:asym")

# 伍
sub(u"""\\subsection{研究資料}\\label{subsec:data}""",
    u"""\\subsection{研究資料與模擬設定}\\label{subsec:data}""", "伍之一 標題")
sub(u"""第\\ref{sec:valid}節與第\\ref{sec:result}節的所有模擬共用同一組設定，為免重複，
此處一次寫定。資料來源、年齡分組、性別的選擇、真值的定義，
以及曝露數 $E_{x,t}=N\\,w_x$ 的人口規模與年齡結構設定均見前一小節，
以下只寫模擬的執行細節。""",
u"""第\\ref{sec:valid}節與第\\ref{sec:result}節的所有模擬共用同一組設定，
以下一次寫定其執行細節。""", "伍之二 開場")
sub(u"""模擬的重複次數、亂數種子、比較指標與計算環境見下一小節。""",
    u"""模擬的重複次數、亂數種子、比較指標與計算環境如下。""", "伍之一 末句")
merge("subsec:data", "sec:setup")
sub(u"""\\subsection{驗證的設計}\\label{subsec:vdesign}""",
    u"""\\subsection{驗證的設計與閉式預測的對照}\\label{subsec:vdesign}""", "伍之二 標題")
merge("subsec:vdesign", "sec:closedfit")
sub(u"""\\subsection{四個二階解析量的數值核對}\\label{sec:vcheck}""",
    u"""\\subsection{二階解析量的數值核對、逐年齡的變異數與校正代價的有效範圍}\\label{sec:vcheck}""",
    "伍之四 標題")
merge("sec:vcheck", "sec:vsd", "sec:vinfl")

# 陸
sub(u"""\\subsection{檢定的構造}\\label{subsec:feasconstr}""",
    u"""\\subsection{檢定的構造及其與信賴度標準的對照}\\label{subsec:feasconstr}""",
    "陸之一 標題")
merge("subsec:feasconstr", "subsec:cred2")
sub(u"""\\subsection{兩個門檻的實務用途}\\label{subsec:thresholds}""",
    u"""\\subsection{兩個門檻的實務用途與替代值 $c_0$ 的揭露}\\label{subsec:thresholds}""",
    "陸之三 標題")
merge("subsec:thresholds", "subsec:c0disc")

# 柒
sub(u"""\\subsection{修改後的估計流程}\\label{subsec:procedure}""",
    u"""\\subsection{修改後的估計流程，及其與既有修正方法的關係}\\label{subsec:procedure}""",
    "柒之一 標題")
sub(u"最後把本文的校正放回式~\\eqref{eq:general} 的分類中，以釐清它與既有補救的\n關係。",
    u"另須把本文的校正放回式~\\eqref{eq:general} 的分類中，以釐清它與既有補救的\n關係。",
    "柒之一 併段語氣")
move_block(u"另須把本文的校正放回式~\\eqref{eq:general} 的分類中",
           u"\\subsection{加權的最適形式}",
           u"\\subsection{參數偏誤的成分分解}",
           u"柒之一 與既有修正的關係")
merge("subsec:procedure", "subsec:relation")

# 捌
sub(u"""\\subsection{兩參數的聯合改善}\\label{sec:joint}""",
    u"""\\subsection{兩參數的聯合改善，與兩者在平均餘命上的分工}\\label{sec:joint}""",
    "捌之六 標題")
merge("sec:joint", "sec:split")

# 玖
sub(u"""\\subsection{本文修正的限制}\\label{sec:lim}""",
    u"""\\subsection{本文修正的限制、偏誤的方向，與四項須誠實面對的界限}\\label{sec:lim}""",
    "玖之一 標題")
merge("sec:lim", "sec:risk", "sec:bounds")
sub(u"""\\subsection{結論}\\label{subsec:concl}""",
    u"""\\subsection{結論與後續工作}\\label{subsec:concl}""", "玖之二 標題")
merge("subsec:concl", "subsec:future")



# ===================== 參之二、之三：一般形式 → 三種分布 → 應用 =====================
LOG.append(u"\n理論推導改序：")

NEW23 = io.open("報告/tools/已施作/v3/sec3.tex", encoding="utf-8").read()


i = s.index(u"\\subsection{估計方程的一致性檢驗}")
j = s.index(u"\\subsection{偏誤函數的性質}")
LOG.append(u"  參之二、之三改寫（%d 字 -> %d 字）" % (j - i, len(NEW23)))
s = s[:i] + NEW23 + s[j:]

# 柒之五 併入上述，整節刪去
cut(u"\\subsection{分布假設的敏感度分析}\\label{sec:normal}",
    u"\\section{實證結果}", u"柒之五 分布假設的敏感度分析")
sub(u"""式~\\eqref{eq:bform} 的形狀完全可由解析得出，而三項性質各對應一個實質的
結論。""",
u"""式~\\eqref{eq:bform} 的形狀完全可由解析得出，而兩個分支的消長各對應一個
實質的結論。""", "參之四 開場")



# ===================== 緒論：貢獻移入結論，改寫為研究目的與取徑 =====================
LOG.append(u"\n緒論：")

m = re.search(u"\\\\subsection\\{研究目的與貢獻\\}\\\\label\\{subsec:design\\}\n(.*?)(?=\\\\section\\{文獻回顧\\})", s, re.S)
assert m, u"找不到壹之三"
CONTRIB = m.group(1).strip()
NEWAIM = io.open("報告/tools/已施作/v3/aim.tex", encoding="utf-8").read()
s = s[:m.start()] + NEWAIM + s[m.end():]
LOG.append(u"  壹之三 貢獻（%d 字）移入結論，改寫為研究目的與取徑" % len(CONTRIB))

sub(u"""簡言之，圖~\\ref{fig:osc} 所呈現的是：若不處理偏誤，鄉鎮市區生命表所顯示的
地區差異將有相當部分是估計方法的產物，而非真實的健康差異。""",
u"""簡言之，圖~\\ref{fig:osc} 所呈現的是：若不處理偏誤，鄉鎮市區生命表所顯示的
地區差異將有相當部分是估計方法的產物，而非真實的健康差異。

此一偏誤的實務後果在於地區之間的比較。鄉鎮市區生命表的用途多為排序與資源
配置，而排序只需要各地區之間的相對位置正確；但偏誤的量隨人口規模而變，
人口較少的地區偏移較大，故各地區的相對位置會被系統性地改變，
而該改變來自估計方法而非真實的健康差異。
加大人口規模並非可行的補救，蓋因鄉鎮市區的規模由行政區劃決定；
延長觀測期間亦不可行，理由見第\\ref{subsec:consistency}。""",
"壹之一 補實務後果")

# 貢獻移入結論（置於結論末、後續工作之前）
sub(u"""參考母體修勻的限制則在第\\ref{sec:smr}得到定量的確認。""",
u"""參考母體修勻的限制則在第\\ref{sec:smr}得到定量的確認。

概括而論，本文的貢獻有三項。其一，標準 Lee-Carter 的年齡參數 Fisher 不一致，
而該不一致有一般形式，其在 Poisson 之下為閉式；$\\hat\\alpha_x$ 的偏誤是一個
有限樣本的恆等式而非漸近近似，且不隨觀測年數增加而消失。
其二，本文導出三個二階解析量，亦即對數尺度的變異數、校正所付的變異數代價，
以及相對於 Poisson 最大概似的效率；三者合起來界定了這一類校正所能達到的
上限，其中效率的損失全程不超過 $12.2\\%$。
其三，由於偏誤的計算只需要曝露數與一組粗略的死亡率，本文把它包裝成一個
事前的可行性檢定；此一檢定補上了古典信賴度標準所缺的一半，
蓋因後者控制估計的變異而假設估計量無偏，而該假設在本問題上不成立。
本文的範圍則限於年齡水準參數 $\\alpha_x$ 的偏誤，其技術理由在於該參數在
雙線性部分被估計之前即已由列平均定出，因而可以完全在單變數的框架內分析；
時間參數的偏誤與雙線性結構的識別問題不在本文的涵蓋範圍內。""",
"貢獻移入結論")


# ===================== 不使用 LC 簡稱 =====================
n0 = s.count(u"LC")
s = s.replace(u"\\newcommand{\\LCe}{\\mathrm{LC}e^{\\dagger}}", u"@@XMACROX@@")
s = s.replace(u"LC", u"Lee-Carter")
s = s.replace(u"@@XMACROX@@", u"\\newcommand{\\LCe}{\\mathrm{LC}e^{\\dagger}}")
s = s.replace(u"Lee-Carter Lee-Carter", u"Lee-Carter")
s = s.replace(u"Lee-Lee-Carter", u"Lee-Carter")
LOG.append(u"\n不使用簡稱：LC -> Lee-Carter（%d 處）" % n0)


# ===================== 末批四：合併標題後的相對指涉 =====================
LOG.append(u"\n末批四（相對指涉）：")

sub(u"""前三小節的驗證針對一階的閉式，本小節與其後兩小節則針對
第\\ref{sec:second2}節的四個二階解析量。
其中 Proposition~\\ref{prop:v}、Proposition~\\ref{prop:bp} 與
Corollary~\\ref{cor:re} 為恆等式，嚴格說來不需要檢驗，
但三式都含無窮級數，截斷的位置可能影響結果，故先以模擬核對其實作；
Proposition~\\ref{prop:infl} 的代價含一階展開，其有效範圍須以數值界定，
於第\\ref{sec:vinfl}處理。""",
u"""前三小節的驗證針對一階的閉式，本小節則針對第\\ref{sec:second2}節的四個
二階解析量。其中 Proposition~\\ref{prop:v}、Proposition~\\ref{prop:bp} 與
Corollary~\\ref{cor:re} 為恆等式，嚴格說來不需要檢驗，
但三式都含無窮級數，截斷的位置可能影響結果，故先以模擬核對其實作；
Proposition~\\ref{prop:infl} 的代價含一階展開，其有效範圍則於本小節末以
數值界定。""", "伍之四 開場")

sub(u"""Proposition~\\ref{prop:infl} 的代價含一階展開，故其有效範圍須以數值界定。
表~\\ref{tab:v3} 把三個人口規模的 22 個年齡依期望死亡數分段彙總。""",
u"""表~\\ref{tab:v3} 把三個人口規模的 22 個年齡依期望死亡數分段彙總，
比較校正代價的預測與實測。""",
"伍之四 去重複句")

sub(u"""本小節把前述的解析結果翻成這個角度，所用的數字全部來自""",
    u"""以下把前述的解析結果翻成這個角度，所用的數字全部來自""", "玖之一 本小節")
sub(u"""除前一小節所述的限制外，另有四項界限須主動寫明，""",
    u"""除前述三項限制外，另有四項界限須主動寫明，""", "玖之一 除前一小節")
sub(u"""前一小節已指出該項的方向對年金不審慎，而本文的修正並未改善它。""",
    u"""前文已指出該項的方向對年金不審慎，而本文的修正並未改善它。""", "玖之一 前一小節")

s = s.replace(u"\\newcommand{\\LCe}{\\mathrm{LC}e^{\\dagger}}\n", u"")


# 節次層級的 label 若漏了「節」字一併補上（V2 沿襲下來的瑕疵）
SECLABS = ["sec:intro", "sec:lit", "sec:method", "sec:second2", "sec:valid",
           "sec:feas", "sec:corr", "sec:result", "sec:disc"]
_n = 0
for _lab in SECLABS:
    _pat = re.compile(r"第(\\ref\{%s\})(?!\s*節)" % re.escape(_lab))
    s, _k = _pat.subn(r"第\g<1>節", s)
    _n += _k
LOG.append(u"  補上節次參照的「節」字（%d 處）" % _n)





io.open(dst, "w", encoding="utf-8").write(s)
print("\n".join(LOG))
print("\n已寫出 %s" % dst)



io.open(dst, "w", encoding="utf-8").write(s)
print("\n".join(LOG))
print("\n已改寫 %s" % dst)
