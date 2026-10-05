# -*- coding: utf-8 -*-
"""圖 2 改為只用輔助線、不在圖上標文字；數值改由正文與圖註承擔。
並把證明環境的標籤由「證明」改為 Proof，與 Proposition／Corollary／Lemma 一致。"""
import io
P = "報告/進度報告_1005_V2.tex"
s = io.open(P, encoding="utf-8").read()
LOG = []
def sub(a, b, lab=""):
    global s
    n = s.count(a); assert n == 1, u"替換失敗（%d 次）：%s" % (n, (lab or a)[:44])
    s = s.replace(a, b); LOG.append(u"  " + lab)

# ---------- 一、圖 2 的圖註：補上輔助線的約定 ----------
sub(u"""{\\footnotesize 註：(a) $v(\\mu;c_0)$，於 $\\mu=2.04$ 取極大 $0.481$；(b) $A(\\mu)=1+\\mu b'(\\mu)$，於 $\\mu=4.43$ 取極大 $1.126$；(c) $1/A(\\mu)^2$，縱軸亦取對數，$\\mu\\to0$ 時發散、於 $\\mu=4.43$ 取極小 $0.789$；(d) 相對於 Poisson 最大概似的效率，於 $\\mu=4.47$ 取極小 $0.878$。}""",
u"""{\\footnotesize 註：圖中不標文字，輔助線的約定為長虛線標示極限或漸近線、
點線標示極值的位置與高度。各量的數值為：(a) $v(\\mu;c_0)$ 兩端趨於零，
於 $\\mu=2.0365$ 取極大 $0.4805$；(b) $A(\\mu)=1+\\mu b'(\\mu)$ 的兩端極限為
$0$ 與 $1$，於 $\\mu=4.4288$ 取極大 $1.1259$；(c) $1/A(\\mu)^2$ 的縱軸亦取對數，
$\\mu\\to0$ 時發散、$\\mu\\to\\infty$ 時趨於 $1$，於 $\\mu=4.4288$ 取極小 $0.7888$；
(d) 相對於 Poisson 最大概似的效率兩端皆趨近 $1$，
於 $\\mu=4.4697$ 取極小 $0.8776$。}""", u"圖 2 圖註補輔助線約定")

# ---------- 二、正文：說明圖的讀法 ----------
sub(u"""連同第\\ref{subsec:bprops}的偏誤 $b(\\mu;c_0)$，構成對該估計量的完整描述。
圖~\\ref{fig:v} 把其中四者並列。""",
u"""連同第\\ref{subsec:bprops}的偏誤 $b(\\mu;c_0)$，構成對該估計量的完整描述。
圖~\\ref{fig:v} 把其中四者並列。該圖不在圖面上標注文字，
各量的極限與漸近線以長虛線標示、極值的位置與高度以點線標示，
其數值如下：$v$ 的兩端極限為零、於 $\\mu=2.0365$ 取極大 $0.4805$；
$A$ 的兩端極限為 $0$ 與 $1$、於 $\\mu=4.4288$ 取極大 $1.1259$；
$1/A^2$ 在 $\\mu\\to0$ 時發散、在 $\\mu\\to\\infty$ 時趨於 $1$、
於 $\\mu=4.4288$ 取極小 $0.7888$；$\\mathrm{RE}$ 的兩端極限皆為 $1$、
於 $\\mu=4.4697$ 取極小 $0.8776$。""", u"圖 2 正文補讀法")

# ---------- 三、證明環境的標籤改為 Proof ----------
sub(u"""\\newenvironment{proof}{\\par\\noindent\\textit{證明}.\\ }{\\hfill$\\square$\\par\\vspace{4pt}}""",
u"""\\newenvironment{proof}{\\par\\noindent\\textit{Proof}.\\ }{\\hfill$\\square$\\par\\vspace{4pt}}""",
u"證明環境改為 Proof")

io.open(P, "w", encoding="utf-8").write(s)
print(u"\n".join(LOG)); print(u"\n已改寫 %s" % P)
