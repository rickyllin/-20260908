# -*- coding: utf-8 -*-
"""依 1004 建議對 V2 做結構與表格說明的調整（就地，只跑一次）。"""
import io, re
P = "報告/進度報告_1005_V2.tex"
s = io.open(P, encoding="utf-8").read()
LOG = []
def sub(a, b, lab=""):
    global s
    n = s.count(a); assert n == 1, u"替換失敗（%d 次）：%s" % (n, (lab or a)[:46])
    s = s.replace(a, b); LOG.append(u"  " + lab)

# ---------- 一、把「加權的最適形式」由柒節移到肆節 ----------
i = s.index(u"\\subsection{加權的最適形式}\\label{sec:weight}")
j = s.index(u"\\section{實證結果}\\label{sec:result}")
WBLOCK = s[i:j]
s = s[:i] + s[j:]
k = s.index(u"\\subsection{相對於 Poisson 最大概似的效率，與四個量的綜述}")
s = s[:k] + WBLOCK + s[k:]
LOG.append(u"  加權的最適形式由柒節移為肆之四（%d 字）" % len(WBLOCK))

sub(u"""式~\\eqref{eq:general} 的三項設定中，權重 $w$ 是唯一在本文最終建議中被實際
調動的一個，故其最適形式值得單獨推導。設損失為 $\\rho$、$\\psi=\\rho'$，""",
u"""權重是式~\\eqref{eq:general} 的三項設定中唯一在本文最終建議中被實際調動的
一個，而它與前三小節的三個量同屬二階：由三明治公式，最適權重由
$\\Ex[\\psi']$ 與 $\\Ex[\\psi^2]$ 決定。設損失為 $\\rho$、$\\psi=\\rho'$，""",
u"肆之四 開場改寫")

sub(u"""上一小節給出校正所付的變異數代價，但尚未涵蓋一個更基本的比較，
亦即留在對數尺度上作最小平方相對於直接在計數尺度上配適 Poisson 的效率損失。""",
u"""前三小節的三個量描述的是標準 Lee-Carter 與其校正各自的行為，
但尚未涵蓋一個更基本的比較，亦即留在對數尺度上作最小平方相對於直接在計數
尺度上配適 Poisson 的效率損失。""", u"肆之五 開場改寫")

sub(u"""標準 Lee-Carter 在對數尺度上的二階行為至此可由 $\\mu$ 的三個函數描述：
變異數 $v(\\mu;c_0)$、靈敏度 $A(\\mu)$ 與相對效率 $\\mathrm{RE}(\\mu)$，
連同校正的變異數代價 $1/A(\\mu)^2$ 與第\\ref{subsec:bprops}的偏誤
$b(\\mu;c_0)$，構成對該估計量的完整描述；圖~\\ref{fig:v} 把後四者並列。""",
u"""標準 Lee-Carter 在對數尺度上的二階行為至此可由 $\\mu$ 的五個函數描述：
變異數 $v(\\mu;c_0)$、靈敏度 $A(\\mu)$、校正的變異數代價 $1/A(\\mu)^2$、
最適權重 $w(\\mu)$ 與相對效率 $\\mathrm{RE}(\\mu)$；
連同第\\ref{subsec:bprops}的偏誤 $b(\\mu;c_0)$，構成對該估計量的完整描述。
圖~\\ref{fig:v} 把其中四者並列。""", u"肆之五 綜述改為五個量")

sub(u"""就本節與下一節的接續而言，四個量都只透過期望死亡數進入，""",
u"""就本節與下一節的接續而言，五個量都只透過期望死亡數進入，""",
u"肆之五 五個量")

# ---------- 二、表格說明：caption 與表註分工 ----------
sub(u"""\\caption{三種分布下的兩項與偏誤（$c_0=0.5$）。Poisson 與負二項皆取
$\\Ex[D]=\\mu$，負二項另取 $\\operatorname{Var}(D)=1.5\\mu$；
常態欄為式~\\eqref{eq:bnormal}，其 $S$ 恆為零}""",
u"""\\caption{三種分布下的 $S$ 項、Jensen 項與偏誤（$c_0=0.5$）}""",
u"tab:normal caption 縮短")
sub(u"""\\bottomrule
\\end{tabular}
\\label{tab:normal}
\\end{table}""",
u"""\\bottomrule
\\end{tabular}

{\\footnotesize Poisson 與負二項皆取 $\\Ex[D]=\\mu$，負二項另取
$\\operatorname{Var}(D)=1.5\\mu$；常態欄為式~\\eqref{eq:bnormal}，其 $S$ 恆為零。}
\\label{tab:normal}
\\end{table}""", u"tab:normal 加表註")

sub(u"""\\caption{零歲平均餘命的偏誤與離散。中位偏誤刻畫估計方程的中心，
中間四項刻畫其離散，均方根誤差與命中率則為兩者的綜合
（重複 100 次，單位：歲）}""",
u"""\\caption{配適末年零歲平均餘命的偏誤與離散（重複 100 次，單位：歲）}""",
u"tab:e0var caption 縮短")

# ---------- 三、伍節：驗證理由壓成一句、計算環境入腳註、女性理由縮短 ----------
sub(u"""Corollary~\\ref{cor:alpha} 是一個恆等式，因此嚴格說來並不需要驗證：
只要式~\\eqref{eq:pois} 的 Poisson 假設成立、且 $\\hat\\alpha_x$ 依
式~\\eqref{eq:svd} 定義為列平均，該式就必然成立。但這兩個前提在實際的估計
流程中未必如字面那樣單純。一方面，實際的 Lee-Carter 配適含奇異值分解與識別條件的
正規化，而 $\\sum_x\\hat\\beta_x=1$ 與 $\\sum_t\\hat\\kappa_t=0$ 這兩個約束是在
分解之後才施加的，其重新尺度化對 $\\hat\\alpha_x$ 的牽動，
單看推導並不明顯；另一方面，推導把 $\\beta_x\\kappa_t$ 當成已知的真值，
而實際上它們也是估計出來的，其誤差向 $\\hat\\alpha_x$ 的滲透，
同樣須以數值確認。本節因此不是在檢驗數學，而是在確認\\textbf{推導所描述的
那個量，與實際估計流程所產生的那個量為同一個量}。""",
u"""Corollary~\\ref{cor:alpha} 是恆等式，驗證的對象因而不是數學本身，
而是\\textbf{推導所描述的那個量與實際估計流程所產生的那個量是否為同一個}：
後者含奇異值分解、識別條件的正規化，以及 $\\hat\\beta_x\\hat\\kappa_t$ 的估計
誤差向 $\\hat\\alpha_x$ 的滲透。""", u"伍之二 驗證理由壓成一句")

sub(u"""\\textbf{計算環境。} 全部計算在 R 4.2.1（x86\\_64-apple-darwin17.0）下完成。
除讀取 Excel 檔的 \\texttt{readxl} 與稀疏矩陣的 \\texttt{Matrix} 外，
不依賴任何外部套件。本文所比較的各估計量（標準 Lee-Carter、Poisson 最大概似、
Firth 預先修正、中心化校正、資訊加權、partial SMR 與 Whittaker 比值修勻）
皆為自行實作，程式置於 \\texttt{R/core/} 與 \\texttt{R/studies/}，
各節所引的表格編號與產生該表的程式檔名一一對應。
奇異值分解採 R 內建的 \\texttt{svd}（LAPACK）， Poisson 配適採自行撰寫的
疊代重加權最小平方，收斂門檻為參數變動的最大絕對值小於 $10^{-10}$。""",
u"""\\textbf{計算環境。} 本文所比較的各估計量皆以 R 自行實作，不依賴外部套件。%
\\footnote{R 4.2.1（x86\\_64-apple-darwin17.0）；除讀取 Excel 檔的
\\texttt{readxl} 與稀疏矩陣的 \\texttt{Matrix} 外不依賴外部套件。
奇異值分解採內建的 \\texttt{svd}（LAPACK），Poisson 配適採自行撰寫的疊代重
加權最小平方，收斂門檻為參數變動的最大絕對值小於 $10^{-10}$。
程式置於 \\texttt{R/core/} 與 \\texttt{R/studies/}，
各節所引的表格編號與產生該表的程式檔名一一對應。}""",
u"計算環境入腳註")

sub(u"""全部分析只用女性，理由有二。首先，女性的死亡率曲線較平滑、偏離秩一雙線性
結構的程度較小，可使驗證聚焦於推導本身而不混入模型誤設。
其次，本文的閉式只透過期望死亡數進入偏誤，性別本身不是閉式中的變數，
男性的差異僅反映為期望死亡數水準的不同，故以女性為例不失一般性。""",
u"""全部分析只用女性。女性的死亡率曲線較平滑，可使驗證不混入模型誤設；
更要緊的是本文的閉式只透過期望死亡數進入偏誤，性別本身不是閉式中的變數，
男性的差異僅反映為期望死亡數水準的不同，故以女性為例不失一般性。""",
u"只用女性的理由縮短")

io.open(P, "w", encoding="utf-8").write(s)
print(u"\n".join(LOG)); print(u"\n已改寫 %s" % P)
