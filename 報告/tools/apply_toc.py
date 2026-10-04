# -*- coding: utf-8 -*-
"""加入目錄。封面、摘要與目錄編為羅馬頁碼，正文自第 1 頁起以阿拉伯頁碼重新
計數，故目錄不佔正文頁數。三個部次（第一部至第三部）一併進入目錄。"""
import io
P = "報告/進度報告_1005_V2.tex"
s = io.open(P, encoding="utf-8").read()
LOG = []
def sub(a, b, lab=""):
    global s
    n = s.count(a); assert n == 1, u"替換失敗（%d 次）：%s" % (n, (lab or a)[:44])
    s = s.replace(a, b); LOG.append(u"  " + lab)

# ---------- 一、部次一併寫入目錄 ----------
sub(u"""\\newcommand{\\partline}[2]{%
  \\par\\vspace{10pt}\\noindent\\rule{\\textwidth}{0.8pt}\\par\\vspace{2pt}
  \\noindent{\\zihao{4}\\bfseries 第#1部\\quad #2}\\par
  \\vspace{2pt}\\noindent\\rule{\\textwidth}{0.8pt}\\par\\vspace{6pt}}""",
u"""\\newcommand{\\partline}[2]{%
  \\par\\vspace{10pt}\\noindent\\rule{\\textwidth}{0.8pt}\\par\\vspace{2pt}
  \\noindent{\\zihao{4}\\bfseries 第#1部\\quad #2}\\par
  \\vspace{2pt}\\noindent\\rule{\\textwidth}{0.8pt}\\par\\vspace{6pt}%
  \\phantomsection\\addcontentsline{toc}{part}{第#1部\\quad #2}}

% 目錄：只收到小節，部次以粗體與一行間距與節次區隔。
\\setcounter{tocdepth}{2}
\\makeatletter
\\renewcommand{\\l@part}[2]{\\par\\vspace{8pt}%
  \\noindent{\\bfseries #1}\\par\\vspace{2pt}}
\\renewcommand{\\l@section}[2]{\\@dottedtocline{1}{0em}{5.2em}{\\bfseries #1}{\\bfseries #2}}
\\renewcommand{\\l@subsection}[2]{\\@dottedtocline{2}{5.2em}{4.6em}{#1}{#2}}
\\makeatother""", u"部次進目錄、目錄版式")

# ---------- 二、插入目錄，並把前置部分改為羅馬頁碼 ----------
sub(u"""\\begin{document}

\\begin{center}
{\\zihao{3}\\bfseries 小人口 Lee-Carter 參數偏誤的閉式}\\\\[4pt]""",
u"""\\begin{document}
\\pagenumbering{roman}

\\begin{center}
{\\zihao{3}\\bfseries 小人口 Lee-Carter 參數偏誤的閉式}\\\\[4pt]""",
u"前置改羅馬頁碼")

sub(u"""\\end{minipage}
\\end{center}

\\vspace{10pt}

\\section{緒論}\\label{sec:intro}""",
u"""\\end{minipage}
\\end{center}

\\clearpage
{\\zihao{4}\\bfseries 目\\quad 錄}\\par\\vspace{8pt}
\\makeatletter\\@starttoc{toc}\\makeatother

\\clearpage
\\pagenumbering{arabic}

\\section{緒論}\\label{sec:intro}""", u"插入目錄並重啟阿拉伯頁碼")

io.open(P, "w", encoding="utf-8").write(s)
print(u"\n".join(LOG)); print(u"\n已改寫 %s" % P)
