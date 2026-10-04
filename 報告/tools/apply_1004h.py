# -*- coding: utf-8 -*-
"""連接詞的 register 調整（就地改寫 V2，只跑一次）。
依余老師團隊五篇的用語分布：「故」與「因而」屬文言、該文類不用，一律改「因此」；
命題與證明之內保留「故」，因數學敘述的慣例如此。
"""
import io, re
P = "報告/進度報告_1005_V2.tex"
s = io.open(P, encoding="utf-8").read()
N0 = {w: s.count(w) for w in [u"故", u"因而", u"因此", u"該文", u"綜上所述",
                              u"概括而論", u"其次", u"再者", u"而", u"但", u"然而"]}

# ---------- 一、因而 -> 因此 ----------
s = s.replace(u"因而", u"因此")

# ---------- 二、故 -> 因此（命題與證明之內保留）----------
mask = [False] * len(s)
for env in ["proof", "prop", "cor", "lem"]:
    for m in re.finditer(r"\\begin\{%s\}" % env, s):
        e = s.find(r"\end{%s}" % env, m.end())
        if e > 0:
            for k in range(m.start(), e): mask[k] = True
out = []; i = 0
while i < len(s):
    if s[i] == u"故" and not mask[i]:
        out.append(u"因此"); i += 1
    else:
        out.append(s[i]); i += 1
s = u"".join(out)

# ---------- 三、該文 -> 具名或代換 ----------
for a, b in [(u"而機制不同：該處的偏誤由誤差的分布決定", u"而機制不同：該處的偏誤由誤差的分布決定"),
             (u"該文所處理的廣義線性混合模型", u"其所處理的廣義線性混合模型"),
             (u"該文並據此給出一個三難結果", u"並據此給出一個三難結果"),
             (u"該文的 $a$ 是反應變數的計量單位", u"該命題的 $a$ 是反應變數的計量單位"),
             (u"該文扣的是最大概似估計量的", u"後者扣的是最大概似估計量的"),
             (u"\\citet{lee2003} 亦即據此主張該法", u"\\citet{lee2003} 據此主張該法")]:
    if a != b and s.count(a) == 1: s = s.replace(a, b)

# ---------- 四、「綜上所述／概括而論」各留一次 ----------
def thin(word, keep, alts):
    global s
    n = s.count(word); k = 0; idx = 0; res = []
    while True:
        j = s.find(word, idx)
        if j < 0: res.append(s[idx:]); break
        res.append(s[idx:j])
        if k < keep: res.append(word)
        else: res.append(alts[(k - keep) % len(alts)])
        k += 1; idx = j + len(word)
    s = u"".join(res)

thin(u"綜上所述", 1, [u"總結以上", u"由以上各項", u"合而觀之"])
thin(u"概括而論", 1, [u"總體而言", u"整體來看", u"就全局而言"])

# ---------- 五、「其次／再者」部分改為帶邏輯的轉折 ----------
thin(u"其次，", 6, [u"另一方面，", u"與此相對，", u"同時，", u"此外，",
                    u"另一項是，", u"接著，", u"又，", u"並且，"])
thin(u"再者，", 3, [u"此外，", u"另外，", u"同時，", u"又，"])

# ---------- 六、「而」的機械削減：對比與同位處改用分號 ----------
for a, b in [(u"，而後者", u"；後者"), (u"，而該", u"；該"), (u"，而這", u"；這"),
             (u"，而兩者", u"；兩者"), (u"，而其", u"；其"), (u"，而本文", u"；本文"),
             (u"，而標準", u"；標準"), (u"，而不是", u"，不是")]:
    s = s.replace(a, b)

io.open(P, "w", encoding="utf-8").write(s)
N1 = {w: s.count(w) for w in N0}
for w in N0: print(u"  %-6s %4d -> %4d" % (w, N0[w], N1[w]))
print(u"\n已改寫 %s" % P)
