# -*- coding: utf-8 -*-
"""把因果連接詞的密度降下來（就地改寫 V2，只跑一次）。
五篇參考論文的「故／因而／因此」合計約每千字 0.5 至 1.0；把「故」一律改成
「因此」之後，本文的密度升到 3.8，必須回頭刪減。規則有二：
前文已有蓋因／由於／既然／因為／源自等因果標記者，刪去重複的連接詞；
其餘每三處取一處改用「於是」，以免單一詞彙過度重複。"""
import io, re
P = "報告/進度報告_1005_V2.tex"
s = io.open(P, encoding="utf-8").read()
n0 = s.count(u"因此")
MARK = re.compile(u"蓋因|由於|既然|因為|源自|之所以|理由|原因")
res = []; idx = 0; k = 0; dropped = 0; swapped = 0
while True:
    j = s.find(u"因此", idx)
    if j < 0:
        res.append(s[idx:]); break
    res.append(s[idx:j])
    ctx = s[max(0, j - 70):j]
    if MARK.search(ctx):
        dropped += 1              # 前文已有因果標記，刪去重複的連接詞
    else:
        k += 1
        if k % 3 == 0:
            res.append(u"於是"); swapped += 1
        else:
            res.append(u"因此")
    idx = j + 2
s = u"".join(res)
# 刪去連接詞後可能留下的「，，」或行首逗號
s = s.replace(u"，，", u"，").replace(u"\n，", u"\n")
io.open(P, "w", encoding="utf-8").write(s)
print(u"  因此 %d -> %d（刪去 %d、改為於是 %d）" % (n0, s.count(u"因此"), dropped, swapped))
for w in [u"故", u"因而", u"因此", u"於是"]:
    print(u"  %-4s %d" % (w, s.count(w)))
