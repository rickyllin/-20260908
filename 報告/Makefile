# 以 xelatex 編譯中文報告。TinyTeX 裝在 ~/Library/TinyTeX，不在 PATH 上時用下面這行補。
TEXBIN  := $(HOME)/Library/TinyTeX/bin/universal-darwin
XELATEX := PATH="$(TEXBIN):$$PATH" xelatex -interaction=nonstopmode
BIBTEX  := PATH="$(TEXBIN):$$PATH" bibtex

# 目前維護中的三份文件。
#   v1012   進度報告_1012_V1       完整版（由 tools/make_1012_v1.py 自 V2 產生）
#   brief   進度報告_1012_摘要版   十三頁摘要（問題、解決方案、策略）
#   v2      進度報告_1005_V2       1005 的完整版，為 1012 V1 的母本，只編輯這一份
#
# 已被取代者置於 封存/：1005 的 V3 與無精算概念版（兩者皆由 V2 衍生，
# 其角色現由 1012 摘要版取代），以及更早的各版。
# 理論初探兩份已移至 ../研究筆記/理論初探/，自備 Makefile。
.PHONY: all v1012 brief v2 gen clean distclean

all: v1012 brief v2

v1012:  進度報告_1012_V1.pdf
brief:  進度報告_1012_摘要版.pdf
v2:     進度報告_1005_V2.pdf

# 由母本重新產生 1012 V1 的 .tex（改完 V2 之後先跑這個）
gen:
	cd .. && python3 報告/tools/make_1012_v1.py

# 書目由 references.bib 經 bibtex 產生，故每份都要四道
%.pdf: %.tex references.bib
	$(XELATEX) $<
	$(BIBTEX)  $(basename $<)
	$(XELATEX) $<
	$(XELATEX) $<

clean:
	rm -f *.aux *.log *.out *.toc *.bbl *.blg *.nav *.snm *.vrb *.synctex.gz

distclean: clean
	rm -f 進度報告_1012_V1.pdf 進度報告_1012_摘要版.pdf 進度報告_1005_V2.pdf
