# 以 xelatex 編譯中文報告。TinyTeX 裝在 ~/Library/TinyTeX，不在 PATH 上時用下面這行補。
TEXBIN  := $(HOME)/Library/TinyTeX/bin/universal-darwin
XELATEX := PATH="$(TEXBIN):$$PATH" xelatex -interaction=nonstopmode
BIBTEX  := PATH="$(TEXBIN):$$PATH" bibtex

# 目前維護中的三份文件。歷次已被取代的版本置於 封存/，不再編譯。
#   v2      進度報告_1005_V2             完整版（母本，只編輯這一份）
#   v3      進度報告_1005_V3             精簡版（主文四十頁以內，由 tools/make_v3.py 產生）
#   nonact  進度報告_1005_V2_無精算概念   統計版（由 tools/make_statistics_version.py 產生）
#   theory  進度報告_1012_理論初探         理論初探（獨立撰寫，書目另加 references_1012.bib）
#   quasi   進度報告_理論初探_擬概似…      理論初探（不取對數，書目另加 references_1013.bib）
.PHONY: all v2 v3 nonact theory quasi gen clean distclean

all: v2 v3 nonact theory quasi

v2:     進度報告_1005_V2.pdf
v3:     進度報告_1005_V3.pdf
nonact: 進度報告_1005_V2_無精算概念.pdf
theory: 進度報告_1012_理論初探.pdf
quasi:  進度報告_理論初探_擬概似與變異安定轉換.pdf

# 理論初探另載入書目增補檔，故不走下方的通用規則
進度報告_1012_理論初探.pdf: 進度報告_1012_理論初探.tex references.bib references_1012.bib
	$(XELATEX) $<
	$(BIBTEX)  $(basename $<)
	$(XELATEX) $<
	$(XELATEX) $<

進度報告_理論初探_擬概似與變異安定轉換.pdf: 進度報告_理論初探_擬概似與變異安定轉換.tex references.bib references_1013.bib
	$(XELATEX) $<
	$(BIBTEX)  $(basename $<)
	$(XELATEX) $<
	$(XELATEX) $<

# 由母本重新產生兩個衍生版的 .tex（改完 V2 之後先跑這個）
gen:
	cd .. && python3 報告/tools/make_v3.py
	cd .. && python3 報告/tools/make_statistics_version.py

# 書目由 references.bib 經 bibtex 產生，故每份都要四道
%.pdf: %.tex references.bib
	$(XELATEX) $<
	$(BIBTEX)  $(basename $<)
	$(XELATEX) $<
	$(XELATEX) $<

clean:
	rm -f *.aux *.log *.out *.toc *.bbl *.blg *.nav *.snm *.vrb *.synctex.gz

distclean: clean
	rm -f 進度報告_1005_V2.pdf 進度報告_1005_V3.pdf 進度報告_1005_V2_無精算概念.pdf \
	      進度報告_1012_理論初探.pdf 進度報告_理論初探_擬概似與變異安定轉換.pdf
