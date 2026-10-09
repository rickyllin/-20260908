# 理論初探兩份的編譯。書目同時載入 ../../報告/references.bib 與本地的增補檔，
# 故 BIBINPUTS 須包含報告目錄。
TEXBIN  := $(HOME)/Library/TinyTeX/bin/universal-darwin
XELATEX := PATH="$(TEXBIN):$$PATH" xelatex -interaction=nonstopmode
BIBTEX  := PATH="$(TEXBIN):$$PATH" BIBINPUTS=".:../../報告:" bibtex

.PHONY: all algebra quasi clean
all: algebra quasi
algebra: 理論初探_修勻加權與中心化的代數關係.pdf
quasi:   理論初探_擬概似與變異安定轉換.pdf

%.pdf: %.tex
	$(XELATEX) $<
	$(BIBTEX)  $(basename $<)
	$(XELATEX) $<
	$(XELATEX) $<

clean:
	rm -f *.aux *.log *.out *.toc *.bbl *.blg
