# 以 xelatex 編譯中文報告與簡報
# TinyTeX 裝在 ~/Library/TinyTeX，不在 PATH 上時用下面這行補
TEXBIN  := $(HOME)/Library/TinyTeX/bin/universal-darwin
XELATEX := PATH="$(TEXBIN):$$PATH" xelatex -interaction=nonstopmode
BIBTEX  := PATH="$(TEXBIN):$$PATH" bibtex

DOCS := 進度報告 進度報告_v2 進度報告_0929 進度報告_1006 進度報告_簡版 進度報告_簡報
PDFS := $(addsuffix .pdf,$(DOCS))

.PHONY: all full v2 v3 v4 v5 part1 brief slides clean distclean

all: $(PDFS)

full:   進度報告.pdf
v2:     進度報告_v2.pdf
v3:     進度報告_0929.pdf
v4:     進度報告_1006.pdf
v5:     進度報告_1005_初版.pdf
part1:  進度報告_1005_深入研究_Part1.pdf
brief:  進度報告_簡版.pdf
slides: 進度報告_簡報.pdf

# 1005 初版的書目由 references.bib 經 bibtex 產生，故須四道
進度報告_1005_初版.pdf: 進度報告_1005_初版.tex references.bib
	$(XELATEX) 進度報告_1005_初版.tex
	$(BIBTEX)  進度報告_1005_初版
	$(XELATEX) 進度報告_1005_初版.tex
	$(XELATEX) 進度報告_1005_初版.tex

進度報告_1005_深入研究_Part1.pdf: 進度報告_1005_深入研究_Part1.tex references.bib
	$(XELATEX) 進度報告_1005_深入研究_Part1.tex
	$(BIBTEX)  進度報告_1005_深入研究_Part1
	$(XELATEX) 進度報告_1005_深入研究_Part1.tex
	$(XELATEX) 進度報告_1005_深入研究_Part1.tex

# 其餘各版無 bibtex，跑兩次以解析交互參照
%.pdf: %.tex
	$(XELATEX) $<
	$(XELATEX) $<

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb *.bbl *.blg

distclean: clean
	rm -f $(PDFS)
