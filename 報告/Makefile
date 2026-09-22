# 以 xelatex 編譯中文報告與簡報
# TinyTeX 裝在 ~/Library/TinyTeX，不在 PATH 上時用下面這行補
TEXBIN  := $(HOME)/Library/TinyTeX/bin/universal-darwin
XELATEX := PATH="$(TEXBIN):$$PATH" xelatex -interaction=nonstopmode

DOCS := 進度報告 進度報告_v2 進度報告_0929 進度報告_簡版 進度報告_簡報
PDFS := $(addsuffix .pdf,$(DOCS))

.PHONY: all full v2 v3 brief slides clean distclean

all: $(PDFS)

full:   進度報告.pdf
v2:     進度報告_v2.pdf
v3:     進度報告_0929.pdf
brief:  進度報告_簡版.pdf
slides: 進度報告_簡報.pdf

# 跑兩次以解析交互參照
%.pdf: %.tex
	$(XELATEX) $<
	$(XELATEX) $<

clean:
	rm -f *.aux *.log *.nav *.out *.snm *.toc *.vrb

distclean: clean
	rm -f $(PDFS)
