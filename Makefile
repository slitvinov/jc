M=main
# arguments
a=-interaction=nonstopmode
# tex
t=sty/def.tex  sty/notheme.tex t/title.tex
i=

$M.pdf: $M.tex $t $i
	pdflatex $a $< && \
	pdflatex $a $<

.PHONY: clean
clean:
	@echo clean
	@rm -f $M.aux $M.log $M.nav $M.out $M.pdf $M.snm $M.toc
