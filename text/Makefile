.PHONY: bibtex compile1 compile2 compile3 open clean s ss log fix

PDFFILE=main.pdf
MAINFILE=main
SECTIONS_FOLDER=sections

all: clean compile1 bibtex compile2 compile3 log
s: compile1 log # simple
ss: compile1 compile2 log

bibtex:
	bibtex ${MAINFILE}.aux

compile1 compile2 compile3:
	pdflatex ${MAINFILE}.tex

open:
	open ${PDFFILE}

fix:
	python tables/fix.py

clean:
	rm -f *.aux *.bbl *.bcf *.blg *.fdb_latexmk *.fls *.lof *.log *.lot *.out *.run.xml *.tex.bbl *.tex.blg *.toc *.tdo

log:
	echo "`date` - \
	Pages: `pdfinfo ${PDFFILE} | grep Pages | sed 's/[^0-9]*//'`\
	, Words: `texcount -total sections/*.tex | grep Words\ in\ text | sed 's/[^0-9]*//'`" >> log.txt