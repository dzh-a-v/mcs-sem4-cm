latexmk -xelatex -outdir=out report.tex
rd .\report.pdf
move .\out\report.pdf .\report.pdf