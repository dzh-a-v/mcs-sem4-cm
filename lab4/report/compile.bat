latexmk -xelatex -outdir=out report.tex
if exist .\report.pdf del .\report.pdf
move /Y .\out\report.pdf .\report.pdf
