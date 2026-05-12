latexmk -xelatex -outdir=out lab5_report.tex
if exist .\lab5_report.pdf del .\lab5_report.pdf
move /Y .\out\lab5_report.pdf .\lab5_report.pdf
