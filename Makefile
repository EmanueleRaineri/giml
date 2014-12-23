gimli: greedy.c
	gcc -Wall  -o $@ greedy.c -lm

gimli_optimized : greedy.c
	gcc -Wall  -o $@ greedy.c -lm -O3
	

gimli_profile : greedy.c
	gcc -Wall  -o $@ greedy.c -lm -g -pg


gimli_static: greedy.c
	gcc -Wall  -o $@ greedy.c -lm -static	

out.gimli.2: gimli 
	./gimli G199.sample > out.gimli
	awk '$$NF==2' out.gimli > out.gimli.2

test: out.gimli.2.ref out.gimli.2 
	diff $^


G199.G202.20.200.dmr.eps : example2.R
	Rscript example2.R

figures: out.correla.eps G199.G202.chr1.gimli.eps G199.G202.20.200.dmr.eps  



gimli_paper.dvi: gimli_paper.tex gimli_paper.bib figures
	latex gimli_paper.tex
	bibtex gimli_paper
	latex gimli_paper.tex
	latex gimli_paper.tex

gimli_paper.pdf: gimli_paper.dvi
	dvipdf $^

clean:
	rm -f gimli gimli_profile gimli_static gimli_optimized out.gimli.2

