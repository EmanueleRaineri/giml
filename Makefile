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


correla: correla.c
	gcc -o correla correla.c -lm


G199.G202.20.200.dmr.eps : example2.R
	Rscript example2.R


G199.chr1.correla.txt:
	zcat ~/Desktop/meth_data/G199_cpg.chr1.txt.gz | awk '{print $2,$4}' | ./correla 1 2>/dev/null > G199.chr1.correla.txt

G200.chr1.correla.txt:
	zcat ~/Desktop/meth_data/G200_cpg.chr1.txt.gz | awk '{print $2,$4}' | ./correla 1 2>/dev/null > G200.chr1.correla.txt

G201.chr1.correla.txt:
	zcat ~/Desktop/meth_data/G201_cpg.chr1.txt.gz | awk '{print $2,$4}' | ./correla 1 2>/dev/null > G201.chr1.correla.txt

G202.chr1.correla.txt:
	zcat ~/Desktop/meth_data/G202_cpg.chr1.txt.gz | awk '{print $2,$4}' | ./correla 1 2>/dev/null > G202.chr1.correla.txt

out.correla.eps: G199.chr1.correla.txt G200.chr1.correla.txt G201.chr1.correla.txt G202.chr1.correla.txt
	Rscript correla_fig.R

figures: out.correla.eps G199.G202.chr1.gimli.eps G199.G202.20.200.dmr.eps  



gimli_paper.dvi: gimli_paper.tex gimli_paper.bib figures
	latex gimli_paper.tex
	bibtex gimli_paper
	latex gimli_paper.tex
	latex gimli_paper.tex

gimli_paper.pdf: gimli_paper.dvi
	dvipdf $^

clean:
	rm -f gimli gimli_profile gimli_static gimli_optimized out.gimli.2 correla
	rm -f gimli_paper.dvi gimli_paper.pdf out.correla.eps
