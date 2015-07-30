PLACE=home

ifeq ($(PLACE),work)
BED=/usr/bin
else
BED=/Users/emanueleraineri/bin/bedtools2/bin
endif

DATA  = ~/Desktop/meth_data
METH = $(DATA)/G199_cpg.chr1.txt.gz 
GIMLICMD = ./gimli - 1:10:100:1000

lik_of_counts: lik_of_counts.ml
	ocamlopt.opt -ccopt -static -o $@ str.cmxa $<

gimli: gimli.c
	gcc -Wall  -o $@ $^ -lm

gimli_optimized : gimli.c
	gcc -Wall  -o $@ $^ -lm -O3

gimli_profile : gimli.c
	gcc -Wall  -o $@ $^ -lm -g -pg

gimli_static: greedy.c
	gcc -Wall  -o $@ $^ -lm -static	

out.gimli.2: gimli 
	./gimli G199.sample > out.gimli
	awk '$$NF==2' out.gimli > out.gimli.2

test: out.gimli.2.ref out.gimli.2 
	diff $^

.PHONY : clean 

clean:
	rm -f gimli gimli_profile gimli_static gimli_optimized out.gimli.2 correla
	rm -f gimli_paper.dvi gimli_paper.pdf out.correla.eps boxplot1.eps boxplot2.eps boxplot_example_3.eps
	rm -f G199.G200.G201.G202.chr1.gimli.eps G199.G202.100.200.dmr
