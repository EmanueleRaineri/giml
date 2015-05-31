PLACE=home

ifeq ($(PLACE),work)
BED=/usr/bin
else
BED=/Users/emanueleraineri/bin/bedtools2/bin
endif

DATA  =~/Desktop/meth_data
METH = $(DATA)/G199_cpg.chr1.txt.gz $(DATA)/G200_cpg.chr1.txt.gz $(DATA)/G201_cpg.chr1.txt.gz $(DATA)/G202_cpg.chr1.txt.gz
GIMLI1000 = $(DATA)/G199_cpg.chr1.gimli.1000 $(DATA)/G200_cpg.chr1.gimli.1000 $(DATA)/G201_cpg.chr1.gimli.1000 $(DATA)/G202_cpg.chr1.gimli.1000 
C004GDH1GIMLI=$(DATA)/C004GD51_cpg.chr1.gimli.gz
GIMLICMD = ./gimli - 0.1:0.2:0.5:1:2:5:10:20:50:100:200:500:1000:2000:5000

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

G199.chr1.correla.txt:
	zcat $(DATA)/G199_cpg.chr1.txt.gz | awk '{print $$2,$$4}' | ./correla 1 2>/dev/null > G199.chr1.correla.txt

G200.chr1.correla.txt:
	zcat $(DATA)/G200_cpg.chr1.txt.gz | awk '{print $$2,$$4}' | ./correla 1 2>/dev/null > G200.chr1.correla.txt

G201.chr1.correla.txt:
	zcat $(DATA)/G201_cpg.chr1.txt.gz | awk '{print $$2,$$4}' | ./correla 1 2>/dev/null > G201.chr1.correla.txt

G202.chr1.correla.txt:
	zcat $(DATA)/G202_cpg.chr1.txt.gz | awk '{print $$2,$$4}' | ./correla 1 2>/dev/null > G202.chr1.correla.txt

out.correla.eps: G199.chr1.correla.txt G200.chr1.correla.txt G201.chr1.correla.txt G202.chr1.correla.txt
	Rscript correla_fig.R

$(C004GDH1GIMLI): $(DATA)/C004GD51_cpg.chr1.txt.gz
	zcat $^ | awk '{print $$1,$$2,$$6,$$7}' | $(GIMLICMD) 2> /dev/null | gzip -c > $@	


$(DATA)/G199_cpg.chr1.gimli.gz: $(DATA)/G199_cpg.chr1.txt.gz
	zcat $(DATA)/G199_cpg.chr1.txt.gz | awk '{print $$1,$$2,$$6,$$7}' | $(GIMLICMD) 2> $(DATA)/G199_cpg.chr1.gimli.log | gzip -c > $(DATA)/G199_cpg.chr1.gimli.gz
	
$(DATA)/G202_cpg.chr1.gimli.gz: $(DATA)/G202_cpg.chr1.txt.gz
	zcat $(DATA)/G202_cpg.chr1.txt.gz | awk '{print $$1,$$2,$$6,$$7}' | $(GIMLICMD) 2> /dev/null | gzip -c > $(DATA)/G202_cpg.chr1.gimli.gz

$(DATA)/G200_cpg.chr1.gimli.gz : $(DATA)/G200_cpg.chr1.txt.gz
	zcat $(DATA)/G200_cpg.chr1.txt.gz | awk '{print $$1,$$2,$$6,$$7}' | $(GIMLICMD) 2> /dev/null | gzip -c > $(DATA)/G200_cpg.chr1.gimli.gz

$(DATA)/G201_cpg.chr1.gimli.gz : $(DATA)/G201_cpg.chr1.txt.gz	
	zcat $(DATA)/G201_cpg.chr1.txt.gz | awk '{print $$1,$$2,$$6,$$7}' | $(GIMLICMD) 2> /dev/null | gzip -c > $(DATA)/G201_cpg.chr1.gimli.gz

$(GIMLI1000) : $(DATA)/G199_cpg.chr1.gimli.gz $(DATA)/G200_cpg.chr1.gimli.gz $(DATA)/G201_cpg.chr1.gimli.gz $(DATA)/G202_cpg.chr1.gimli.gz
	zcat $(DATA)/G199_cpg.chr1.gimli.gz | awk '$$NF==1000' > $(DATA)/G199_cpg.chr1.gimli.1000
	zcat $(DATA)/G200_cpg.chr1.gimli.gz | awk '$$NF==1000' > $(DATA)/G200_cpg.chr1.gimli.1000
	zcat $(DATA)/G201_cpg.chr1.gimli.gz | awk '$$NF==1000' > $(DATA)/G201_cpg.chr1.gimli.1000
	zcat $(DATA)/G202_cpg.chr1.gimli.gz | awk '$$NF==1000' > $(DATA)/G202_cpg.chr1.gimli.1000

$(DATA)/G199_cpg.chr1.gimli.100: $(DATA)/G199_cpg.chr1.gimli.gz
	zcat $(DATA)/G199_cpg.chr1.gimli.gz | awk '$$NF==100' > $(DATA)/G199_cpg.chr1.gimli.100

$(DATA)/G202_cpg.chr1.gimli.100: $(DATA)/G202_cpg.chr1.gimli.gz
	zcat $(DATA)/G202_cpg.chr1.gimli.gz | awk '$$NF==100' > $(DATA)/G202_cpg.chr1.gimli.100

$(DATA)/G199_cpg.chr1.gimli.10: $(DATA)/G199_cpg.chr1.gimli.gz
	zcat $(DATA)/G199_cpg.chr1.gimli.gz | awk '$$NF==10' > $(DATA)/G199_cpg.chr1.gimli.10

$(DATA)/G202_cpg.chr1.gimli.10: $(DATA)/G202_cpg.chr1.gimli.gz
	zcat $(DATA)/G202_cpg.chr1.gimli.gz | awk '$$NF==10' > $(DATA)/G202_cpg.chr1.gimli.10

$(DATA)/G199_cpg.chr1.gimli.100.filtered: $(DATA)/G199_cpg.chr1.gimli.100
	awk '$$4/($$3-$$2+1)>0.0' $^ > $@ 

$(DATA)/G202_cpg.chr1.gimli.100.filtered: $(DATA)/G202_cpg.chr1.gimli.100
	awk '$$4/($$3-$$2+1)>0.0' $^ > $@ 


gimli1000: $(GIMLI1000)

#boxplots
fig_g199_pmd_chr10.eps: G199_cpg.chr10.slice.txt G199_cpg.chr10.slice.gimli
	Rscript plot_pmd.R $^ $@

fig_boxplot_size.eps fig_boxplot_lik.eps: make_boxplots.R
	Rscript make_boxplots.R


delta.vs.cov.txt: cov.effect.R
	Rscript cov.effect.R > delta.vs.cov.txt


delta.vs.cov.mean.txt : delta.vs.cov.txt
	cat $^  | datamash -g 1 mean 5 sstdev 5 > $@


fig_delta_cov.eps: delta.vs.cov.txt delta.vs.cov.R
	Rscript delta.vs.cov.R


#fig_variance random segments vs gimli segments

random_mean_var_le15.txt: $(DATA)/G199_cpg.chr1.txt.gz
	zcat $^ | ocaml str.cma random.cpgs.ml > $@


gimli_mean_var_le15.txt : $(DATA)/G199_cpg.chr1.gimli.gz
		zcat /home/emanuele/Desktop/meth_data/G199_cpg.chr1.gimli.gz | awk '$$4=15' | awk 'BEGIN{c=0}{if (rand()>0.1) {c=c+1; print $$0;}; if (c==10000) exit 0}' | awk '{print $$2"\t"$$3"\t"$$8"\t"$$9}' > $@

fig_variance.eps: random.vs.gimli.R
	Rscript random.vs.gimli.R


#fig_jumps
fig_jumps.eps: C000S5A1bs_cpg.chr1.1655618.1656083.txt C0010KA2bs_cpg.chr1.1655618.1656083.txt C001UYA3bs_cpg.chr1.1655618.1656083.txt C004SQ51_cpg.chr1.1655618.1656083.txt plot_jumps.R
	Rscript plot_jumps.R

#figrho
figrho.eps: plot_rho.R
	Rscript $^

#figures: out.correla.eps figrho.eps fig2.eps fig_variance.eps fig6.eps G199.G200.G201.G202.chr1.gimli.eps G199.G202.100.200.dmr.eps fig_boxplot_size.eps fig_boxplot_lik.eps boxplot_example_3.eps

figures: out.correla.eps figrho.eps fig_delta_cov.eps fig_variance.eps fig_jumps.eps fig_boxplot_size.eps fig_boxplot_lik.eps 

#G199.G200.G201.G202.chr1.gimli.eps: $(DATA)/G199_cpg.chr1.gimli.100.filtered $(DATA)/G202_cpg.chr1.gimli.100.filtered
#	Rscript figure2.R $(DATA)	

#boxplot1.eps boxplot2.eps : $(DATA)/C004GD51_cpg.chr1.gimli.gz
#	Rscript make_boxplot1.R $(DATA)

#G199.G202.10.200.dmr.eps : G199.G202.10.200.dmr example2.R
#	Rscript example2.R $(DATA) G199.G202.10.200.dmr

#G199.G202.100.200.dmr.eps : G199.G202.100.200.dmr example2.R
#	Rscript example2.R $(DATA) G199.G202.100.200.dmr

#G199.G202.1000.200.dmr.eps : G199.G202.1000.200.dmr example2.R
#	Rscript example2.R $(DATA) G199.G202.1000.200.dmr

#boxplot_example_3.eps : rpkm.vs.met
	Rscript boxplot_example_3.R	

#example2

#G199.G202.10.200.dmr : $(DATA)/G199_cpg.chr1.gimli.10  $(DATA)/G202_cpg.chr1.gimli.10
#	$(BED)/bedtools intersect -a $(DATA)/G199_cpg.chr1.gimli.10 -b $(DATA)/G202_cpg.chr1.gimli.10 -wao | awk '$$NF>=200 && ($$7<$$15 || $$17<$$5)' > $@


#G199.G202.100.200.dmr : $(DATA)/G199_cpg.chr1.gimli.100  $(DATA)/G202_cpg.chr1.gimli.100
#	$(BED)/bedtools intersect -a $(DATA)/G199_cpg.chr1.gimli.100 -b $(DATA)/G202_cpg.chr1.gimli.100 -wao | awk '$$NF>=200 && ($$7<$$15 || $$17<$$5)' > $@

#G199.G202.1000.200.dmr : $(DATA)/G199_cpg.chr1.gimli.1000  $(DATA)/G202_cpg.chr1.gimli.1000
#	$(BED)/bedtools intersect -a $(DATA)/G199_cpg.chr1.gimli.1000 -b $(DATA)/G202_cpg.chr1.gimli.1000 -wao | awk '$$NF>=200 && ($$7<$$15 || $$17<$$5)' > $@

###example 3###

#S000RD13.tss.bed : S000RD13.gene_quantification.gem_grape_crg.20130415.gff
#	cat $^ | awk 'BEGIN{FS="\t";OFS="\t"}{if ($$7=="+") {print $$1,$$4,$$4,$$7} else {print $$1,$$5,$$5,$$7}}' | awk '$$1=="chr1"'  > $@ 

#$(DATA)/C004GD51_cpg.chr1.gimli.tss.bed : $(DATA)/C004GD51_cpg.chr1.gimli.gz S000RD13.tss.bed
#	$(BED)/bedtools intersect -b $(DATA)/C004GD51_cpg.chr1.gimli.gz -a S000RD13.tss.bed -wao > $@ 

#C004GD12.rpkm: C004GD12.gene_quantification.gem_grape_crg.20130415.gff 
#	awk 'BEGIN{FS=";"}{print $$1,$$3}' $^ | awk '$$1=="chr1"' | awk 'BEGIN{FS="\t";OFS="\t"}{print $$1,$$4,$$5,$$7,$$9}' | sed -e "s/gene_id//g " | sed -e "s/\"//g" | sed -e "s/RPKM//g" | tr -s ' ' '\t' > $@

#C004GD51_cpg.chr1.gimli.tss.filtered.bed : filter.tss.gimli.py $(DATA)/C004GD51_cpg.chr1.gimli.tss.bed
#	python filter.tss.gimli.py $(DATA) > $@

#rpkm.vs.met : join.tss.rpkm.py C004GD51_cpg.chr1.gimli.tss.filtered.bed C004GD12.rpkm
#	python join.tss.rpkm.py > $@ 

############################

gimli_paper.dvi: gimli_paper.tex gimli_paper.bib figures
	latex gimli_paper.tex
	bibtex gimli_paper
	latex gimli_paper.tex
	latex gimli_paper.tex

gimli_paper.pdf: gimli_paper.dvi
	dvipdf $^

##################

.PHONY : clean 

clean:
	rm -f gimli gimli_profile gimli_static gimli_optimized out.gimli.2 correla
	rm -f gimli_paper.dvi gimli_paper.pdf out.correla.eps boxplot1.eps boxplot2.eps boxplot_example_3.eps
	rm -f G199.G200.G201.G202.chr1.gimli.eps G199.G202.100.200.dmr
	rm -f $(DATA)/G199_cpg.chr1.gimli.100 $(DATA)/G202_cpg.chr1.gimli.100
	rm -f $(DATA)/G199_cpg.chr1.gimli.gz $(DATA)/G202_cpg.chr1.gimli.gz 
