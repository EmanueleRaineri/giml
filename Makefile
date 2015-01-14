DATA  =~/Desktop/meth_data
METH = $(DATA)/G199_cpg.chr1.txt.gz $(DATA)/G200_cpg.chr1.txt.gz $(DATA)/G201_cpg.chr1.txt.gz $(DATA)/G202_cpg.chr1.txt.gz
GIMLI1000 = $(DATA)/G199_cpg.chr1.gimli.1000 $(DATA)/G200_cpg.chr1.gimli.1000 $(DATA)/G201_cpg.chr1.gimli.1000 $(DATA)/G202_cpg.chr1.gimli.1000 
BED=/Users/emanueleraineri/bedtools2/bin
C004GDH1GIMLI:$(DATA)/C004GD51_cpg.chr1.gimli.gz

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
	zcat $(DATA)/G199_cpg.chr1.txt.gz | awk '{print $2,$4}' | ./correla 1 2>/dev/null > G199.chr1.correla.txt

G200.chr1.correla.txt:
	zcat $(DATA)/G200_cpg.chr1.txt.gz | awk '{print $2,$4}' | ./correla 1 2>/dev/null > G200.chr1.correla.txt

G201.chr1.correla.txt:
	zcat $(DATA)/G201_cpg.chr1.txt.gz | awk '{print $2,$4}' | ./correla 1 2>/dev/null > G201.chr1.correla.txt

G202.chr1.correla.txt:
	zcat $(DATA)/G202_cpg.chr1.txt.gz | awk '{print $2,$4}' | ./correla 1 2>/dev/null > G202.chr1.correla.txt

out.correla.eps: G199.chr1.correla.txt G200.chr1.correla.txt G201.chr1.correla.txt G202.chr1.correla.txt
	Rscript correla_fig.R

G199.G200.G201.G202.chr1.gimli.eps: $(GIMLI1000)
	Rscript figure2.R	

gencode_promoters: gencode.v19.TSS.notlow.chr1.promoters

gencode.v19.TSS.notlow.chr1.promoters: $(DATA)/gencode.v19.TSS.notlow.chr1.gff
	awk 'BEGIN{FS="\t";OFS="\t"}{if ($$7=="+") {print $$1,$$4,$$4+1000,$$7}; if ($$7=="-") { print $$1,$$4-1000,$$4,$$7} if ($$7!="+" && $$7!="-") {print "NaN"} }' $^ > $@

gimli1000: $(GIMLI1000)

$(GIMLI1000) : $(METH) 
	zcat $(DATA)/G199_cpg.chr1.txt.gz | awk '{print $$1,$$2,$$6,$$7}' | ./gimli 2> /dev/null | gzip -c > $(DATA)/G199_cpg.chr1.gimli.gz
	zcat $(DATA)/G200_cpg.chr1.txt.gz | awk '{print $$1,$$2,$$6,$$7}' | ./gimli 2> /dev/null | gzip -c > $(DATA)/G200_cpg.chr1.gimli.gz
	zcat $(DATA)/G201_cpg.chr1.txt.gz | awk '{print $$1,$$2,$$6,$$7}' | ./gimli 2> /dev/null | gzip -c > $(DATA)/G201_cpg.chr1.gimli.gz
	zcat $(DATA)/G202_cpg.chr1.txt.gz | awk '{print $$1,$$2,$$6,$$7}' | ./gimli 2> /dev/null | gzip -c > $(DATA)/G202_cpg.chr1.gimli.gz
	zcat ~/Desktop/meth_data/G199_cpg.chr1.gimli.gz | awk '$$NF==1000' > $(DATA)/G199_cpg.chr1.gimli.1000
	zcat ~/Desktop/meth_data/G200_cpg.chr1.gimli.gz | awk '$$NF==1000' > $(DATA)/G200_cpg.chr1.gimli.1000
	zcat ~/Desktop/meth_data/G201_cpg.chr1.gimli.gz | awk '$$NF==1000' > $(DATA)/G201_cpg.chr1.gimli.1000
	zcat ~/Desktop/meth_data/G202_cpg.chr1.gimli.gz | awk '$$NF==1000' > $(DATA)/G202_cpg.chr1.gimli.1000


boxplot1.eps boxplot2.eps : $(DATA)/C004GD51_cpg.chr1.gimli.gz
	Rscript make_boxplot1.R

figures: out.correla.eps G199.G200.G201.G202.chr1.gimli.eps G199.G202.20.200.dmr.eps boxplot1.eps boxplot2.eps

gimli_paper.dvi: gimli_paper.tex gimli_paper.bib figures
	latex gimli_paper.tex
	bibtex gimli_paper
	latex gimli_paper.tex
	latex gimli_paper.tex

gimli_paper.pdf: gimli_paper.dvi
	dvipdf $^

clean:
	rm -f gimli gimli_profile gimli_static gimli_optimized out.gimli.2 correla
	rm -f gimli_paper.dvi gimli_paper.pdf out.correla.eps boxplot1.eps boxplot2.eps

active_promoters: C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed

C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed : C004GDH1_12_Blueprint_release_082014_segments.chr1.bed
	awk '$$NF=="E5" || $$NF=="E6" || $$NF=="E7"' $^ > $@
active_promoters_intersect_gencode.txt : active_promoters gencode.v19.TSS.notlow.chr1.promoters
	$(BED)/bedtools intersect -a C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed -b gencode.v19.TSS.notlow.chr1.promoters > active_promoters_intersect_gencode.txt


C004GDH1.active.C004GDH1_cpg.chr1.gimli.1000 : C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed C004GD51_cpg.chr1.gimli.1000
	$(BED)/bedtools intersect -a C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed -b C004GD51_cpg.chr1.gimli.1000 -wao > C004GDH1.active.C004GDH1_cpg.chr1.gimli.1000	

C004GDH1.active.C004GDH1_cpg.chr1.gimli.100 : C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed C004GD51_cpg.chr1.gimli.100
	$(BED)/bedtools intersect -a C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed -b C004GD51_cpg.chr1.gimli.100 -wao > C004GDH1.active.C004GDH1_cpg.chr1.gimli.100	

#zcat ~/Desktop/meth_data/G199_cpg.chr1.gimli.gz | awk '$3-$2>100 && $4/($3-$2)>=0.05 && $5>=0.3 && $7<=0.7' | awk '{print $0,$3-$2}' > G199.pmd
#/Users/emanueleraineri/bedtools2/bin/bedtools intersect -a G199.pmd -b G200.pmd

gencode.chr1.genes.unique.tss :  gencode.v19.TSS.notlow.chr1.gff
	cat $^ | awk '{print $$10}' | sort | uniq -c | awk '{if ($$1==1) print $$2}' > $@

uniq.tss.coords: filter_gff.py gencode.chr1.genes.unique.tss gencode.v19.TSS.notlow.chr1.gff
	python filter_gff.py > uniq.tss.coords

#/Users/emanueleraineri/bedtools2/bin/bedtools intersect -a uniq.tss.coords -b $(C004GDH1GIMLI) -wao |  awk '$5!="." && $11-$9<0.3' | awk '{print $10}' | Rscript -e "summary(as.numeric(readLines(file('stdin'))))"

#/Users/emanueleraineri/bedtools2/bin/bedtools intersect -a uniq.tss.coords -b ~/Desktop/meth_data/C004GD51_cpg.chr1.gimli.gz  -wao |  awk '$5!="." && $11-$9<0.3' > uniq.tss.C004GD51.gimli

#awk '{print $0,$2$3}'  | awk '{if (a[$NF]==1) next; else {a[$NF]=1;print $0}}'

#/Users/emanueleraineri/bedtools2/bin/bedtools intersect -a C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed -b ~/Desktop/meth_data/C004GD51_cpg.chr1.gimli.gz -wao | awk '$5!="."' | awk '$7-$6>100 && $8/($7-$6)>0.05'  | awk '{print $0,$2$3}'  | awk '{if (a[$NF]==1) next; else {a[$NF]=1;print $0}}'

#/Users/emanueleraineri/bedtools2/bin/bedtools intersect -a C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed -b ~/Desktop/meth_data/C004GD51_cpg.chr1.gimli.gz -wao | awk '$5!="."' | awk '$8/($7-$6+1)>0.05'  | awk '{print $0,$2$3}'  | awk '{if (a[$NF]==1) next; else {a[$NF]=1;print $0}}' > active.vs.gimli.txt

#/Users/emanueleraineri/bedtools2/bin/bedtools intersect -a C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed -b ~/Desktop/meth_data/C004GD51_cpg.chr1.gimli.gz -wao | awk '$5!="."' | awk '$7-$6+1>100 && $8/($7-$6+1)>0.05'  | awk '{print $0,$2$3}'  | awk '{if (a[$NF]==1) next; else {a[$NF]=1;print $0}}' > active.vs.gimli.txt

uniq.tss.extended.100.coords : uniq.tss.coords
	awk 'BEGIN{OFS="\t"}{print $$1,$$2-100,$$3+100,$$4}' uniq.tss.coords > uniq.tss.extended.100.coords



S000RD13.gene.body.txt : S000RD13.gene_quantification.gem_grape_crg.20130415.gff
	cat $^ | awk 'BEGIN{FS="\t"}{print $$1,$$4,$$5,$$7}' | awk '$$1=="chr1"'  > $@ 

###example 3###

S000RD13.tss.bed : S000RD13.gene_quantification.gem_grape_crg.20130415.gff
	cat $^ | awk 'BEGIN{FS="\t";OFS="\t"}{if ($$7=="+") {print $$1,$$4,$$4,$$7} else {print $$1,$$5,$$5,$$7}}' | awk '$$1=="chr1"'  > $@ 


C004GD51_cpg.chr1.gimli.tss.bed : $(DATA)/C004GD51_cpg.chr1.gimli.gz S000RD13.tss.bed
	/Users/emanueleraineri/bedtools2/bin/bedtools intersect -b $(DATA)/C004GD51_cpg.chr1.gimli.gz -a S000RD13.tss.bed > $@ 


C004GD12.rpkm: C004GD12.gene_quantification.gem_grape_crg.20130415.gff 
	awk 'BEGIN{FS=";"}{print $$1,$$3}' $^ | awk '$$1=="chr1"' | awk 'BEGIN{FS="\t";OFS="\t"}{print $$1,$$4,$$5,$$7,$$9}' | sed -e "s/gene_id//g " | sed -e "s/\"//g" | sed -e "s/RPKM//g" | tr -s ' ' '\t' > $@

C004GD51_cpg.chr1.gimli.tss.filtered.bed : filter.tss.gimli.py C004GD51_cpg.chr1.gimli.tss.bed
	python filter.tss.gimli.py > $@

rpkm.vs.met : join.tss.rpkm.py C004GD51_cpg.chr1.gimli.tss.filtered.bed C004GD12.rpkm
	python join.tss.rpkm.py > $@ 



#S000RD13.promoters.1000.bed : S000RD13.gene_quantification.gem_grape_crg.20130415.gff
#	cat $^ | awk 'BEGIN{FS="\t";OFS="\t"} {if ($$7=="+") {print $$1,$$4-1000,$$4,$$7} else {print $$1,$$5,$$5+1000,$$7} }' | awk '$$1=="chr1"' > $@ 

