DDIR:=~/Desktop/meth_data
BED1:=$(DDIR)/fake1.bed
BED2:=$(DDIR)/fake2.bed
STRIPPED:=$(DDIR)/fake%_cpg.stripped.txt.gz


all : $(DDIR)/fake1_cpg.stripped.txt.gz $(DDIR)/fake2_cpg.stripped.txt.gz

$(BED1) $(BED2): make_fakes_bed.py
	python $< $(BED1) $(BED2) 

$(STRIPPED): $(DDIR)/C000S5A1bs_cpg.txt.gz $(DDIR)/fake%.bed
	zcat $<  | awk '$$1=="chr1"' | awk '{print $$2}' \
	| python unroll_meth_bed.py $(DDIR)/fake$*.bed \
	| python sim_meth.py | gzip -c >\
	$@

.PHONY: clean

clean: 
	rm -f fake?.bed fake?_cpg.stripped.txt.gz
