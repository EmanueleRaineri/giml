DDIR:=~/Desktop/meth_data
BED1:=$(DDIR)/fake1.bed
BED500:=$(DDIR)/fake2.500.bed
BED5000:=$(DDIR)/fake2.5000.bed
BED50000:=$(DDIR)/fake2.50000.bed
BED500000:=$(DDIR)/fake2.500000.bed
STRIPPED500:=$(DDIR)/fake2.500_cpg.stripped.txt.gz
STRIPPED5000:=$(DDIR)/fake2.5000_cpg.stripped.txt.gz
STRIPPED50000:=$(DDIR)/fake2.50000_cpg.stripped.txt.gz
STRIPPED500000:=$(DDIR)/fake2.500000_cpg.stripped.txt.gz

all : $(DDIR)/fake1_cpg.stripped.txt.gz\
	$(STRIPPED500)\
	$(STRIPPED5000)\
	$(STRIPPED50000)\
	$(STRIPPED500000)

$(BED1):
	echo "chr1\t1\t250000000\t0.9\t30" > $@

$(BED500): make_fake2_bed.py
	python $< $(BED500) 500 

$(BED5000): make_fake2_bed.py
	python $< $(BED5000) 5000

$(BED50000): make_fake2_bed.py
	python $< $(BED50000) 50000

$(BED500000): make_fake2_bed.py
	python $< $(BED500000) 500000

$(DDIR)/fake1_cpg.stripped.txt.gz : $(DDIR)/C000S5A1bs_cpg.txt.gz $(DDIR)/fake1.bed
	zcat $<  | awk '$$1=="chr1"' | awk '{print $$2}' \
	| python unroll_meth_bed.py $(DDIR)/fake1.bed \
	| python sim_meth.py | gzip -c >\
	$@

$(STRIPPED500): $(DDIR)/C000S5A1bs_cpg.txt.gz $(DDIR)/fake2.500.bed
	zcat $<  | awk '$$1=="chr1"' | awk '{print $$2}' \
	| python unroll_meth_bed.py $(DDIR)/fake2.500.bed \
	| python sim_meth.py | gzip -c >\
	$@

$(STRIPPED5000): $(DDIR)/C000S5A1bs_cpg.txt.gz $(DDIR)/fake2.5000.bed
	zcat $<  | awk '$$1=="chr1"' | awk '{print $$2}' \
	| python unroll_meth_bed.py $(DDIR)/fake2.5000.bed \
	| python sim_meth.py | gzip -c >\
	$@

$(STRIPPED50000): $(DDIR)/C000S5A1bs_cpg.txt.gz $(DDIR)/fake2.50000.bed
	zcat $<  | awk '$$1=="chr1"' | awk '{print $$2}' \
	| python unroll_meth_bed.py $(DDIR)/fake2.50000.bed \
	| python sim_meth.py | gzip -c >\
	$@

$(STRIPPED500000): $(DDIR)/C000S5A1bs_cpg.txt.gz $(DDIR)/fake2.500000.bed
	zcat $<  | awk '$$1=="chr1"' | awk '{print $$2}' \
	| python unroll_meth_bed.py $(DDIR)/fake2.500000.bed \
	| python sim_meth.py | gzip -c >\
	$@



.PHONY: clean

clean: 
	rm -f $(DDIR)/fake*.bed $(DDIR)/fake*_cpg.stripped.txt.gz
