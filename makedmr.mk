SHELL=/bin/bash
DDIR:=~/Desktop/meth_data
NAME1:=C001UYA3bs
NAME2:=C000S5A1bs
CPG1:=$(DDIR)/$(NAME1)_cpg.txt.gz
CPG2:=$(DDIR)/$(NAME2)_cpg.txt.gz
STRIPPED1:=$(DDIR)/$(NAME1)_cpg.stripped.txt.gz
STRIPPED2:=$(DDIR)/$(NAME2)_cpg.stripped.txt.gz
GIMLI1:=$(DDIR)/$(NAME1)_cpg.gimli
GIMLI2:=$(DDIR)/$(NAME2)_cpg.gimli
INTERSECT1000:=$(DDIR)/$(NAME1).$(NAME2).gimli.1000.intersect
INTERSECT100:=$(DDIR)/$(NAME1).$(NAME2).gimli.100.intersect
INTERSECT10:=$(DDIR)/$(NAME1).$(NAME2).gimli.10.intersect
WIN:=../../cpg_pipeline/gimli.windows

.SECONDARY:



dmrs: $(DDIR)/$(NAME1).$(NAME2).gimli.1000.dmr $(DDIR)/$(NAME1).$(NAME2).gimli.100.dmr

$(STRIPPED1): $(CPG1)
	zcat $(CPG1) | \
	awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}' |\
	gzip -c > $@

$(STRIPPED2): $(CPG2)
	zcat $(CPG2) | \
	awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}' |\
	gzip -c > $@

$(GIMLI1): $(STRIPPED1)
	zcat $(STRIPPED1) | \
	./gimli - 1:10:100:1000 > $@ 2>$(GIMLI1).log

$(GIMLI2): $(STRIPPED2)
	zcat $(STRIPPED2) | \
	./gimli - 1:10:100:1000 > $@ 2>$(GIMLI2).log

$(INTERSECT1000): $(GIMLI1) $(GIMLI2)
	bedtools intersect \
	-a <(awk '{if ($$NF==1000) print $$1"\t"$$2"\t"$$3}' $(GIMLI1)) \
	-b <(awk '{if ($$NF==1000) print $$1"\t"$$2"\t"$$3}' $(GIMLI2)) -wao \
	| awk '{lb=($$2>$$5)?$$2:$$5;ub=($$3>$$6)?$$6:$$3;print $$0"\t"lb"\t"ub}' \
	> $@ 2>/dev/null

$(INTERSECT100): $(GIMLI1) $(GIMLI2)
	bedtools intersect \
	-a <(awk '{if ($$NF==100) print $$1"\t"$$2"\t"$$3}' $(GIMLI1)) \
	-b <(awk '{if ($$NF==100) print $$1"\t"$$2"\t"$$3}' $(GIMLI2)) -wao \
	| awk '{lb=($$2>$$5)?$$2:$$5;ub=($$3>$$6)?$$6:$$3;print $$0"\t"lb"\t"ub}' \
	> $@ 2>/dev/null

$(INTERSECT10): $(GIMLI1) $(GIMLI2)
	bedtools intersect \
	-a <(awk '{if ($$NF==10) print $$1"\t"$$2"\t"$$3}' $(GIMLI1)) \
	-b <(awk '{if ($$NF==10) print $$1"\t"$$2"\t"$$3}' $(GIMLI2)) -wao \
	| awk '{lb=($$2>$$5)?$$2:$$5;ub=($$3>$$6)?$$6:$$3;print $$0"\t"lb"\t"ub}' \
	> $@ 2>/dev/null

$(DDIR)/$(NAME1).$(NAME2).gimli.100.counts.1 :  $(STRIPPED1) $(INTERSECT100)
	zcat $< | \
	$(WIN) <(awk '{print $$1"\t"$$(NF-1)"\t"$$NF}' $(INTERSECT100)) \
	> $@

$(DDIR)/$(NAME1).$(NAME2).gimli.100.counts.2 :  $(STRIPPED2) $(INTERSECT100)
	zcat $< |  \
	$(WIN) <(awk '{print $$1"\t"$$(NF-1)"\t"$$NF}' $(INTERSECT100)) \
	> $@

$(DDIR)/$(NAME1).$(NAME2).gimli.1000.counts.1 :  $(STRIPPED1) $(INTERSECT1000)
	zcat $< | \
	$(WIN) <(awk '{print $$1"\t"$$(NF-1)"\t"$$NF}' $(INTERSECT1000)) \
	> $@

$(DDIR)/$(NAME1).$(NAME2).gimli.1000.counts.2 :  $(STRIPPED2) $(INTERSECT1000)
	zcat $< |  \
	$(WIN) <(awk '{print $$1"\t"$$(NF-1)"\t"$$NF}' $(INTERSECT1000)) \
	> $@

%.counts: %.counts.1 %.counts.2
	paste $^ | grep -v "-" |\
	awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$3,$$6,$$7,$$17,$$18}' > $@

%.lrt: %.counts
	Rscript likratiotest.R $< > $@


%.dmr: %.lrt
	awk  '$$NF<0.01' $< > $@

.PHONY: clean

clean:
	rm -f $(INTERSECT100) $(DDIR)/$(NAME1).$(NAME2).gimli.*.counts.1 
	rm -f $(INTERSECT1000) $(DDIR)/$(NAME1).$(NAME2).gimli.*.counts.2
	rm -f $(DDIR)/$(NAME1).$(NAME2).gimli.*.lrt 
	rm -f $(GIMLI1) $(GIMLI2) $(DDIR)/$(NAME1).$(NAME2).gimli.*.counts
	rm -f $(STRIPPED1) $(STRIPPED2)
	rm -f $(DDIR)/$(NAME1).$(NAME2).gimli.*.dmr
	rm -f $(GIMLI1).log $(GIMLI2).log
