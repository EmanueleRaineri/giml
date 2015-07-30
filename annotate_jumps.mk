SHELL=/bin/bash
SCRATCH=/scratch/devel/eraineri
GTEST=$(SCRATCH)/gimli_test
ANNOTATION=/scratch/devel/heath/Blueprint/gencode.v15.annotation.gtf.gz
######example 1 : variable regions

######example 2 : jumps


#%_cpg.gimli.jumps: %_cpg.gimli
#awk 'BEGIN{c=1}{a[c]=$0;c++;if (c==4) {split(a[1],l1,"\t");split(a[2],l2,"\t");split(a[3],l3,"\t");if (l1[5]>5 && l2[5]>5 && l3[5]>5) {print a[1]"\n"a[2]"\n"a[3]};c=1;}}'
#take central part of the jump
%.jumps.2: %.jumps
	awk '$$1==2' $< | cut -f 1 --complement > $@
#for i in /scratch/devel/eraineri/gimli_test/*.jumps ; do echo ${i}; awk '$1==2' ${i} | cut -f 1 --complement > ${i}.2; done
#intersect jumps from different samples

monocytes.jumps.2.intersection.txt: $(GTEST)/C000S5A1bs_cpg.gimli.jumps.2 $(GTEST)/C0010KA2bs_cpg.gimli.jumps.2 $(GTEST)/C001UYA3bs_cpg.gimli.jumps.2 $(GTEST)/C004SQ51_cpg.gimli.jumps.2
	bedtools multiinter -i <(awk '$$NF==10' $(GTEST)/C000S5A1bs_cpg.gimli.jumps.2) <(awk '$$NF==10' $(GTEST)/C0010KA2bs_cpg.gimli.jumps.2) <(awk '$$NF==10' $(GTEST)/C001UYA3bs_cpg.gimli.jumps.2) <(awk '$$NF==10' $(GTEST)/C004SQ51_cpg.gimli.jumps.2)  | awk '$$4==4' > $@ 
#intersect with annotation

gencode.v15.tss : $(ANNOTATION)
	zcat $(ANNOTATION) | grep -ve "^#" |\
	awk '$$3=="gene"' |\
	awk 'BEGIN{OFS="\t"}{c=($$7=="+")?$$4:$$5;print $$1,c,c,$$7}' >$@

monocytes.jumps.2.intersection.tss.txt: monocytes.jumps.2.intersection.txt gencode.v15.tss
	bedtools intersect -a $< -b gencode.v15.tss  -wao |\
	awk '$10!="."' |\
	awk 'BEGIN{OFS="\t"}{print $0,$1"_"$2"_"$3}' |\
	sort -k15,15 |\
	uniq -f 14  > $@ 

monocytes.jumps.2.intersection.annotated.txt:  monocytes.jumps.2.intersection.txt $(ANNOTATION)
	bedtools intersect -a monocytes.jumps.2.intersection.txt -b <(zcat $(ANNOTATION) | awk 'BEGIN{OFS="\t"}{print $$1,$$4,$$5,$$3,$$14}') -wao > $@ 
#look for unique intersections
#cat $SCRATCH/gimli_test/monocytes.jumps.2.intersection.annotated.txt | awk '{print $0"\t"$2"_"$3"_"$11"_"$12}' | uniq -f 14
#cat $SCRATCH/gimli_test/monocytes.jumps.2.intersection.bkp.annotated.txt | awk 'BEGIN{OFS="\t"}{print $0,$2"_"$3}' | uniq -f 9 > $SCRATCH/gimli_test/monocytes.jumps.2.intersection.bkp.annotated.uniq.txt

######example 3 : dmrs 

MACROM0=$(SCRATCH)/henk/macrophages/M0
MONO=$(SCRATCH)/henk/monocytes


#####joining cmd

#monocytes

$(MONO)/mono.joined.txt:  $(MONO)/C000S5A1bs_cpg.txt.gz $(MONO)/C0010KA2bs_cpg.txt.gz $(MONO)/C004SQ51_cpg.txt.gz $(MONO)/C005PS51_cpg.txt.gz $(MONO)/S000RD54_cpg.txt.gz $(MONO)/S007G756_cpg.txt.gz
	gimli.join  <( zcat $(MONO)/C000S5A1bs_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}')  <( zcat $(MONO)/C0010KA2bs_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}')  <( zcat $(MONO)/C004SQ51_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}')  <( zcat $(MONO)/C005PS51_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}')  <( zcat $(MONO)/S000RD54_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}')  <( zcat $(MONO)/S007G756_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}') > $@


#macrophages M0

$(MACROM0)/M0.joined.txt: $(MACROM0)/C005VG51_cpg.txt.gz $(MACROM0)/S001S751_cpg.txt.gz $(MACROM0)/S0039051_cpg.txt.gz $(MACROM0)/S00BHQ51_cpg.txt.gz $(MACROM0)/S00DVR51_cpg.txt.gz
	
	gimli.join  <( zcat $(MACROM0)/C005VG51_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}')  <( zcat $(MACROM0)/S001S751_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}')  <( zcat $(MACROM0)/S0039051_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}')  <( zcat $(MACROM0)/S00BHQ51_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}')  <( zcat $(MACROM0)/S00DVR51_cpg.txt.gz | awk 'BEGIN{OFS="\t"}{print $$1,$$2,$$6,$$7}') > $@ 

##########

MONOGEC:=$(MONO)/mono.joined.gec.txt
M0GEC:=$(MACROM0)/M0.joined.gec.txt
GIMLIMONO:=$(MONO)/mono.joined.gec.gimli
GIMLIM0:=$(MACROM0)/M0.joined.gec.gimli
OLAP10:=$(SCRATCH)/henk/mono.M0.l10.olap
OLAP100:=$(SCRATCH)/henk/mono.M0.l100.olap
MONOWINDOWS10:=$(SCRATCH)/henk/mono.l10.olap.windows
M0WINDOWS10:=$(SCRATCH)/henk/M0.l10.olap.windows
DMR10:=$(SCRATCH)/henk/mono.M0.l10.dmr
MONOWINDOWS100:=$(SCRATCH)/henk/mono.l100.olap.windows
M0WINDOWS100:=$(SCRATCH)/henk/M0.l100.olap.windows
DMR100:=$(SCRATCH)/henk/mono.M0.l100.dmr

$(MONOGEC): /scratch/devel/eraineri/henk/monocytes/mono.joined.txt
	cat $^  | gimli.gec > $@ 2>$(MONOGEC).log

$(M0GEC): /scratch/devel/eraineri/henk/macrophages/M0/M0.joined.txt
	cat $^ | gimli.gec > $@ 2>$(M0GEC).log

$(GIMLIMONO): $(MONOGEC)
	gimli $(MONOGEC) 1:10:100:1000 > $@ 2>$(GIMLIMONO).log

$(GIMLIM0): $(M0GEC)
	gimli $(M0GEC) 1:10:100:1000 > $@ 2>$(GIMLIM0).log
	
$(OLAP10): $(GIMLIMONO) $(GIMLIM0)
	bedtools intersect -a <(awk '$$NF==10' $(GIMLIMONO) | awk '{print $$1"\t"$$2"\t"$$3}') -b <(awk '$$NF==10' $(GIMLIM0) | awk '{print $$1"\t"$$2"\t"$$3}' ) -wao | awk '{lb=($$2>$$5)?$$2:$$5;rb=($$3>$$6)?$$6:$$3;print $$0"\t"lb"\t"rb}'   > $@ 

$(OLAP100): $(GIMLIMONO) $(GIMLIM0)
	bedtools intersect -a <(awk '$$NF==100' $(GIMLIMONO) | awk '{print $$1"\t"$$2"\t"$$3}') -b <(awk '$$NF==100' $(GIMLIM0) | awk '{print $$1"\t"$$2"\t"$$3}' ) -wao | awk '{lb=($$2>$$5)?$$2:$$5;rb=($$3>$$6)?$$6:$$3;print $$0"\t"lb"\t"rb}'   > $@ 

$(MONOWINDOWS10): $(MONOGEC) $(OLAP10)
	cat $(MONOGEC) | gimli.windows <(awk '$$NF>0' $(OLAP10) | awk '{print $$1"\t"$$8"\t"$$9}'  ) 2>$@.log > $@ 

$(M0WINDOWS10): $(M0GEC) $(OLAP10)
	cat $(M0GEC) | gimli.windows <(awk '$$NF>0' $(OLAP10) | awk '{print $$1"\t"$$8"\t"$$9}' ) 2>$@.log > $@ 

$(DMR10): $(MONOWINDOWS10)  $(M0WINDOWS10)
	paste $^ | awk '$$4!="-" && $$13!="-" && $$4>4 && ($$8+$$17)>0' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$7"\t"$$8"\t"$$13"\t"$$16"\t"$$17"\t"($$7-$$16)/sqrt($$8/$$4+$$17/$$13)}' > $@

$(MONOWINDOWS100): $(MONOGEC) $(OLAP100)
	cat $(MONOGEC) | gimli.windows <(awk '$$NF>0' $(OLAP100) | awk '{print $$1"\t"$$8"\t"$$9}'  ) 2>$@.log > $@ 

$(M0WINDOWS100): $(M0GEC) $(OLAP100)
	cat $(M0GEC) | gimli.windows <(awk '$$NF>0' $(OLAP100) | awk '{print $$1"\t"$$8"\t"$$9}' ) 2>$@.log > $@ 

$(DMR100): $(MONOWINDOWS100)  $(M0WINDOWS100)
	paste $^ | awk '$$4!="-" && $$13!="-" && $$4>4 && ($$8+$$17)>0' | awk '{print $$1"\t"$$2"\t"$$3"\t"$$4"\t"$$7"\t"$$8"\t"$$13"\t"$$16"\t"$$17"\t"($$7-$$16)/sqrt($$8/$$4+$$17/$$13)}' > $@



~/Desktop/meth_data/C001UYA3bs.C000S5A1bs.gimli.1000.counts.1 :  /home/emanuele/Desktop/meth_data/C001UYA3bs_cpg.txt.gz ~/Desktop/meth_data/C001UYA3bs.C000S5A1bs.gimli.1000.olap
	zcat /home/emanuele/Desktop/meth_data/C001UYA3bs_cpg.txt.gz | awk '{print
	$$1"\t"$$2"\t"$$6"\t"$$7}' |  ./gimli.windows <(awk '{print
	$$1"\t"$$(NF-1)"\t"$$NF}'
	~/Desktop/meth_data/C001UYA3bs.C000S5A1bs.gimli.1000.olap) > $@ 


~/Desktop/meth_data/C001UYA3bs.C000S5A1bs.gimli.1000.counts: ~/Desktop/meth_data/C001UYA3bs.C000S5A1bs.gimli.1000.counts.1 ~/Desktop/meth_data/C001UYA3bs.C000S5A1bs.gimli.1000.counts.2
	paste $^ | grep -v "-" | awk '{print $$6"\t"$$7"\t"$$17"\t"$$18}' > $@ 


%.lrt: %.counts
	Rscript likratiotest.R $< > $@

.PHONY : clean

clean: 
	rm -f $(MONO)/mono.joined.txt $(MACROM0)/M0.joined.txt $(MONOGEC) $(M0GEC) $(GIMLIMONO) $(GIMLIM0) $(OLAP10) $(MONOWINDOWS10) $(M0WINDOWS10) $(DMR10)   
	rm -f $(OLAP100) $(MONOWINDOWS100) $(M0WINDOWS100) $(DMR100)
