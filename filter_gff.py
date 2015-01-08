import sys

genes1tss=open("gencode.chr1.genes.unique.tss").readlines()
genes1tss=[g.strip() for g in genes1tss]
gff = open("gencode.v19.TSS.notlow.chr1.gff")

for line in gff:
	line=line.strip()
	fields=line.split()
	if (fields[9] in genes1tss):
		print "%s\t%s\t%s\t%s"%(fields[0],fields[3],fields[4],fields[6]) 
