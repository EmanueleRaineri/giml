filetss=open("C004GD51_cpg.chr1.gimli.tss.filtered.bed")
filerpkm=open("C004GD12.rpkm")

tablerpkm={}

"""chr1	41157320	41237275	+	ENSG00000066136.14	10.5251"""

for line in filerpkm:
	line=line.strip()
	fields=line.split()
	
	if (fields[3]=="+"):
		tss=int(fields[1])
	else:
		tss=int(fields[2])
	
	rpkm = float(fields[-1])

	tablerpkm[tss]=rpkm

filerpkm.close()

for line in filetss:
	
	line   =  line.strip()
	fields =  line.split()
	
	tss= int(fields[1])
	mle=float(fields[9])
	print tss,mle,tablerpkm[tss]

filetss.close()
