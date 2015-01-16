import sys
datadir=sys.argv[1]
ifile = open(datadir+"/C004GD51_cpg.chr1.gimli.tss.bed")

table={}
# values = [ width, ncpg, dm ]

for line in ifile:
	line=line.strip()
	fields=line.split()
	if (fields[4]=="."):continue
	tss=int(fields[1])
	#print tss
	width = int(fields[6])-int(fields[5])
	dm = float(fields[10])-float(fields[8])
	ncpg = int (fields[7])
	if ( tss in table ):
		if (dm<table[tss][3]
			and width>table[tss][1]):
			table[tss]=[line,width,ncpg,dm]
	else:
		table[tss]=[line,width,ncpg,dm]
	
ifile.close()

for tss in table.keys():
	print table[tss][0]	
