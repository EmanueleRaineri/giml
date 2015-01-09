import sys

"""chr1	10471	713512	296	0.4000	0.9122	1.0000	-433.2	-3559	100.0000"""

#all = open("C004GD51_cpg.chr1.gimli.100").readlines()
all = open(sys.argv[1]).readlines()

for i in range(len(all)):
	all[i]=all[i].strip()



for i in range(len(all)-2):
	min1 = all[i].split()[4]
	max2 = all[i+1].split()[6]
	min3 = all[i+2].split()[4]
	if (max2<min1 and max2<min3):
		print >>sys.stderr,"***\n", all[i],"\n",all[i+1],"\n",all[i+2]
		print >>sys.stdout,all[i+1]


