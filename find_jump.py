import sys

all = open("C004GD51_cpg.chr1.gimli.100").readlines()
print all

for i in range(len(all)-2):
	print all[i],all[i+1],all[i+2]
