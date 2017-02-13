#!/usr/bin/env python
import numpy as np
import sys
# takes in input a  methylation counts file
# and a giml output file and a chrom and lambda

methfn=sys.argv[1]
gimlfn=sys.argv[2]
chrom=sys.argv[3]


meth=np.loadtxt(fname=methfn,dtype="S10,int,int,int")

meth=filter(lambda e : (e[2]+e[3] >0 and e[0]==chrom),meth)


pos=np.array([e[1] for e in meth])
theta=np.array([float(e[2])/(e[2]+e[3]) for e in meth])

gimlf = open(gimlfn)

for line in gimlf:
	line   = line.strip()
	fields = line.split("\t")
	start=int(fields[4])
	stop=int(fields[5])+1
	avg=np.mean(theta[start:stop])
	cnt=np.size(theta[start:stop])
	print line,"%.4f"%(avg)

