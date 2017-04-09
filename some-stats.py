#!/usr/bin/env python
import numpy as np
import sys
# takes in input a  methylation counts file
# and a giml output file and a chrom 

try:
	methfn=sys.argv[1]
	gimlfn=sys.argv[2]
	chrom=sys.argv[3]
	lam=float(sys.argv[4])
except:
	sys.stderr.write("Usage: <methfn> <gimlfn> <chrom> <lambda>\n")
	sys.exit(1)

sys.stderr.write("reading up meth file %s\n"%(methfn))
meth=np.loadtxt(fname=methfn,dtype="S10,int,int,int")
sys.stderr.write("filtering...\n")

meth=filter(lambda e : (e[2]+e[3] >0 and e[0]==chrom),meth)

pos=np.array([e[1] for e in meth])
theta=np.array([float(e[2])/(e[2]+e[3]) for e in meth])

gimlf = open(gimlfn)
sys.stderr.write("processing giml file %s...\n"%(gimlfn))

for line in gimlf:
	line   = line.strip()
	fields = line.split("\t")
	l=float(fields[3])
	if (l!=lam):
		continue
	start=int(fields[4])
	stop=int(fields[5])+1
	sl=theta[start:stop]
	avg=np.mean(sl)
	var=np.var(sl)
	cnt=np.size(sl)
	sys.stdout.write("%s\t%d\t%.4f\t%.4f\n"%(line,cnt,avg,var))

