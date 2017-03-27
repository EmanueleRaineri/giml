#!/usr/bin/env python
import sys
# takes in input a  methylation counts file
# and a giml output file

try:
	methfn=sys.argv[1]
	gimlfn=sys.argv[2]
	lam=sys.argv[3]
except:
	sys.stderr.write("usage: <meth-file> <giml-file> <lambda>\n")
	sys.exit(1)

gimlf=open(gimlfn)
arrows=""
for line in gimlf:
	fields=line.split("\t")
	if (fields[3]!=lam):continue
	level=float(fields[6])
	col="blue"
	arrows=arrows+"set arrow from %s,%s to %s,%s nohead lw 1.5 lc rgb '%s'\n"%\
	(fields[1],level,fields[2],level,col)

sys.stdout.write(arrows)
sys.stdout.write("set key bmargin box\nset yra[0:1]\nset grid\n")
cmd="plot \"%s\" u 2:($3/($3+$4)) pt 7 ps 0.7  lc rgb \"%s\"\n"%(methfn,"red")
sys.stdout.write(cmd)
