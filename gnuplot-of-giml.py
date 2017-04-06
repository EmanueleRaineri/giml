#!/usr/bin/env python
import sys
# takes in input a  methylation counts file
# and a giml output file

try:
	methfn=sys.argv[1]
	gimlfn=sys.argv[2]
	lam=sys.argv[3]
	mode=sys.argv[4]
except:
	sys.stderr.write(
		"usage: <meth-file> <giml-file> <lambda> <mode(all/segs)>\n")
	sys.exit(1)

gimlf=open(gimlfn)
arrows=""
xmin=float('+inf')
xmax=float("-inf")
for line in gimlf:
	fields=line.split("\t")
	if (fields[3]!=lam):continue
	level=float(fields[6])
	col="blue"
	if (float(fields[1])<xmin):
		xmin=float(fields[1])
	if (float(fields[2])>xmax):
		xmax=float(fields[2])
	arrows=arrows+"set arrow from %s,%s to %s,%s nohead lw 1.5 lc rgb '%s'\n"%\
	(fields[1],level,fields[2],level,col)

sys.stdout.write(arrows)
sys.stdout.write("set key bmargin box\nset yra[0:1]\nset grid\n")
sys.stdout.write("set xra[%f:%f]\n"%(xmin,xmax))
if ("all"==mode):
	cmd="plot \"%s\" u 2:($3/($3+$4)) pt 7 ps 0.7  lc rgb \"%s\" w dots\n"%\
		(methfn,"red")
elif ("segs"==mode):
	cmd="plot 1/0\n"
else:
	sys.stderr.write("bon't know about mode %s\n"%(mode))
	sys.exit(1)
sys.stdout.write(cmd)
