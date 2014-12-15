import sys

metfile=sys.argv[1]
segfile=sys.argv[2]

if (segfile=="-"):
	segfile=sys.stdin
else:
	segfile=open(segfile)

print "set yra[-0.1:1.1]"

c=0
for line in segfile:
	line=line.strip()
	fields=line.split()
	pos1=int(fields[1])
	pos2=int(fields[2])
	mle=float(fields[5])
	print "set arrow %d from %d,%.4f to %d,%.4f nohead lc 2 lw 5"%(c+1,pos1,mle,pos2,mle)
	c+=1

print "plot \"%s\" u 2:($3/($3+$4)) w d lc rgb\"#ffffff\""%(metfile)

