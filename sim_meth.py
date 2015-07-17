import random
import sys
#chrom,pos,meth,rd
for line in sys.stdin:
	line=line.strip()
	fields=line.split()
	meth=float(fields[2])
	rd=int(fields[3])
	if (meth<0 or meth >1):
		sys.stderr.write("invalide methylation %f"%meth)
		continue
	nc = 0
	for i in range(rd):
		t = random.random()
		if (t<=meth):
			nc=nc+1
	print "%s\t%s\t%d\t%d"%(fields[0],fields[1],nc,rd-nc)
