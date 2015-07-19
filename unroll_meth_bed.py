import sys

def bin_search(v,x):
	#returns the index of x, or an index close by
	#if x is not in v.
	uidx=len(v)-1
	lidx=0
	while True:
		midx=(uidx+lidx)/2
		if ( x > v[midx] ):
			lidx=midx
		elif (x < v[midx]):
			uidx=midx
		else:
			return midx
		if ((uidx-lidx)<=1):
			return midx

bedfile=open(sys.argv[1])

lefts  = []
rights = []
meth   = []
rd     = []

for line in bedfile:
	line      = line.strip()
	fields    = line.split()
	chrom     = fields[0]
	lefts.append(int( fields[1] )) 
	rights.append(int( fields[2] ))
	meth.append(float(fields[3]))
	rd.append( int(fields[4]))
bedfile.close()

for line in sys.stdin:
	line = line.strip()
	x    = int(line)
	idx  = bin_search(lefts,x)
	if (idx<(len(lefts)-1)):
		if (lefts[idx+1]<=x):
			idx=idx+1
		
	while (lefts[idx]>x):
		idx=idx-1
		if (idx==0):break
	if ( lefts[idx] <= x and rights[idx]>=x  ):
		print chrom,x,meth[idx],rd[idx]
	else:
		print "%d is not covered by the bed file"%x
		sys.exit(1)

