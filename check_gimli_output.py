import sys

met_file=open(sys.argv[1])
gimli_file=open(sys.argv[2])

linec=0

sys.stderr.write("uploading methylation data...\n")
for line in met_file:
	linec=+1

met_file.close()
rsum=0

for line in gimli_file:
	fields = line.split()
	lam = float(fields[-1])
	if ( rsum == 0 ): #this is true for first line only
		oldlam = lam
	nocpg = int(fields[3])
	if ( nocpg <= 0 ):
		print "nocpg<=0"
		sys.exit(1)
	if (lam == oldlam):
		rsum += nocpg
	else:
		print "%.4f\t%d\t%d"%(oldlam,rsum,linec)
		rsum = nocpg
	oldlam = lam

print "%.4f\t%d\t%d"%(oldlam,rsum,linec)

