import sys

f1=open(sys.argv[1],'w')
ldmr=int(sys.argv[2])



print >>f1,"chr1\t%d\t%d\t%f\t%d"%(1,int(1e6),0.9,30)
a=1
b=2
for i in range(a,b):
	print >>f1, "chr1\t%d\t%d\t%f\t%d"%\
		(int(i*1e6)+1,int(i*1e6)+ldmr,0.7,30)
	print >>f1, "chr1\t%d\t%d\t%f\t%d"%\
	(int(i*1e6)+ldmr+1,(i+1)*int(1e6),0.9,30)	
print >>f1, "chr1\t%d\t%d\t%f\t%d"%(b*int(1e6)+1,int(250*1e6),0.9,30)

f1.close()
