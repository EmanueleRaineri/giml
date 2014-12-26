g199<-read.table("G199.chr1.correla.txt",stringsAsFactors=F)
g200<-read.table("G200.chr1.correla.txt",stringsAsFactors=F)
g201<-read.table("G201.chr1.correla.txt",stringsAsFactors=F)
g202<-read.table("G202.chr1.correla.txt",stringsAsFactors=F)

cols=rainbow(4)
postscript("out.correla.eps")
plot((g199$V1+g199$V2)/2,g199$V8,type='b',col=cols[1],
	xlab="distance",ylab="correlation")
lines((g200$V1+g200$V2)/2,g200$V8,type='b',col=cols[2])
lines((g201$V1+g201$V2)/2,g201$V8,type='b',col=cols[3])
lines((g202$V1+g202$V2)/2,g202$V8,type='b',col=cols[4])
grid()
dev.off()
