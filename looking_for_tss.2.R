source("libgimli.R")
if (!exists("meth")){
	meth<-libgimli$load.meth("~/Desktop/meth_data/C004GD51_cpg.chr1.txt.gz")
}

bed<-read.table("S000RD13.tss.bed",stringsAsFactors=F)
bed<-bed[1:10,]

meth.rd<-libgimli$rd.of.meth(meth,2,c(6,7))


colors=rainbow(10)
#par(mfcol=c(5,2))

for ( i in 1:2 ) {
	tss      <- bed[i,2]
	tss.meth <- libgimli$take.slice(meth.rd,tss-1000,tss+1000)
	if (nrow(tss.meth)==0) next
	xcoord <- start(tss.meth)
	plot(xcoord,
		tss.meth$nc/(tss.meth$nc+tss.meth$c),
		type='b',
		col="black",
		ylim=c(0,1),main=paste(i,tss),
	)
	abline(v=tss,col="red")
	readline()
}

tss<-226249552
tss.meth <- libgimli$take.slice(meth.rd,tss-1000,tss+1000)
	xcoord <- start(tss.meth)
	plot(xcoord,
		tss.meth$nc/(tss.meth$nc+tss.meth$c),
		type='b',
		col="black",
		ylim=c(0,1),main=paste(i,tss),
	)
	abline(v=tss,col="red")


