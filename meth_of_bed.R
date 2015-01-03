source("libgimli.R")

bed.fname  <- "gencode.v19.TSS.notlow.chr1.promoters"
meth.fname <- "~/Desktop/meth_data/C004GD51_cpg.chr1.txt.gz"

if (!exists("meth")){
	meth<-libgimli$load.meth(meth.fname)
	bed<-read.table( bed.fname , stringsAsFactors=F )
}

meth.rd    <- libgimli$rd.of.meth(meth,2,c(6,7))
#lb<-bed[2,2]
#ub<-bed[2,3]
#print(lb)
#print(ub)
#lb<-lb-1000
#ub<-ub+1000
lb<-64431711-500
ub<-64432082+500
meth.slice <- libgimli$take.slice(meth.rd,lb,ub)
plot(start(meth.slice$ranges),meth.slice$nc/(meth.slice$nc+meth.slice$c),
	col="red",type='b',ylim=c(0,1))
abline(v=c(64431711,64432082),col="blue")	
