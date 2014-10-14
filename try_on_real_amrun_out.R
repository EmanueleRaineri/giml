options(echo=TRUE)
source("segmentlib.R")


amrun<-read.table("G199.G200.cpg.chr1.2000.amrun")
amrun<-amrun[1:1000,2:3]

d<-data.frame(chrom=rep("chr1",1000),pos=amrun[,1],amrun=amrun[,2])
seg.list<-slib$seg.list.of.data.frame(d)


all.lambda<-c(0.01,0.02,0.05,0.1,0.2,0.5)
segmentation<-slib$loop.over.lambda(seg.list,all.lambda)

plot.segments<-function(amrun , lin.coords, xlim){
	plot(amrun[,1],amrun[,2],xlim=xlim,type='b')
	for (i in 1:nrow(lin.coords)){
		cat("line ",lin.coords[i,1],lin.coords[i,2],"\n")
		#lines(c(lin.coords[i,1],lin.coords[i,2]),
		#c(lin.coords[i,4],lin.coords[i,4]),col="red")
		#
		lines(c(lin.coords[i,1],lin.coords[i,2]),
		c(lin.coords[i,4]+lin.coords[i,5],lin.coords[i,4]+lin.coords[i,5]),col="red")
		#
		lines(c(lin.coords[i,1],lin.coords[i,2]),
		c(lin.coords[i,4]-lin.coords[i,5],lin.coords[i,4]-lin.coords[i,5]),col="red")
		#
		lines(c(lin.coords[i,1],lin.coords[i,1]),
		c(lin.coords[i,4]-lin.coords[i,5],lin.coords[i,4]+lin.coords[i,5]),col="red")
		#
		lines(c(lin.coords[i,2],lin.coords[i,2]),
		c(lin.coords[i,4]-lin.coords[i,5],lin.coords[i,4]+lin.coords[i,5]),col="red")
		#
		points(amrun[,1],amrun[,2],xlim=c(lin.coords[i,1],lin.coords[i,2]),col="green")
	}
}

xlim<-c(760000,770000)

seg01<-segmentation[segmentation$lambda==0.1 & segmentation$from>=xlim[1] & segmentation$to<=xlim[2],]
seg02<-segmentation[segmentation$lambda==0.2 & segmentation$from>=xlim[1] & segmentation$to<=xlim[2],]
seg05<-segmentation[segmentation$lambda==0.5 & segmentation$from>=xlim[1] & segmentation$to<=xlim[2],]

png("seg01.png",height=1000,width=1000)
plot.segments(amrun , seg01, xlim)
dev.off()
png("seg02.png",height=1000,width=1000)
plot.segments(amrun , seg02, xlim)
dev.off()
png("seg05.png",height=1000,width=1000)
plot.segments(amrun , seg05, xlim)
dev.off()


