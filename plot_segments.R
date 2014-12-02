infile<-"G199_cpg.segmeth.chr10.in.pmd"
infile.segs<-"G199_cpg.segmeth.chr10.in.pmd.out"
methyl<-read.table(infile,stringsAsFactors=F,sep=" ")
segments<-read.table(infile.segs,stringsAsFactors=F,sep="\t")
names(segments)<-c("chrom","from","to","ncpg","mintheta","mletheta","maxtheta","loglik","delta","lambda")

plot.segmentation<-function( methyl , segments , lambda , lb , ub, colour ){
	seg<-segments[segments$lambda==lambda,]
	#lb<-min(seg$from)
	#ub<-max(seg$to)
	la<-lambda
	plot(methyl$V2,methyl$V3/(methyl$V3+methyl$V4),
	xlim=c(lb,ub),pch=19,cex=0.2,xlab=NULL,
	ylab=expression(theta),xaxt="n",main=paste("lambda=",la,sep=""))
	#abline(v=c(lb,ub),col="blue")
	for (i in 1:nrow(seg)){
		#vector of xs followed by vector of ys
		lines(c(seg[i,2],seg[i,3]),c(seg[i,"mletheta"],seg[i,"mletheta"]),col=colour,lwd=2)
	}
}
colors=rainbow(3)

par(mfrow=c(3,1),mai=c(0.3,0.8,0.3,0.5))
plot.segmentation(methyl,segments,10,1.1885e8,1.1893e8,col=colors[1])

plot.segmentation(methyl,segments,100,1.1885e8,1.1893e8,col=colors[2])
plot.segmentation(methyl,segments,1000,1.1885e8,1.1893e8,col=colors[3])
axis(1,at=seq(118850000,118930000,by=10000),labels=T)
