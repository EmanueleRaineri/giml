infile<-"G199_cpg.segmeth.chr10.in.pmd"
infile.segs<-"G199_cpg.segmeth.chr10.in.pmd.gimli.mle"
methyl<-read.table(infile,stringsAsFactors=F,sep=" ")
segments<-read.table(infile.segs,stringsAsFactors=F,sep="\t")
names(segments)<-c("from","to","ncpg","mintheta","mletheta","maxtheta","loglik","delta","lambda")

plot.segmentation<-function( methyl , segments , lambda , lb , ub, colour ){
	seg<-segments[segments$lambda==lambda,]
	#lb<-min(seg$from)
	#ub<-max(seg$to)
	plot(methyl$V2,methyl$V3/(methyl$V3+methyl$V4),xlim=c(lb,ub),pch=19,cex=0.2,main=paste("lambda=",lambda))
	abline(v=c(lb,ub),col="blue")
	for (i in 1:nrow(seg)){
		#vector of xs followed by vector of ys
		lines(c(seg[i,1],seg[i,2]),c(seg[i,"mletheta"],seg[i,"mletheta"]),col=colour,lwd=2)
	}
}
colors=rainbow(3)

par(mfrow=c(3,1))
plot.segmentation(methyl,segments,10,1.1885e8,1.1893e8,col=colors[1])
plot.segmentation(methyl,segments,100,1.1885e8,1.1893e8,col=colors[2])
plot.segmentation(methyl,segments,1000,1.1885e8,1.1893e8,col=colors[3])


