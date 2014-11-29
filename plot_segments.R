source("segmentlib.R")
infile<-"G199_cpg.segmeth.chr10.in.pmd"
infile.segs<-"G199_cpg.segmeth.chr10.in.pmd.gimli"
methyl<-read.table(infile,stringsAsFactors=F,sep=" ")
segments

slib$plot.segmentation<-function(methyl,segments,lambda,lb,ub){
	seg<-segments[segments$lambda==lambda,]
	#lb<-min(seg$from)
	#ub<-max(seg$to)
	plot(methyl$V2,methyl$V3/(methyl$V3+methyl$V4),xlim=c(lb,ub))
	abline(v=c(lb,ub),col="blue")
	for (i in 1:nrow(seg)){
		lines(c(seg[i,1],seg[i,2]),c(seg[i,"theta"],seg[i,"theta"]),col="red",lwd=2)
	}
}
