library(IRanges)
libgimli = new.env()
#
libgimli$load.segments<-function(fname){
	#returns a data frame
	seg<-read.table(fname,stringsAsFactors=F)
	names(seg)<-c("chrom","start","end","ncpgs","min","mle","max","loglik",
		"delta","lambda")
	return(seg)
}
#
libgimli$plot.segments<-function( seg ,  lb , ub, col ){
	plot(1, type="n", xlab="", ylab="", 
	xlim=c(lb, ub), 
	ylim=c(-0.1, 1.1)
	)
	for (i in 1:nrow(seg)){
		#vector of xs followed by vector of ys
		lines(c(seg[i,2],seg[i,3]),c(seg[i,"mle"],seg[i,"mle"]),col=col,lwd=2)
	}
}
#
libgimli$plot.dmr<-function( dmr ,  lb , ub, col1 , col2 ){
	
}
#
libgimli$plot.multiple.segments<-function( list_of_seg ,  lb , ub, colors ){
	plot(1, type="n", xlab="", ylab="", 
	xlim=c(lb, ub), 
	ylim=c(-0.1, 1.1)
	)
	for(j in 1:length(list_of_seg)){
		seg=list_of_seg[[j]]
		for (i in 1:nrow(seg)){
			#vector of xs followed by vector of ys
			lines(c(seg[i,2],seg[i,3]),c(seg[i,"mle"],seg[i,"mle"]),
			col=colors[j],lwd=2)
		}
	}
}
#
libgimli$theta.segment<-function(nc,c){
	sum.nc<-sum(ns)
	sum.c<-sum(c)
	if (sum.nc + sum.c==0) {stop("undefined position")}
	theta<-sum.nc/(sum.nc+sum.c)
	return(theta)	
}
#
libgimli$load.meth<-function(fname){
	meth<-read.table(fname,stringsAsFactors=F)
	names(meth)<-c("chrom","pos"
	"phred","meth","var","nc","c",
	"cov1","cov2"
	)
	return(meth)
}
#
libgimli$load.dmr<-function(fname){
	dmr<-read.table(fname,stringsAsFactors=F)
	names(dmr)<-c("chrom",
	"start1","end1","ncpgs1",
	"min1","mle1","max1",
	"start2","end2","ncpgs2",
	"min2","mle2","max2"
	)
	return(dmr)
}
#
libgimli$lik.segment<-function(nc,c){
	theta<-	theta.segment(nc,c);
	lik<-sum(binom(x+nc,size=nc+c,p=theta,log=T))
	return(lik)
}
#
libgimli$merged.lik<-function(nc1,c1,nc2,c2){
	return(lik.segment(nc1+nc2,c1+c2))
}
#
libgimli$delta.lik<-function(nc1,c1,nc2,c2){
	ml<-merged.lik(nc1,c1,nc2,c2)
	l1<-lik.segment(nc1,c1)
	l2<-lik.segment(nc2,c2)
	return (ml-l1-l2)
}
#
libgimli$take.slice<-function(rd,lb,ub){
	q<-IRanges(lb,ub)
	ov<-findOverlaps(RangedData(q),rd)
	rd[as.matrix(ov)[,2],]
}
#
libgimli$rd.of.meth<-function(df,col.pos,col.values){
	#second columns of df must be a genomic coordinate
	pos<-IRanges(df[,col.pos],df[,col.pos])
	values<-df[,col.values]
	RangedData(pos,values)
}
#
while("libgimli" %in% search())
  detach("libgimli")
attach(libgimli)
