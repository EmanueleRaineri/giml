library(IRanges)
libgimli = new.env()
#math
libgimli$loglik<-function(nonconv , conv){
	nc<-sum(nonconv)
	co<sum(conv)
	theta<-nc/(co+nc)
	dbinom(x=nc,size=nc+co,p=theta,log=T)+
	nc*theta+co*(1-theta)
}

libgimli$totloglik<-function(nonconv,conv,ls){
	#ls is list of segments
}
#I/O
libgimli$load.segments<-function(fname){
	#returns a data frame
	seg<-read.table(fname,stringsAsFactors=F)
	names(seg)<-c("chrom","start","end","ncpgs","min","mle","max","loglik",
		"delta","lambda")
	return(seg)
}
#plotting
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

libgimli$plot.dmr<-function( dmr ,  lb , ub, color1 , color2 ){	
	plot( 1 , type="n", xlab="", ylab="", 
	xlim=c(lb, ub), 
	ylim=c(-0.1, 1.1)
	)
	for (i in 1:nrow(dmr)){
		if ( dmr[i,"start1"] > ub || dmr[i,"end1"] < lb ) next
		lines(c(dmr[i,"start1"],dmr[i,"end1"]),
		c(dmr[i,"mle1"],dmr[i,"mle1"]),
		col=color1,lwd=2)
		lines(c(dmr[i,"start2"],dmr[i,"end2"]),
		c(dmr[i,"mle2"],dmr[i,"mle2"]),
		col=color2,lwd=2)
	}
}
#
libgimli$plot.intersected.dmr<-function( dmr ,  lb , ub, color1, color2 ){	
	plot( 1 , type="n", xlab="", ylab="", 
	xlim=c(lb, ub), 
	ylim=c(-0.1, 1.1),
	bty="n"
	)
	grid()
	for (i in 1:nrow(dmr)){
		if ( dmr[i,"start1"] > ub || dmr[i,"end1"] < lb ) next
		a<-	dmr[i,"start1"]
		b<- dmr[i,"end1"]
		c<- dmr[i,"start2"]
		d<- dmr[i,"end2"]
		xl<-max(a,c)
		xr<-min(b,d)
		y1<-dmr[i,"mle1"]
		y2<-dmr[i,"mle2"]
		lines(c(xl,xr), c(y1,y1), col=color1, lwd=3)
		lines(c(xl,xr), c(y2,y2), col=color2, lwd=3)
		lines(c(xl,xl) , c(0,max(y1,y2)), col="lightgreen", lty="dotted",lwd=1)
		lines(c(xr,xr) , c(0,max(y1,y2)), col="lightgreen", lty="dotted",lwd=1)
	}
}
#
libgimli$plot.multiple.segments<-function( list_of_seg ,  lb , ub, colors ){
	plot( 1 , type="n", xlab="", ylab="", 
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
	if ( sum.nc + sum.c == 0 ) {stop("undefined position")}
	theta<-sum.nc/(sum.nc+sum.c)
	return(theta)	
}
#
libgimli$load.meth<-function(fname){
	meth<-read.table(fname,stringsAsFactors=F)
	names(meth)<-c("chrom","pos",
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
libgimli$rd.of.bed<-function(df){
	pos<-IRanges(df[,2],df[,3])
	RangedData(pos)
}
#
libgimli$rd.of.meth<-function(df,col.pos,col.values){
	#second columns of df must be a genomic coordinate
	pos<-IRanges(df[,col.pos],df[,col.pos])
	values<-df[,col.values]
	RangedData(pos,values)
}
#
libgimli$intersect.with.bed<-function( bed.rd , meth.rd ){
	res<-matrix(1:(nrow(bed.rd)*6),nrow=nrow(bed.rd))
	ov.bed.meth     <- findOverlaps(bed.rd,meth.rd)
	ov.bed.meth.mat <- as.matrix(ov.bed.meth)
	lb  <- min(ov.bed.meth.mat[,1])	
	ub  <- max(ov.bed.meth.mat[,1])
	keys<-unique(ov.bed.meth.mat[,1])
	for ( i in 1:nrow(bed.rd)){
		if (i %in% keys){
			#print(i)
			slice<-meth.rd[ov.bed.meth.mat[ov.bed.meth.mat[,1]==i,2],]
			res[i,1] <- bed[i,2]
			res[i,2] <- bed[i,3]
			res[i,3] <- nrow(slice)
			mv <- slice$nc/(slice$nc+slice$c)
			res[i,4]<-min(mv)
			res[i,5]<-median(mv)
			res[i,6]<-max(mv)
		} else {
			res[i,1] <- bed[i,2]
			res[i,2] <- bed[i,3]
			res[i,3] <- 0 
			res[i,4]<- NA 
			res[i,5]<- NA
			res[i,6]<- NA
		}
	}
	names(res)<-c("start","end","ncpg","min","median","max")
	res
}
#
while("libgimli" %in% search())
  detach("libgimli")
attach(libgimli)
