#library(IRanges)
#library(GenomicRanges)
#library(plyr)

libgimli = new.env()

#math
libgimli$loglik<-function(nonconv , conv){
	nc<-sum(nonconv)
	co<-sum(conv)
	theta<-nc/(co+nc)
	dbinom(x=nc,size=nc+co,p=theta,log=T)
}

libgimli$llik<-function(nonconv , conv, theta){
	dbinom(x=nonconv,size=nonconv+conv,p=theta,log=T)
}


libgimli$smooth.meth<-function( meth , bw=1000 , n.points=1000) {
	res<-ksmooth(meth$pos, 
		meth$meth, 
		kernel = "normal", 
		bandwidth = bw, n.points = n.points)
	return(res)
}

libgimli$cpg.density<-function(meth, lb=meth[1,2] , ub=meth[nrow(meth),2] ){
	tmp<-meth[meth$pos>=lb && meth$pos<=ub,]
	return(nrow(tmp)/(tmp[nrow(tmp),2]-tmp[1,2]))
}

libgimli$beta.diff<-function(a1,b1,a2,b2){
	integrate(function(x) dbeta(x,a1,b1)*pbeta(x,a2,b2),0,1)[[1]]
}

libgimli$methyl.diff<-function(nc1,c1,nc2,c2){
	libgimli$beta.diff (nc1+1 , c1+1 , nc2+1 , c2+1 )
}

libgimli$sim.counts<-function(d1,d2,size){
	theta <- runif( size )
	nc1   <- mapply( function( t ){ rbinom( 1, d1, t ) }, theta )
	c1    <- ( d1 - nc1 )
	nc2   <- mapply( function( t ){ rbinom( 1, d2, t ) }, theta )
	c2    <- ( d2 - nc2 )
	d<-data.frame(nc1=nc1,c1=c1,nc2=nc2,c2=c2)
}

libgimli$fisher.diff<-function(nc1,c1,nc2,c2){
		ft<-fisher.test(matrix(c(nc1,c1,nc2,c2),nrow=2),alternative="greater")
		ft$p.value
}

libgimli$theta.segment<-function( nc , c ){
	sum.nc<-sum(ns)
	sum.c<-sum(c)
	if ( sum.nc + sum.c == 0 ) {stop("undefined position")}
	theta<-sum.nc/(sum.nc+sum.c)
	return(theta)	
}

libgimli$lik.segment<-function(nc,c){
	theta<-	theta.segment(nc,c);
	lik<-sum(binom(x+nc,size=nc+c,p=theta,log=T))
	return(lik)
}

libgimli$merged.lik<-function(nc1,c1,nc2,c2){
	return(lik.segment(nc1+nc2,c1+c2))
}

libgimli$delta.lik<-function(nc1,c1,nc2,c2){
	ml<-merged.lik(nc1,c1,nc2,c2)
	l1<-lik.segment(nc1,c1)
	l2<-lik.segment(nc2,c2)
	return ( ml - l1 - l2 )
}

#I/O

libgimli$load.segments<-function(fname){
	#returns a data frame
	seg<-read.table( fname , stringsAsFactors=F )
	names(seg)<-c( "chrom","start","end","ncpgs",
		"min","mle","max","mean","var","loglik",
		"delta","lambda" )
	return(seg)
}

libgimli$load.meth<-function(fname){
	meth<-read.table(fname,stringsAsFactors=F)
	names(meth)<-c( "chrom" , "pos" ,
	"phred" , "meth" , "var" , "nc", "c" ,
	"cov1" , "cov2"
	)
	return( meth )
}

libgimli$load.stripped<-function(fname){
	meth.counts<-read.table(fname,stringsAsFactors=F)
	names(meth.counts)<-c("chrom" , "pos" ,"nc", "c")
	return ( meth.counts )
}


libgimli$load.dmr<-function( fname ){
	dmr<-read.table( fname , stringsAsFactors=F )
	names( dmr ) <- c( "chrom",
	"start", "end", "ncpgs1","ncpgs2",
	"mle1", "mle2", "mle12",
	"lik1", "lik2", "lik12",
	"deltalik", "pval" )
	return( dmr )
}

libgimli$read.binary.table <- function( fname ){
	bin.table <- read.table( fname, skip=1, header=T, stringsAsFactors=F )
	return( bin.table )
}

libgimli$count.conf <- function(df){
	counts<-plyr::count(df)
	return(counts)
}

libgimli$dec.of.bin <- function(df){
	decoded.bin <- apply(bin,MARGIN=1,
		FUN=function(x) {strtoi(paste(x,collapse=""),base=2)})
	return (decoded.bin)
}

libgimli$count.trans <- function(numvect){
	tr2<-table(paste(head(numvect,-1),tail(numvect,-1)))
	return (tr2)
}
#plotting
libgimli$plot.segments<-function( seg ,  lb , ub, col ){
	plot(1, type="n", xlab="", ylab=expression(theta), xlim=c(lb, ub), ylim=c(-0.1, 1.1))
	for (i in 1:nrow(seg)){
		#vector of xs followed by vector of ys
		if (seg[i,2]>ub) break
		lines(c(seg[i,2],seg[i,3]),c(seg[i,"mle"],seg[i,"mle"]),xlim=c(lb,ub),
			col=col,lwd=2)
	}
}

libgimli$plot.segments.box<-function( seg ,  lb , ub, col ){
	#note : postscript does not support partial transparency
	plot(1, type="n", xlab="", ylab="", 
	xlim=c(lb, ub), 
	ylim=c(-0.1, 1.1)
	)
	for (i in 1:nrow(seg)){
		#vector of xs followed by vector of ys
		#lines(c(seg[i,2],seg[i,3]),c(seg[i,"mle"],seg[i,"mle"]),col=col,lwd=2)
		if (seg[i,3]<lb) next	
		if (seg[i,2]>ub) break
		rect(seg[i,2],seg[i,"min"],seg[i,3],seg[i,"max"],
		density=30,col=col)
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
		lines(c(xl,xl) , c(0,max(y1,y2)), col="green", lty="solid",lwd=2)
		lines(c(xr,xr) , c(0,max(y1,y2)), col="green", lty="solid",lwd=2)
	}
}

libgimli$plot.intersected.dmr.box<-function( dmr ,  lb , ub, color1, color2 ){	
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
		ybottom1<-dmr[i,"min1"]
		ytop1<-dmr[i,"max1"]
		ybottom2<-dmr[i,"min2"]
		ytop2<-dmr[i,"max2"]
		#y1<-dmr[i,"mle1"]
		#y2<-dmr[i,"mle2"]
		#lines(c(xl,xr), c(y1,y1), col=color1, lwd=3)
		#lines(c(xl,xr), c(y2,y2), col=color2, lwd=3)
		#lines(c(xl,xl) , c(0,max(y1,y2)), col="green", lty="solid",lwd=2)
		#lines(c(xr,xr) , c(0,max(y1,y2)), col="green", lty="solid",lwd=2)
		rect(xl,ybottom1,xr,ytop1,density=10,col=color1)
		rect(xl,ybottom2,xr,ytop2,density=10,col=color2)
	}
}



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


libgimli$take.slice<-function(rd,lb,ub){
	q<-IRanges(lb,ub)
	ov<-findOverlaps(RangedData(q),rd)
	rd[as.matrix(ov)[,2],]
}

libgimli$rd.of.bed<-function(df){
	pos<-IRanges(df[,2],df[,3])
	RangedData(pos)
}

libgimli$gr.of.bed<-function(df){
	cat("GRanges of bed\n")
	GRanges(df[[1]],IRanges(df[,2],df[,3]))		
}

libgimli$rd.of.meth<-function(df,col.pos,col.values){
	#second columns of df must be a genomic coordinate
	pos<-IRanges(df[,col.pos],df[,col.pos])
	values<-df[,col.values]
	RangedData(pos,values)
}

libgimli$gr.of.meth<-function(df,col.pos,col.values){
	cat("GRanges of methylation values\n")
	pos<-IRanges(df[,col.pos],df[,col.pos])
	values<-df[[col.values]]
	GRanges(df[[1]],pos,score=values)
}

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

libgimli$genomic.intersect.with.bed<-function( bed, bed.gr , meth.gr ){
	res<-matrix(0,nrow=length(bed.gr),ncol=7)
	ov.bed.meth     <- findOverlaps(bed.gr,meth.gr)
	if (length(ov.bed.meth)==0){
		stop("no overlaps\n")
	}
	ov.bed.meth.mat <- as.matrix(ov.bed.meth)
	lb  <- min(ov.bed.meth.mat[,1])	
	ub  <- max(ov.bed.meth.mat[,1])
	keys<-unique(ov.bed.meth.mat[,1])
	for ( i in 1:length(bed.gr)){
		if ( i %% 100 == 0) {
			cat("processing interval:",i,"\n")
		}
		res[i,1]<-bed[i,1]	
		res[i,2]<-bed[i,2]	
		res[i,3]<-bed[i,3]	
		if (i %in% keys){
			sl<-meth.gr[ov.bed.meth.mat[ov.bed.meth.mat[,1]==i,2],]
			res[i,4] <- length(sl$score)
			mv <- sl$score
			res[i,5]<-min(mv)
			res[i,6]<-median(mv)
			res[i,7]<-max(mv)
		} else {
			res[i,4] <- 0 
			res[i,5]<- NA 
			res[i,6]<- NA
			res[i,7]<- NA
		}
	}
	res<-data.frame(res)
	names(res)<-c("chr","start","end","ncpg","min","median","max")
	res
}

while("libgimli" %in% search())
  detach("libgimli")
attach(libgimli)
