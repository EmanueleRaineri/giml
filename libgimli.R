load.segments<-function(fname){
	#returns a data frame
	seg<-read.table(fname,stringsAsFactors=F)
	return(seg)
}

plot.segments<-function( seg ,  lb , ub, col ){
	for (i in 1:nrow(seg)){
		#vector of xs followed by vector of ys
		lines(c(seg[i,2],seg[i,3]),c(seg[i,"mletheta"],seg[i,"mletheta"]),col=colour,lwd=2)
	}
}

theta.segment<-function(nc,c){
	sum.nc<-sum(ns)
	sum.c<-sum(c)
	if (sum.nc + sum.c==0) {stop("undefined position")}
	theta<-sum.nc/(sum.nc+sum.c)
	return(theta)	
}

lik.segment<-function(nc,c){
	theta<-	theta.segment(nc,c);
	lik<-sum(binom(x+nc,size=nc+c,p=theta,log=T))
	return(lik)
}

merged.lik<-function(nc1,c1,nc2,c2){
	return(lik.segment(nc1+nc2,c1+c2))
}

delta.lik<-function(nc1,c1,nc2,c2){
	ml<-merged.lik(nc1,c1,nc2,c2)
	l1<-lik.segment(nc1,c1)
	l2<-lik.segment(nc2,c2)
	return (ml-l1-l2)
}

