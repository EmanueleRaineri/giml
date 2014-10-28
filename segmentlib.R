library(stats4)

slib = new.env()

#### beta ####

ll<-function(theta){
	r<-dbinom(x=nc,size=d,p=theta)
	-sum(log(r))
}

slib$moment.match <- function(mu,var){
	alpha <- (mu^2-mu^3-mu*var)/var
	beta  <- alpha*(1-mu)/mu
	c(alpha,beta)
}

#### L2 ####

slib$sdist<-function(score){
	max(score)-min(score)
}

slib$l2.dist<-function(score){
	mean.score<-mean(score)
	return(sum((score-mean.score)^2))
}

slib$make.segment.l2<-function( pos , score ){
	# pos, score are vectors
	mean.score<-mean(score)
	l2dist<-sum((score-mean.score)^2)
	list(pos = pos, 
	score = score, 
	mean.score=mean.score,
	l2dist = l2dist
	)
}

slib$join.segments.l2<- function( seg1 , seg2 ){
	pos<-c( seg1$pos , seg2$pos )
	score<-c( seg1$score , seg2$score )
	mean.score<-mean(score)
	l2dist<-sum((score-mean.score)^2)
	list(pos = pos, 
	score = score, 
	mean.score=mean.score,
	l2dist = l2dist
	)
}

slib$delta.segments.l2 <- function ( seg1, seg2, lambda ){
	#want delta to be negative
	score<-c( seg1$score , seg2$score )
	mean.score<-mean(score)
	joinedl2dist<-sum((score-mean.score)^2)
	#tmp.join<-join.segments.l2( seg1 , seg2 )
	#tmp.join$l2dist-seg1$l2dist-seg2$l2dist-lambda
	joinedl2dist-seg1$l2dist-seg2$l2dist-lambda
}

slib$all.pairs.l2 <- function(seg.list,lambda){
	les<-length(seg.list)
	seg.pairs<-rep(0,les-1)
	seg.pairs<-
		vapply( seq( 1 , les-1 ),
		FUN=function(i){slib$delta.segments.l2(seg.list[[i]],seg.list[[i+1]],lambda)},
		FUN.VALUE=1)
	return(seg.pairs)	
}

slib$update.segmentation.l2<-function(seg.list,all.pairs){
	idx <- which.min(all.pairs)
	if (all.pairs[[idx]]<0){
		seg.list[[idx]]<-slib$join.segments.l2(seg.list[[idx]],seg.list[[idx+1]])
		seg.list<-seg.list[-(idx+1)]
	}
	return(seg.list)
}

slib$seg.list.of.data.frame<-function(d){
	l<-list()
	for(i in 1:nrow(d)){
		pos <-d[i,2]
		score  <- d[i,3]
		l[[i]] <- slib$make.segment.l2( pos, score)   	
	}
	return(l)
}

slib$print.seg.list<-function(seg.list,lambda){
	les<-length(seg.list)
	#cat(les," segments at lambda:",lambda,"\n")
	for ( i in 1:les ){
		l<-seg.list[[i]]
		cat( min(l$pos), max(l$pos), l$mean.score, l$l2dist, lambda , "\n" )	
	}
}

slib$loop.over.lambda<-function(seg.list,all.lambda){
	segmentation<-data.frame()
	for (lambda in all.lambda){
		cat("lambda:",lambda,"\n")
		while(TRUE){
			if (length(seg.list)==1) break
			all.pairs<-slib$all.pairs.l2(seg.list,lambda)
			new.seg.list<-slib$update.segmentation.l2( seg.list , all.pairs )
			if (length( new.seg.list ) == length( seg.list ) ) break
			seg.list<-new.seg.list
		}
		idx<-nrow(segmentation)
		for (i in 1:length(seg.list)){
			li<-seg.list[[i]]
			idx<- idx+1
			segmentation[idx,1]<-min(li$pos)
			segmentation[idx,2]<-max(li$pos)
			segmentation[idx,3]<-length(li$pos)
			segmentation[idx,4]<-li$mean.score
			segmentation[idx,5]<-li$l2dist
			segmentation[idx,6]<-lambda
		}
	}
	names(segmentation)<-c("from","to","npoints","mean","l2dist","lambda")
	return(segmentation)
}

###############likelyhood############################

slib$make.segment<-function( pos , nconv, conv ){
	# pos, nonconv, conv are vectors
	d<- nconv + conv
	theta<-nconv/d
	p <- dbinom( x=nconv, size=d, p=median(theta))
	loglik <- sum(log(p/(1-p)))
	if (loglik==Inf) loglik<-100
	list(
		pos = pos, 
		nconv = nconv,
		conv = conv,
		depth = d,
		theta = theta,
		loglik = loglik)
}

slib$join.segments <- function( seg1 , seg2 ){
	pos<-c( seg1$pos , seg2$pos )
	nconv<-c( seg1$nconv , seg2$nconv )
	conv<-c( seg1$conv , seg2$conv )
	d<-c( seg1$d , seg2$d )
	theta<-c( seg1$theta , seg2$theta )
	p<-dbinom( x=nconv, size=d, p=median(theta))
	loglik <- log(p/(1-p)) 
	loglik[loglik==Inf]<-100
	loglik[loglik==-Inf]<- -100
	loglik<-sum(loglik)
	if (loglik==Inf) loglik<-100
	list( pos=pos , nconv=nconv , conv=conv , 
	depth=d , theta=theta , loglik=loglik)
}

slib$delta.segments <- function ( seg1, seg2, lambda ){
	#want delta to be positive
	tmp.join<-join.segments( seg1 , seg2 )
	tmp.join$loglik-seg1$loglik-seg2$loglik+lambda
}

slib$all.pairs <- function(seg.list,lambda){
	les<-length(seg.list)
	seg.pairs<-rep(0,les-1)
	tmp.max<--Inf
	idx.max<--1
	for ( i in 1:(les-1)){
		seg.pairs[[i]]<- slib$delta.segments(seg.list[[i]],seg.list[[i+1]],lambda)
		#cat(i," ",seg.pairs[[i]],"\n")
		if (seg.pairs[[i]]>=tmp.max) {
			tmp.max<-seg.pairs[[i]]
			idx.max<-i
		} 	
	}
	list( seg.pairs , idx.max )	
}

slib$update.segmentation<-function(seg.list,all.pairs, idx.max){
	if (all.pairs[[idx.max]]>0 && idx.max<length(seg.list)){
		seg.list[[idx.max]]<-slib$join.segments(seg.list[[idx.max]],seg.list[[idx.max+1]])
		seg.list<-seg.list[-(idx.max+1)]
	}
	return(seg.list)
}

slib$print.seg.list.lik<-function(seg.list,lambda){
	les<-length(seg.list)
	#cat(les," segments at lambda:",lambda,"\n")
	for ( i in 1:les ){
		l<-seg.list[[i]]
		cat( min(l$pos), 
		max(l$pos), 
		l$nconv, 
		l$conv, 
		l$d,
		min(l$theta),
		l$theta, 
		max(l$theta),
		l$loglik,
		lambda, "\n" )	
	}
	cat("*******\n")
}

slib$loop.over.lambda.lik<-function(seg.list,all.lambda){
	segmentation<-data.frame()
	for (lambda in all.lambda){
		cat("lambda:",lambda,"\n")
		all.pairs.max<-slib$all.pairs(seg.list,lambda)
		all.pairs<-all.pairs.max[[1]]
		idx.max<-all.pairs.max[[2]]
		tmp.max<-all.pairs[[idx.max]]
		#cat("all pairs initial:",all.pairs,"\n")
		while(TRUE){
			#print(all.pairs)
			if (length(seg.list)==1) break
			#cat("idx.max:",idx.max,"\n")
			#slib$print.seg.list.lik( seg.list , lambda )
			new.seg.list<-slib$update.segmentation( seg.list, all.pairs, idx.max )
			#slib$print.seg.list.lik( new.seg.list , lambda )
			if (length(new.seg.list)==length(seg.list) ) break
			seg.list<-new.seg.list
			all.pairs<-all.pairs[(-idx.max)]
			#print(all.pairs)
			#cat("idx.max:",idx.max,"\n")
			if (idx.max>1){
				all.pairs[[idx.max-1]]<-slib$delta.segments(seg.list[[idx.max-1]],seg.list[[idx.max]],lambda)
			}
			#print(all.pairs)
			if (idx.max<length(seg.list)){
				all.pairs[[idx.max]]<-slib$delta.segments(seg.list[[idx.max]],seg.list[[idx.max+1]],lambda)
			}
			#print(all.pairs)
			#cat("all pairs:",all.pairs,"\n")
			idx.max <- which.max(all.pairs)
		}
		idx<-nrow(segmentation)
		for (i in 1:length(seg.list)){
			li<-seg.list[[i]]
			idx<- idx+1
			segmentation[idx,1]<- min(li$pos)
			segmentation[idx,2]<- max(li$pos)
			segmentation[idx,3]<- length(li$pos)
			segmentation[idx,4]<- min(li$theta)
			segmentation[idx,5]<- median(li$theta)
			segmentation[idx,6]<- max(li$theta)
			segmentation[idx,7]<- li$loglik
			segmentation[idx,8]<- lambda
		}
	}
	names(segmentation)<-c("from","to","npoints","mintheta","theta","maxtheta", "loglik","lambda")
	return(segmentation)
}

slib$seg.list.of.data.frame.lik<-function(d){
	l<-list()
	for(i in 1:nrow(d)){
		pos <-d[i,2]
		nonconv  <- d[i,3]
		conv <- d[i,4]
		l[[i]] <- slib$make.segment( pos, nonconv, conv )   	
	}
	return(l)
}

slib$plot.segmentation<-function(methyl,segments,lambda,lb,ub){
	seg<-segments[segments$lambda==lambda,]
	#lb<-min(seg$from)
	#ub<-max(seg$to)
	plot(methyl$V2,methyl$V3/(methyl$V3+methyl$V4),xlim=c(lb,ub))
	for (i in 1:nrow(seg)){
		lines(c(seg[i,1],seg[i,2]),c(seg[i,4],seg[i,4]),col="red",lwd=2)
	}
}

###########################################

while("slib" %in% search())
  detach("slib")
attach(slib)
