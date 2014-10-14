slib = new.env()
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
		l[[i]] <- make.segment.l2( pos, score)   	
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
			new.seg.list<-slib$update.segmentation.l2(seg.list,all.pairs)
			if (length(new.seg.list)==length(seg.list) ) break
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
slib$make.segment <- function(pos,nconv,conv){
	# pos,nconv,conv are all vectors
	d      <- (nconv+conv)
	theta  <- nconv/(nconv+conv)
	loglik <- dbinom(x=nconv,size=d,p=median(theta),log=TRUE)
	list( pos=pos , nconv=nconv , conv=conv , 
	depth=d , theta=theta , loglik=sum(loglik))
}

slib$join.segments <- function( seg1 , seg2 ){
	pos<-c( seg1$pos , seg2$pos )
	nconv<-c( seg1$nconv , seg2$nconv )
	conv<-c( seg1$conv , seg2$conv )
	d<-c( seg1$d , seg2$d )
	theta<-c( seg1$theta , seg2$theta )
	loglik <- dbinom( x=nconv, size=d, p=median(theta),log=TRUE) 
	list( pos=pos , nconv=nconv , conv=conv , 
	depth=d , theta=theta , loglik=sum(loglik))
}

slib$delta.segments <- function ( seg1, seg2, lambda ){
	#want delta to be positive
	tmp.join<-join.segments( seg1 , seg2 )
	tmp.join$loglik-seg1$loglik-seg2$loglik+lambda
}




while("slib" %in% search())
  detach("slib")
attach(slib)
