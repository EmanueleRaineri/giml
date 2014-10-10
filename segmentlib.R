slib = new.env()




#segments do not overlap
slib$make.segment <- function(pos,nconv,conv){
	# pos,nconv,conv are all vectors
	d      <- (nconv+conv)
	theta  <- nconv/(nconv+conv)
	loglik <- log(dbinom(x=nconv,size=d,p=median(theta)))
	list( pos=pos , nconv=nconv , conv=conv , 
	depth=d , theta=theta , loglik=sum(loglik))
}

slib$join.segments <- function( seg1 , seg2 ){
	pos<-c( seg1$pos , seg2$pos )
	nconv<-c( seg1$nconv , seg2$nconv )
	conv<-c( seg1$conv , seg2$conv )
	d<-c( seg1$d , seg2$d )
	theta<-c( seg1$theta , seg2$theta )
	loglik <- log(dbinom( x=nconv, size=d, p=median(theta)) )
	list( pos=pos , nconv=nconv , conv=conv , 
	depth=d , theta=theta , loglik=sum(loglik))
}

slib$delta.segments <- function ( seg1, seg2, lambda ){
	tmp.join<-join.segments( seg1 , seg2 )
	tmp.join$loglik-seg1$loglik-seg2$loglik+lambda
}

while("slib" %in% search())
  detach("slib")
attach(slib)

