var.beta<-function(a,b){
	a*b/((a+b)^2*(a+b+1))
}

nc<-rep(28)
co<-rep(4)

make.fn<-function(nc,co){
	fn<-function(theta){
		v<-dbinom(x=nc,size=nc+co,p=theta,log=T)
		sv<-sum(v)
		return(-sv)
	}
	return(fn)
}

fn <- make.fn( nc , co )

opt <- nlm( fn , 0.5 , hessian=T )
opt.sdev <- sqrt( diag( solve( opt$hessian ) ) )

cat ("opt$estimate:", opt$estimate, "\n" )
cat ( "opt.sdev:",opt.sdev, "\n" )
#print ( var.beta( nc+1 , co+1 ) )
