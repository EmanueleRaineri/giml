var.beta<-function(a,b){
	a*b/((a+b)^2*(a+b+1))
}

#nc<-rep(28)
#co<-rep(4)

nc<-c(5,5)
co<-c(5,5)

make.fn<-function(nc,co){
	fn<-function(theta){
		v<-dbinom(x=nc,size=nc+co,p=theta,log=T)
		sv<-sum(v)
		return(-sv)
	}
	return(fn)
}

fn <- make.fn( nc , co )

make.fn.beta<-function(nc,co){
	est.theta<-nc/(nc+co)
	fn<-function(par){
		a<-par[1]
		b<-par[2]
		v<-log(dbeta(est.theta,a,b))
		sv<-sum(v)
		return (-sv)
	}
	return (fn)
}

mm.beta<-function(m,v){
	a <- ( (1.0 - m ) * m^2 - m * v ) / v;
	b <- a * (1.0 - m) / m ;
	c(a,b)
}

opt <- nlm( fn , 0.5 , hessian=T )
opt.sdev <- sqrt( diag( solve( opt$hessian ) ) )

fn.beta<-make.fn.beta(nc,co)
opt.beta<-nlm( fn.beta, c(nc[1]+1,co[1]+1), hessian=T)
cat("nc,co\n")
str(nc)
str(co)
cat ("opt$estimate:", opt$estimate, "\n" )
cat ( "opt.sdev:",opt.sdev, "\n" )
cat ("opt.beta$estimate:",opt.beta$estimate,"\n")
est.theta<-nc/(nc+co)
mm<-mm.beta(mean(est.theta),var(est.theta))
cat("est.theta:",est.theta,"\n")
cat("mm-a:",mm[1],"\tmm-b:",mm[2],"\n");

make.pab<-function(nc,co){
	pab<-function(a,b){
		v<-gamma(a+b)/(gamma(a)*gamma(b))*(gamma(a+nc)*gamma(b+co))/gamma(a+b+nc+co)
		(a+b)^-0.5*prod(v)
	}
	pab	
}
pab<-make.pab(nc,co)
ma<-0;mb<-0;m<--1;
for(a in 1:50){
	for(b in 1:50){
		if (pab(a,b)>=m){
			ma<-a; mb<-b;
			m<-pab(a,b)
		}
	}
}
cat("ma:",ma," mb:",mb,"\n")
#print ( var.beta( nc+1 , co+1 ) )
