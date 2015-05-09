options(digits=5)
var.beta<-function(a,b){
	a*b/((a+b)^2*(a+b+1))
}
#nc<-c(20,40,40)
#co<-c(10,20,20)
#nc<-33
#co<-17
#nc<-20
#co<-10
nc<-100
co<-5

fn<-function(theta){
	v<-dbinom(x=nc,size=nc+co,p=theta,log=T)
	sv<-sum(v)
	return(-sv)
}
#
make.fn<-function(nc,co){
	fn<-function(theta){
		v<-dbinom(x=nc,size=nc+co,p=theta,log=T)
		sv<-sum(v)
		return(-sv)
	}
	return(fn)
}
#
opt<-nlm(fn,0.5,hessian=T)
opt.sdev<-sqrt(diag(solve(opt$hessian)))
#print (opt$estimate)
#print (opt.sdev)
#print (var.beta(nc+1,co+1))
#
cov<-40
cat("#","nc","co","est_theta","ncsum","cosum",
	"est_theta_sum","opt$estimate", "opt.sdev",
		"opt.sum$estimate", "opt.sum.sdev","beta.sum.sdev","\n")
for ( i in 1:1000 ){
	theta<-runif(1)
	nc <- rbinom(n=3,size=cov,prob=theta)
	co <- ( cov - nc )
	fn<-function( theta ){
		v<-dbinom( x=nc, size=nc+co, p=theta, log=T)
		sv<-sum(v)
		return(-sv)
	}
	fn.vanilla<-make.fn(nc,co)
	ncsum<-sum(nc)
	cosum<-sum(co)
	cat(paste(nc,collapse=","),"\t",paste(co,collapse=","),"\t")
	opt<-nlm( f=fn.vanilla, p=0.5, hessian=T )
	if (any(co>0)){
		opt.sdev<-sqrt( diag(solve(opt$hessian)) )
	} else{
		opt.sdev<-0
	}
	fn.sum<-make.fn(ncsum,cosum)
	opt.sum<-nlm(f=fn.sum,p=0.5,hessian=T)
	if (cosum>0){
		opt.sum.sdev<-sqrt(diag(solve(opt.sum$hessian)))
	} else {
		opt.sum.sdev<-0
	}
	cat(paste(nc/(co+nc),collapse=","), ncsum, cosum, ncsum/(cosum+ncsum),
		opt$estimate, opt.sdev,
		opt.sum$estimate, opt.sum.sdev,
		sqrt(var.beta(ncsum+1,cosum+1)), sep="\t",fill=T)
}
