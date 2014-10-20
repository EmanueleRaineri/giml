library(stats4)

d<-20

diff<-rep(0,1000)
for (i in 1:1000){
	theta<-runif(1)
	nc<-rbinom(100,size=d,p=theta)
	ll<-function(theta){
		r<-dbinom(x=nc,size=d,p=theta)
		-sum(log(r))
	}

	ml<-stats4::mle(minuslog=ll,
		start=list(theta=mean(nc/d)))

	cat(ml@coef,"\t",median(nc/d),"\t",ml@min,"\n")
	diff[[i]]<-ml@coef-median(nc/d)
}
# ml@min is the max likelyhood
