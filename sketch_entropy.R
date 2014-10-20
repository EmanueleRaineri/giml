library(entropy)


nsamples<-1000

ent<-rep(0,1000)
va<-rep(0,1000)

for (i in 1:1000){
	v<-rbeta(50,20,20)
	ent[[i]]<-entropy::entropy.plugin(entropy::discretize(v,numBins=10,r=c(0,1)),unit="log2")
	va[[i]]<-var(v)	
}
fit<-lm(va~ent)
