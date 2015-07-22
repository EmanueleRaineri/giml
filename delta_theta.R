nc1<-2
c1<-2
nc2<-6
c2<-15

dseq<-seq(-1,1,by=0.01)
deltapdf<-rep(0,length(dseq))
for (i in 1:length(dseq)){
	d<-dseq[i]
	st<-system.time(deltapdf[i]<-integrate(
		function(x){dbeta(x,nc1+1,c1+1)*dbeta(x+d,nc2+1,c2+1)},0,1)[[1]])
	print(st)
}
