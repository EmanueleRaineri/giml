source("libgimli.R")

rdepth<-seq(5,40,by=5)

l1<-c()
l2<-c()
l12<-c()

for (r in rdepth){
	for ( i in 1:100 ) {
		seg1.nc<-rbinom(n=1,size=r,prob=0.8)
		seg1.co<-r-seg1.nc
		seg2.nc<-rbinom(n=1,size=r,prob=0.45)
		seg2.co<-r-seg2.nc
		l1[i]<-libgimli$loglik(seg1.nc,seg1.co)
		l2[i]<-libgimli$loglik(seg2.nc,seg2.co)
		theta.1<-seg1.nc/r
		theta.2<-seg2.nc/r
		theta.merge <- (seg1.nc + seg2.nc)/(2*r)
		l12[i]<-libgimli$generic.loglik( seg1.nc, seg1.co, theta.merge )+
			libgimli$generic.loglik( seg2.nc, seg2.co, theta.merge )
		cat( r, l1[i], l2[i], theta.1, theta.2, theta.merge,
			l12[i],l12[i]-l1[i]-l2[i],
			sep="\t",fill=TRUE )
	}
}
