source("libgimli.R")

rdepth<-seq(5,40,by=5)

l1<-c()
l2<-c()
l12<-c()

for (r in rdepth){
	for ( i in 1:100 ) {
		seg1.nc<-rbinom(n=10,size=r,prob=0.8)
		seg1.co<-r-seg1.nc
		seg2.nc<-rbinom(n=10,size=r,prob=0.45)
		seg2.co<-r-seg2.nc
		l1[i]<-libgimli$loglik(seg1.nc,seg1.co)
		l2[i]<-libgimli$loglik(seg2.nc,seg2.co)
		seg12.nc<-c(seg1.nc,seg2.nc)
		seg12.co<-c(seg1.co,seg2.co)
		l12[i]<-libgimli$loglik(seg12.nc,seg12.co)
		cat( r,l1[i],l2[i],l12[i],l12[i]-l1[i]-l2[i],
			sep="\t",fill=TRUE )
	}
}
