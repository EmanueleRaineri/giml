d0<-10000
d1<-100
postscript("figrho.eps")
x<-seq(0,1e5,by=1000)
plot(x,d1/(1+exp(-x/d0))-d1/2,col="red")
dev.off()
