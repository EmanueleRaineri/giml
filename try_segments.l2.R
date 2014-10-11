source("segmentlib.R")

simsize<-100

amrun<-rep(0,simsize)
amrun[1:50]<-runif(n=50,min=0,max=0.4)
amrun[51:70]<-runif(n=20,min=0.5,max=1)
amrun[71:100]<-runif(n=30,min=0,max=0.4)

d<-data.frame(chrom=rep("chr1",simsize),pos=seq(1,simsize,by=1),amrun=amrun)
seg.list<-slib$seg.list.of.data.frame(d)
slib$print.seg.list(seg.list,0)


all.lambda<-c(0.01,0.02,0.05,0.1,0.2,0.5)
slib$loop.over.lambda(seg.list,all.lambda)
