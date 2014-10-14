source("segmentlib.R")

simsize<-1000

amrun<-rep(0,simsize)
amrun[1:500]<-runif(n=500,min=0,max=0.4)
amrun[501:1000]<-runif(n=500,min=0.5,max=1)

d<-data.frame(chrom=rep("chr1",simsize),pos=seq(1,simsize,by=1),amrun=amrun)
seg.list<-slib$seg.list.of.data.frame(d)
slib$print.seg.list(seg.list,0)


all.lambda<-c(0.01,0.02,0.05,0.1,0.2,0.5)
slib$loop.over.lambda(seg.list,all.lambda)
