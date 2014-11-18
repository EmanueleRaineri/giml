library(lineprof)
#load some real data
source("segmentlib.R")
options(digits=4)

#g199<-read.table("pluto.txt",stringsAsFactors=F,sep=" ")
#g199<-read.table("example.in",stringsAsFactors=F,sep=" ")
g199<-read.table("G199.sample",stringsAsFactors=F,sep=" ")
#g199<-read.table("ziopaperone.txt",stringsAsFactors=F,sep=" ")
seg.list<-slib$seg.list.of.data.frame.lik(g199)
slib$print.seg.list.lik(seg.list,0)
lambda.list<-c(0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000)

segs<-slib$loop.over.lambda.lik(seg.list,lambda.list)
slib$plot.segmentation(g199,segs,1000,762891, 763814)
#print(segs)
lelist<-list()
liklist<-list()
i<-1
for (l in lambda.list){
	m1<-median(segs[segs$lambda==l,2]-segs[segs$lambda==l,1])
	m2<-mean(segs[segs$lambda==l,2]-segs[segs$lambda==l,1])
	cat("lambda:",l," median:",m1," mean:",m2,"\n")
	lelist[[i]]<-segs[segs$lambda==l,2]-segs[segs$lambda==l,1]
	liklist[[i]]<-segs[segs$lambda==l,"loglik"]
	i<-i+1
}
boxplot(lelist,range=0,names=lambda.list)
boxplot(liklist,range=0,names=lambda.list)
