#library(IRanges)
source("segmentlib.R")
options(digits=4)
#
cat("defining example data\n")
pos1 <- c(10,20,30,40,50)
d1   <- rep(40,5)
nconv1 <- rbinom(n=5,size=d1,p=0.5)
conv1  <- (d1-nconv1)
seg1  <- slib$make.segment(pos1,nconv1,conv1)
print(seg1)
#
pos2 <- c(60,70,80,90,100)
d2   <- rep(40,5)
nconv2 <- rbinom(n=5,size=d2,p=0.5)
conv2  <- (d2-nconv2)
seg2  <- slib$make.segment(pos2,nconv2,conv2)
print(seg2)
seg.join<-slib$join.segments(seg1,seg2)
print(seg.join)
#
pos3   <- seq( 110 , 150 , by = 10 ) 
d3     <- rep( 40 , 5 )
nconv3 <- rbinom(n=5,size=d3,p=0.5)
conv3  <- (d3-nconv3)
seg3   <- slib$make.segment( pos3 , nconv3 , conv3 )
#print( seg3 )
#
seg.list<-list(seg1,seg2,seg3)
les<-length(seg.list)
seg.pairs<-rep(0,les-1)
for(i in 1:(les-1)){
	seg.pairs[i]<-slib$delta.segments(seg.list[[i]],seg.list[[i+1]],3)
}
#compute new segmentation
idx <- which.max(seg.pairs)
tmp.join <- slib$join.segments(seg.list[[idx]],seg.list[[idx+1]])
seg.list[[idx]]<-tmp.join
seg.list<-seg.list[-(idx+1)]

#load some real data

#g199<-read.table("pluto.txt",stringsAsFactors=F,sep=" ")
#g199<-read.table("G199.sample",stringsAsFactors=F,sep=" ")
#seg.list<-slib$seg.list.of.data.frame.lik(g199)
#slib$print.seg.list.lik(seg.list,0)
#lambda.list<-c(0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000)
#segments<-slib$loop.over.lambda.lik(seg.list,lambda.list)
#slib$plot.segmentation(g199,segments,0.1)
#print(segments)
#lelist<-list()
#liklist<-list()
#i<-1
#for (l in lambda.list){
#	m1<-median(segments[segments$lambda==l,2]-segments[segments$lambda==l,1])
#	m2<-mean(segments[segments$lambda==l,2]-segments[segments$lambda==l,1])
#	cat("lambda:",l," median:",m1," mean:",m2,"\n")
#	lelist[[i]]<-segments[segments$lambda==l,2]-segments[segments$lambda==l,1]
#	liklist[[i]]<-segments[segments$lambda==l,"loglik"]
#	i<-i+1
#}
#boxplot(lelist,range=0,names=lambda.list)
#boxplot(liklist,range=0,names=lambda.list)
