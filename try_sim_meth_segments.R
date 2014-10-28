source("segmentlib.R")

s1<-slib$make.segment(100,3,20)
s2<-slib$make.segment(200,20,2)
s3<-slib$make.segment(300,20,2)
s4<-slib$make.segment(400,3,20)

seg.list<-list(s1,s2,s3,s4)
lambda.list<-c(2000)

segments<-slib$loop.over.lambda.lik(seg.list,lambda.list )

