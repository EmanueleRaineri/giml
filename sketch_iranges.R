library(IRanges)
ir<-IRanges::IRanges(1:1000000,1:1000000)

v1<-runif(10)
v2<-runif(10)

rd<-IRanges::RangedData(ir,v1,v2)

query<-IRanges::RangedData(IRanges::IRanges(seq(1,100000,by=5000),seq(1,100000,by=5000)+10000))

ov<-IRanges::findOverlaps(query,rd)
ra<-ranges(ov,ranges(query),ranges(rd))

