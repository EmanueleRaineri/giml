source("segmentlib.R")

d<-data.frame(chrom=rep("chr1",10),pos=seq(1,10,by=1),amrun=runif(10))
seg.list<-slib$seg.list.of.data.frame(d)
slib$print.seg.list(seg.list)
lambda<-0.2
while(TRUE){
	if (length(seg.list)==1) break
	all.pairs<-slib$all.pairs.l2(seg.list,lambda)
	new.seg.list<-slib$update.segmentation.l2(seg.list,all.pairs)
	slib$print.seg.list(new.seg.list)
	if (length(new.seg.list)==length(seg.list) ) break
	seg.list<-new.seg.list
}
