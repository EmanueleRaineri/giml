options(digits=4)
source("libgimli.R")
args<-commandArgs(trailingOnly=T)
fname<-args[1]

counts<-read.table(fname,stringsAsFactors=F)
names(counts)<-c("chrom","start","end","nc1","c1","nc2","c2",
	"theta1","theta2","theta12","l1","l2","l12","dl")


#cat(nrow(counts)," rows\n")

for (i in 1:nrow(counts)){

	dl<-counts[i,"dl"]
	lr<-2*(dl)
	pval = 1.0 - pchisq(q=lr,df=1)
	cat(counts[i,1],counts[i,2],counts[i,3],
	counts[i,"nc1"],counts[i,"nc2"],
	counts[i,"theta1"],counts[i,"theta2"],
	counts[i,"theta12"],
	counts[i,"l1"],counts[i,"l2"],counts[i,"l12"],
        pval,sep="\t",fill=F)
	cat("\n")

}
