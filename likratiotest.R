options(digits=4)
source("libgimli.R")
args<-commandArgs(trailingOnly=T)
fname<-args[1]

counts<-read.table(fname,stringsAsFactors=F)
names(counts)<-c("chrom","start","end","nc1","c1","nc2","c2")

splitvec<-function(v){
	as.numeric(strsplit(v,",")[[1]])
}

#cat(nrow(counts)," rows\n")

for (i in 1:nrow(counts)){

	nc1<-splitvec(counts[i,"nc1"])
	c1<-splitvec(counts[i,"c1"])
	nc2<-splitvec(counts[i,"nc2"])
	c2<-splitvec(counts[i,"c2"])

	theta1 <- sum(nc1)/(sum(nc1)+sum(c1))
	theta2 <- sum(nc2)/(sum(nc2)+sum(c2))

	theta12 <- (sum(nc1)+sum(nc2))/(sum(nc1)+sum(c1)+sum(nc2)+sum(c2))


	l1<-sum(libgimli$llik(nc1 , c1, theta1))
	l2<-sum(libgimli$llik(nc2 , c2, theta2))
	l12<-sum(libgimli$llik(c(nc1,nc2),c(c1,c2),theta12))

	lr<-2*(l1+l2-l12)
	pval = 1.0 - pchisq(q=lr,df=1)
	cat(counts[i,1],counts[i,2],counts[i,3],length(nc1),length(nc2),
	theta1,theta2,theta12,
        l1,l2,l12,lr,pval,sep="\t",fill=F)
	cat("\n")
}
