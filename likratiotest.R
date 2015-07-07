source("libgimli.R")
args<-commandArgs(trailingOnly=T)
fname<-args[1]

counts<-read.table(fname)
names(counts)<-c("nc1","c1","nc2","c2")


#iterate on the rows of counts...

nc1<-counts$nc1
c1<-counts$c1
nc2<-counts$nc2
c2<-counts$c2

theta1 <- sum(nc1)/(sum(nc1)+sum(c1))
theta2 <- sum(nc2)/(sum(nc2)+sum(c2))

theta12 <- sum(nc1+nc2)/(sum(nc1)+sum(c1)+sum(nc2)+sum(c2))


l1<-sum(libgimli$llik(nc1 , c1, theta1))
l2<-sum(libgimli$llik(nc2 , c2, theta2))
l12<-sum(libgimli$llik(c(c1,nc2),c(c1,c2),theta12))

lr<-2*(l1+l2-l12)
pval = 1.0 - pchisq(q=lr,df=1)
cat(fname,nrow(counts),theta1,theta2,theta12,
        l1,l2,l12,lr,pval,sep="\t",fill=F)

