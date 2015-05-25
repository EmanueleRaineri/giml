gimlifile<-read.table("~/Desktop/meth_data/C001UYA3bs_cpg.chr1.gimli")
colors<-rainbow(4)
postscript("fig_boxplot_size.eps")
boxplot((V3-V2+1)~V12,data=gimlifile,border="red",log="y",
	ylab="length",xlab=expression(lambda),xaxt="n",
	main=expression(paste("length vs ",lambda)))
axis(1,at=1:4,labels=c(1,10,100,1000),cex=0.7)
dev.off()
postscript("fig_boxplot_lik.eps")
boxplot(V10~V12,data=gimlifile,border="red",
	ylab="loglik",xlab=expression(lambda),xaxt="n",
	main=expression(paste("loglik vs ",lambda)))
axis(1,at=1:4,labels=c(1,10,100,1000),cex=0.7)
dev.off()
