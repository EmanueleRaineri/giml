deltacov<-read.table("delta.vs.cov.txt")
postscript("fig_delta_cov.eps")
boxplot(V8~V1,data=deltacov,border="red",xlab="read depth",ylab=expression(paste(Delta,"L")))
dev.off()
