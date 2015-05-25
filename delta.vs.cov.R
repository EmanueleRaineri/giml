deltacov<-read.table("delta.vs.cov.txt")
postscript("fig2.eps")
boxplot(V5~V1,data=deltacov,border="red")
dev.off()
