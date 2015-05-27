deltacov<-read.table("delta.vs.cov.txt")
postscript("fig_delta_cov.eps")
boxplot(V5~V1,data=deltacov,border="red")
dev.off()
