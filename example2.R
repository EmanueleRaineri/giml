source("libgimli.R")
args <- commandArgs(trailingOnly = TRUE)

#datadir<-args[1]
datadir<-"~/Desktop/meth_data"
dmr<-libgimli$load.dmr("G199.G202.100.200.dmr")
print(head(dmr))

xlim<-c(45.179e6,45.20e6)

if (!exists("g199")){
	g199 <- libgimli$load.meth(paste(datadir,"/","G199_cpg.chr1.txt.gz",sep=""))
}
if (!exists("g202")){
	g202 <- libgimli$load.meth(paste(datadir,"/","G202_cpg.chr1.txt.gz",sep=""))
}

postscript("G199.G202.100.200.dmr.eps")
libgimli$plot.intersected.dmr(dmr,xlim[1],xlim[2],"red","blue")
points(g199$pos,g199$meth,xlim=xlim,col=rgb(1,0,0,alpha=0.3),bg=rgb(1,0,0,alpha=0.5),pch=20,cex=0.5)
points(g202$pos,g202$meth,xlim=xlim,col=rgb(0,0,1,alpha=0.3),bg=rgb(0,0,1,alpha=0.5),pch=20,cex=0.5)
dev.off()
