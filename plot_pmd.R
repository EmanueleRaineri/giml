source("libgimli.R")
args<-commandArgs(trailingOnly = TRUE)

meth.counts<-libgimli$load.counts(args[1])

segs<-libgimli$load.segments(args[2])
segs.10<-segs[segs$lambda==10 & segs$start>1e8,]
segs.100<-segs[segs$lambda==100 & segs$start>1e8,]
segs.1000<-segs[segs$lambda==1000 & segs$start>1e8,]


cat("# lambda=10:",nrow(segs.10),"\n")
cat("# lambda=100:",nrow(segs.100),"\n")
cat("# lambda=1000:",nrow(segs.1000),"\n")
lb <- 1.1885e8
ub <- 1.1893e8


theta<-meth.counts$nc/(meth.counts$nc+meth.counts$c)
postscript(args[3])
par(mfrow=c(3,1),mai=c(0.3,0.8,0.3,0.5))
libgimli$plot.segments(segs.10,lb,ub,"red")
points(meth.counts$pos,theta, pch=19,cex=0.2)
text(lb+10000,0.2,expression(paste(lambda,"=10")))
libgimli$plot.segments(segs.100,lb,ub,"green")
points(meth.counts$pos,theta, pch=19,cex=0.2)
text(lb+10000,0.2,expression(paste(lambda,"=100")))
libgimli$plot.segments(segs.1000,lb,ub,"blue")
points(meth.counts$pos,theta, pch=19,cex=0.2)
text(lb+10000,0.2,expression(paste(lambda,"=1000")))
dev.off()
