source("libgimli.R")
dmr<-libgimli$load.dmr("G199.G202.20.200.dmr")

print(head(dmr))

postscript("G199.G202.20.200.dmr.eps")
libgimli$plot.intersected.dmr(dmr,950000,1000000,"red","blue")
dev.off()
