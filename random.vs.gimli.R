fnames<-commandArgs(trailingOnly = T )
gimli<-read.table(fnames[1])
rnd<-read.table(fnames[2])
postscript(fnames[3])

plot(rnd$V3,rnd$V4,col="green",log="y",xlab="mean methylation",ylab="variance")
points(gimli$V3,gimli$V4,col="red")

dev.off()


