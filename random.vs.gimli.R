postscript("fig_variance.eps")

gimli<-read.table("./gimli_mean_var_le15.txt")
rnd<-read.table("./random_mean_var_le15.txt")

plot(rnd$V3,rnd$V4,col="green",log="y",xlab="mean methylation",ylab="variance")
points(gimli$V3,gimli$V4,col="red")

dev.off()


