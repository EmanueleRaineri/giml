source("libgimli.R")
args <- commandArgs(trailingOnly = TRUE)
datadir<-args[1]

suffix<-"_cpg.chr1.gimli.100.filtered"
fg199<-paste(datadir,"/G199",suffix,sep="")
#fg200<-paste(data,"/G200",suffix,sep="")
#fg201<-paste(data,"/G201",suffix,sep="")
fg202<-paste(datadir,"/G202",suffix,sep="")

g199<-libgimli$load.segments(fg199)
#g200<-libgimli$load.segments(fg200)
#g201<-libgimli$load.segments(fg201)
g202<-libgimli$load.segments(fg202)

if (!exists("g199meth")){
	g199meth <- libgimli$load.meth(paste(datadir,"/","G199_cpg.chr1.txt.gz",sep=""))
}
if (!exists("g202meth")){
	g202meth <- libgimli$load.meth(paste(datadir,"/","G202_cpg.chr1.txt.gz",sep=""))
}

postscript("G199.G200.G201.G202.chr1.gimli.eps")
#libgimli$plot.multiple.segments( list(g199,g200,g201,g202) , 
#	10000 , 2.5e8, rainbow(4) )
xlim<-c(45.17e6 , 45.20e6)
libgimli$plot.multiple.segments( list(g199,g202) , 
	xlim[1] , xlim[2], c("red","blue") )
cat("points\n")
points(g199meth$pos,g199meth$meth,xlim=xlim,col=rgb(1,0,0),bg=rgb(1,0,0),pch=20,cex=0.3)
points(g202meth$pos,g202meth$meth,xlim=xlim,col=rgb(0,0,1),bg=rgb(0,0,1),pch=20,cex=0.3)
grid()
legend("bottomleft",  c("G199","G202"),  
lty=c(1,1), lwd=c(2.5,2.5),ncol=2,col=c("red","blue")) 
dev.off()
