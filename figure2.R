source("libgimli.R")
args <- commandArgs(trailingOnly = TRUE)
data<-args[1]

suffix<-"_cpg.chr1.gimli.1000"
fg199<-paste(data,"/G199",suffix,sep="")
fg200<-paste(data,"/G200",suffix,sep="")
fg201<-paste(data,"/G201",suffix,sep="")
fg202<-paste(data,"/G202",suffix,sep="")

g199<-libgimli$load.segments(fg199)
g200<-libgimli$load.segments(fg200)
g201<-libgimli$load.segments(fg201)
g202<-libgimli$load.segments(fg202)


postscript("G199.G200.G201.G202.chr1.gimli.eps")
libgimli$plot.multiple.segments( list(g199,g200,g201,g202) , 
	10000 , 2.5e8, rainbow(4) )
grid()
legend("bottomleft",  c("G202","G200","G201","G202"),  
lty=c(1,1,1,1), lwd=c(2.5,2.5,2.5,2.5),ncol=4,col=rainbow(4)) 
dev.off()
