source("libgimli.R")
dirdata<-"/Users/emanueleraineri/Desktop/meth_data"
gfile.10<-"G199_cpg.chr1.gimli.10"

segments.10<-libgimli$load.segments(paste(dirdata,"/",gfile.10,sep=""))
#libgimli$plot.segments.box( segments ,  10000 , 1e6, "red" )
gfile.1000<-"G199_cpg.chr1.gimli.1000"
segments.1000<-libgimli$load.segments(paste(dirdata,"/",gfile.1000,sep=""))
gfile.100<-"G199_cpg.chr1.gimli.100"
segments.100<-libgimli$load.segments(paste(dirdata,"/",gfile.100,sep=""))
deltam.10<-segments.10$max-segments.10$min
deltam.100<-segments.100$max-segments.100$min
deltam.1000<-segments.1000$max-segments.1000$min
boxplot(list(deltam.10,deltam.100,deltam.1000),range=0,col=rainbow(3))
width.10<-segments.10$end-segments.10$start+1
width.100<-segments.100$end-segments.100$start+1
width.1000<-segments.1000$end-segments.1000$start+1
boxplot(list(width.10,width.100,width.1000),range=0,col=rainbow(3),log="y")
dmr1<-libgimli$load.dmr("G199.G202.100.200.dmr")
libgimli$plot.intersected.dmr.box(dmr1,45.179e6,45.20e6,"green","blue")
