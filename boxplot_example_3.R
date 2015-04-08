rpkm.vs.met <- read.table("rpkm.vs.met")

names(rpkm.vs.met)<-c("tss","met","rpkm")

rpkm.low.met<-rpkm.vs.met$rpkm[rpkm.vs.met$met<0.5]
rpkm.high.met<-rpkm.vs.met$rpkm[rpkm.vs.met$met>=0.5]

postscript("boxplot_example_3.eps")
boxplot(rpkm.low.met+1,rpkm.high.met+1,
	range=0,log="y",col=rainbow(2),names=c("low meth","high meth"))
dev.off()
