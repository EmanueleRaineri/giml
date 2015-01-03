source("libgimli.R")

datadir<-"/Users/emanueleraineri/Desktop/meth_data"
fgimli<-paste(datadir,"/C004GD51_cpg.chr1.gimli.gz",sep="")

if (!exists("gimli")){
	gimli<-libgimli$load.segments(fgimli)
}

ul<-unique(gimli$lambda)


lul<-length(ul)
colors<-rainbow(lul)
seg<-list()
for ( l in 1:lul ){
	tmp<-gimli[gimli$lambda==ul[l],]
	seg[[l]]<-tmp[,3]-tmp[,2]+1		
}

postscript("boxplot1.eps")

boxplot(seg,range=0,col=colors,log="y",
	ylab="ln(length)",xlab=expression(lambda),xaxt="n",
	main=expression(paste("length vs ",lambda))
	)

axis(1,at=1:lul,labels=ul,cex=0.7)
grid()

dev.off()

seg<-list()
for ( l in 1:lul ){
	tmp<-gimli[gimli$lambda==ul[l],]
	seg[[l]]<-tmp$loglik
}


postscript("boxplot2.eps")

boxplot( seg,range=0,col=colors,
	ylab="loglik",xlab=expression(lambda),xaxt="n",
	main=expression(paste("loglik vs ",lambda))
	)

axis(1,at=1:lul,labels=ul,cex=0.7)
grid()

dev.off()



