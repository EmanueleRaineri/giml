source("libgimli.R")
dirprefix<-"/home/emanuele/Desktop/meth_data"
f100<-paste(dirprefix,"/","GRAF.5mc.B.GRAF.5mc.es.gimli.100.dmr",sep="")
dmr.100<-libgimli$load.dmr(f100)
f1000<-"/home/emanuele/Desktop/meth_data/GRAF.5mc.B.GRAF.5mc.es.gimli.1000.dmr"
dmr.1000<-libgimli$load.dmr(f1000)
fstripped1<-paste(dirprefix,"/","GRAF.5mc.B_cpg.stripped.txt.gz",sep="")
fstripped2<-paste(dirprefix,"/","GRAF.5mc.es_cpg.stripped.txt.gz",sep="")
if(!exists("str1")){
	cat("loading stripped 1\n")
	str1<-libgimli$load.stripped(fstripped1)
}
if(!exists("str2")){
	cat("loading stripped 2\n")
	str2<-libgimli$load.stripped(fstripped2)
}


dmr.100.sample<-dmr.100[sample(1:nrow(dmr.100),100),]

for ( i in 1:100 ) {
	cat("generating image ",i,"\n")
	chrom<-dmr.100.sample[i,1]
	lb<-dmr.100.sample[i,2]
	ub<-dmr.100.sample[i,3]
	str1.slice<-str1[str1$chrom==chrom & str1$pos>=lb & str1$pos<= ub,]
	str2.slice<-str2[str2$chrom==chrom & str2$pos>=lb & str2$pos<= ub,]
	ptitle=paste(chrom,lb,ub,sep="_")
	png(paste("/home/emanuele/Desktop/meth_data/",ptitle,".png",sep=""))
	par(mfrow=c(2,1))
	plot(str1.slice$pos,str1.slice$nc/(str1.slice$nc+str1.slice$c),type="b",
		ylim=c(0,1) , main = ptitle)
	lines(str2.slice$pos,str2.slice$nc/(str2.slice$nc+str2.slice$c),type="b",
		col="red" )
	boxplot(list(str1.slice$nc/(str1.slice$nc+str1.slice$c),
	str2.slice$nc/(str2.slice$nc+str2.slice$c))
	,names=c("B","es"),
	border=c("black","red")
	)
	dev.off()
}


