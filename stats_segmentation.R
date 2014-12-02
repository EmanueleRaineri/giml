infile.segs<-"G199_cpg.segmeth.chr10.in.pmd.out"
segments<-read.table(infile.segs,stringsAsFactors=F,sep="\t")
names(segments)<-c("chrom","from","to","ncpg","mintheta","mletheta","maxtheta","loglik","delta","lambda")

ul<-unique(segments[,"lambda"])
lul<-length(ul)
colors<-rainbow(lul)
seg<-list()
for ( l in 1:lul ){
	tmp<-segments[segments$lambda==ul[l],]
	seg[[l]]<-tmp[,3]-tmp[,2]+1		
}

boxplot(seg,range=0,col=colors,log="y",
	ylab="ln(length+1)",xlab=expression(lambda),xaxt="n")



axis(1,at=1:13,labels=ul,cex=0.7)
grid()
