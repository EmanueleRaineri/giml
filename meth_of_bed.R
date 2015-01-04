source("libgimli.R")

bed.fname  <- "C004GDH1_12_Blueprint_release_082014_segments.chr1.active_promoter.bed"
meth.fname <- "~/Desktop/meth_data/C004GD51_cpg.chr1.txt.gz"

if (!exists("meth")){
	meth<-libgimli$load.meth(meth.fname)
	bed<-read.table( bed.fname , stringsAsFactors=F )
}

meth.rd    <- libgimli$rd.of.meth(meth,2,c(6,7))
bed.rd <-  
#lb<-bed[2,2]
#ub<-bed[2,3]
#print(lb)
#print(ub)
#lb<-lb-1000
#ub<-ub+1000


for ( b in 1:nrow(bed) ){
	lb <- bed[b,2]
	ub <- bed[b,3]
	meth.slice <- libgimli$take.slice(meth.rd,lb,ub)
	#"Min."    "1st Qu." "Median"  "Mean"    "3rd Qu." "Max."
	m <- meth.slice$nc/(meth.slice$nc+meth.slice$c)
	s <- summary(m)
	cat(lb,ub,length(m),unname(s["Min."]),unname(s["Median"]),unname(s["Max."]),"\n")
}

#plot(start(meth.slice$ranges),meth.slice$nc/(meth.slice$nc+meth.slice$c),
#	col="red",type='b',ylim=c(0,1))
#abline(v=c(64431711,64432082),col="blue")	
