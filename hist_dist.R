if (!exists("counts")){
	counts<-read.table("~/Desktop/meth_data/C000S5A1bs_cpg.counts.chr1")
}
pos<-counts[,2]

pos.dist <- (pos[2:length(pos)]-pos)[1:length(pos)-1]

rho<-function(d){
	d0<-1e4
	d1<-100
	d1/(1+exp(-d/d0))-d1/2+1
}

rho.on.pos.dist<-mapply(rho,pos.dist)


hist(log(rho.on.pos.dist,10),breaks=50,freq=F)
hist(log(pos.dist,10),freq=F)
