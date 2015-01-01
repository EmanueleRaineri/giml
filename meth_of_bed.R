source("libgimli.R")

bed.fname <- "gencode.v19.TSS.notlow.chr1.promoters"
meth.fname<- "~/Desktop/meth_data/C004GD51_cpg.chr1.txt.gz"

meth<-libgimli$load.meth(meth.fname)
bed<-read.table( bed.fname , stringsAsFactors=F )

