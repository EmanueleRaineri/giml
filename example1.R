source("libgimli.R")

data1<-libgimli$load.segments("G199.chr1.gimli.1000")
data2<-libgimli$load.segments("G199.chr1.gimli.20")

#libgimli$plot.segments(data1,1.12e8,1.22e8,"red")
#libgimli$plot.segments(data2[1:100,],1.12e8,1.22e8,"green")
libgimli$plot.multiple.segments(list(data1,data2),1e4,2.5e8,rainbow(2))
