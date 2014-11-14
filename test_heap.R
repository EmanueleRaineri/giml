source("heap.R")

#
size<-as.integer(10^3)
a<-lapply(1:size,function(i) {c(runif(1),i)})
h<-make.heap(size)
h$heap<-a
bmark<-system.time(s<-build.max.heap(h))
bmark.ref<-system.time(control<-max(as.numeric(lapply(a,function(e){e[[1]]}))))
print(control)
print(s$heap[[1]])
print(bmark)
print(bmark.ref)
#

