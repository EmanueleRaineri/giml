#build an heap of 
#(float,int) sorting on the float

parent <- function(i) { floor(i/2) }
left   <- function(i) { 2*i }
right  <- function(i) { 2*i + 1 }

make.heap <- function(heap.size){
	#prepares the data structure, although
	#it can contain an heap which is not already
	#heapyfied
	list(size=heap.size,heap=list())
}

pair.list.of.array<-function(v){
	#this is only to help with testing
	l=list()
	for (i in 1:length(v)){
		l[[i]]<-c(v[i],i)	
	}
	l
}

max.heapify <- function( a , i ){
	# @a heap of (float,int)
	# @i index in the array
	l <- left(i)
	r <- right(i)
	if (l<=a$size && a$heap[[l]][1]> a$heap[[i]][1]){
		largest<-l
	}
	else largest<-i
	if ( r<=a$size && a$heap[[r]][1] > a$heap[[largest]][1] ) {
		largest<-r
	}
	if ( largest != i ){
		tmp<-a$heap[[i]]
		a$heap[[i]]<-a$heap[[largest]]
		a$heap[[largest]]<-tmp
		max.heapify( a , largest )
	}
	else return ( a )
}

build.max.heap <- function( a  ){
	ub<-as.integer(a$size/2)
	for ( i in ub:1 ) {
		a <- max.heapify( a , i )	
	}
	a
}

heap.increase.key <- function( a , i , key ){
	if ( key[1] < a$heap[[i]][1] ) {
		stop("new key smaller than current key")
	}
	a$heap[[i]][1]<-key[1]
	while( i > 1 && a$heap[[parent(i)]][1] < a$heap[[i]][1]  ){
		tmp <- a$heap[[parent(i)]]
		a$heap[[parent(i)]]<-a$heap[[i]]
		a$heap[[i]]<-tmp
		i <- parent( i )
	}
	a
}

max.heap.insert <- function( a, key  ){
	#assumes a max heap to start with
	a$size<-a$size+1
	a$heap[[a$size]] <- c( -Inf , a$size  )
	a <- heap.increase.key( a , a$size , key )
	a
}

heap.extract.max <- function(a){
	m<-a$heap[[1]]		
	a$heap[[1]]<-a$heap[[a$size]]
	a$heap<-a$heap[-a$size]
	a$size<-a$size-1
	a <- max.heapify( a , 1 )
	list(m,a) 
}
