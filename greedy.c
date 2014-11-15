#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>


float dbinom(int x, int size, float p){
	float bcoeff = lgammaf(size+1) - lgammaf(x+1) - lgammaf(size-x+1);
	//printf("%f\n",bcoeff);
	return ( bcoeff + x*log(p) + (size-x)*log(1-p) );
}




int main(){

	const int size=1000000;
	float a[size];
	int i,j;
	float max;
	int maxi;
	const int sizeline=1000;
	char* buffer=malloc(sizeline*sizeof(char));
	int nlines=0;
	while(!feof(stdin)){
		fgets(buffer,sizeline,stdin);
		nlines++;
	}
	nlines--;
	printf("nlines:%d\n",nlines);
	
	/*count lines*/
	

	
	
	/*for(i=0;i<size;i++){
		a[i]=rand()/10;
	}

	for (j=0;j<100;j++){
		max = -FLT_MAX; 
		maxi=-1;
		i=0;
		while(1){
			if (i==size) break;
			if (a[i]>max) {max=a[i]; maxi=i;}
			i++;
		}
	}

	printf( "max:%f at:%d\n", max , maxi );
*/
	for (j=0;j<1000000;j++){
		max=dbinom(4,10,0.5);
		//printf("dbinom(%d,10,0.0)=%f\n",j,dbinom(j,10,0.5));
	}
	

	/*
		data structures 

		float array theta
		int array pos
		int array nc
		int array c
		int segment_id

		list of segments. a node in the list contains:

		from, to ,loglik, segment_id
		prev,next

		from, to refers to indexes for the four arrays mentioned above, not to real positions along 
		the chromosome
	
			
	*/
}
