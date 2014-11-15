#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>


float dbinom(int x, int size, float p){
	if (p==1){
		if (x==size) return 0;
		else return -FLT_MAX ;
	}
	float bcoeff = lgammaf(size+1) - lgammaf(x+1) - lgammaf(size-x+1);
	//printf("%f\n",bcoeff);
	return ( bcoeff + x*log(p) + (size-x)*log(1-p) );
}




int main(int argc, char* argv[]){

	int i,j;
	const int sizeline=1000;
	char* buffer=malloc(sizeline*sizeof(char));
	int nlines=0;
	
	/*count lines*/
	
	FILE* in = fopen(argv[1],"r");
	
	while(!feof(in)){
		fgets(buffer,sizeline,in);
		nlines++;
	}
	nlines--;
	printf("nlines:%d\n",nlines);

	fclose(in);
	/* 5 static data structures */

	float*  theta=malloc(nlines*sizeof(float));
	int*    pos=malloc(nlines*sizeof(int));
	int*    nc=malloc(nlines*sizeof(int));
	int*    c=malloc(nlines*sizeof(int));
	int*    segment_id=malloc(nlines*sizeof(int));


	typedef struct node{
		// from, to refers to indexes for the four arrays mentioned above, 
		// not to real positions along 
		// the chromosome
		int from;
		int to;
		float loglik;
		int segment_id;
		struct node* prev;
		struct node* next;
	} node;


	node* el, *head;
	el=malloc(sizeof(node));
	head=el;
	head->prev=NULL;

	in=fopen(argv[1],"r");
		


	for( i=0; i<nlines; i++ ) {
		buffer =	fgets(buffer,sizeline,in);
		sscanf(buffer,"%*s %d %d %d",pos+i,nc+i,c+i);
		segment_id[i]=i;
		theta[i]=(float)nc[i]/(nc[i]+c[i]);
		el->from=i;
		el->to=i;
		el->loglik = dbinom(nc[i],nc[i]+c[i],theta[i]);  
		el->segment_id=i;
		el->next=malloc(sizeof(node));
		el->next->prev=el;
		el=el->next;
	}
	el->prev->next=NULL;
	fclose(in);
	
	while(1){
		if (head==NULL) break;
	
		printf("%d\t",head->segment_id);
		printf("%d\t",head->from);
		printf("%d\t",head->to);
		printf("%.4f\n",head->loglik);
		
		head = head->next;
	}
}
