#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

typedef struct node{
	// from, to refers to indexes for the four arrays mentioned above, 
	// not to real positions along 
	// the chromosome
	int from;
	int to;
	float loglik;
	float delta;
	int segment_id;
	struct node* prev;
	struct node* next;
} node;

float dbinom(int x, int size, float p){
	if (p==1){
		if (x==size) return 0;
		else {
			fprintf(stderr,"invalid dbinom call x:%d size:%d p:%f\n",
			x,size,p);
			exit(1);
		}
	}
	if (p==0){
		if (x==0) return 0;
		else {
			fprintf(stderr,"invalid dbinom call x:%d size:%d p:%f\n",
			x,size,p);
			exit(1);
		}
	}
	float bcoeff = lgammaf(size+1) - lgammaf(x+1) - lgammaf(size-x+1);
	//printf("%f\n",bcoeff);
	return ( bcoeff + x*log(p) + (size-x)*log(1-p) );
}


void print_list(node* head, float lambda){
	
	node* el = head;
	while(el!=NULL){
	
		printf("%d\t",el->segment_id);
		printf("%d\t",el->from);
		printf("%d\t",el->to);
		printf("%.4f\t",el->loglik);
		printf("%.4f\t",el->delta);
		printf("%.4f\n",lambda);

		
		el = el->next;
	}

}

float update_lik(node* el, float* theta, int* nc, int* c){
	float sumtheta=0,avgtheta;
	int i;
	for( i = el->from; i <= el->to; i++ ){
		sumtheta+=theta[i];		
	}
	avgtheta=sumtheta/(el->to - el->from +1);
	//fprintf(stderr,"avgtheta:%g\n",avgtheta);
	el->loglik=0;
	for( i = el->from; i <= el->to; i++ ){
		el->loglik+=dbinom(nc[i],nc[i]+c[i],avgtheta);
	}
	return (avgtheta);
}


float delta_lik( node* el, float* theta, int* nc, int* c ){
	int i;
	float sumtheta=0,avgtheta;
	for( i = el->from; i <= el->next->to; i++ ){
		sumtheta+=theta[i];		
	}
	avgtheta=sumtheta/(el->next->to - el->from +1);
	el->delta=0;
	for( i = el->from; i <= el->next->to; i++ ){
		el->delta+=dbinom(nc[i],nc[i]+c[i],avgtheta);
	}
	return (el->delta);
}

void list_of_file(char* fname,node* head, int nlines, 
	float* theta, int* pos, int* nc, int* c, int* segment_id ){
	FILE* in=fopen(fname,"r");
	const int sizeline=1000;
	char* buffer=malloc(sizeline*sizeof(char));
	node* el = head;
	int i;
	for( i=0; i<nlines; i++ ) {
		buffer =	fgets(buffer,sizeline,in);
		sscanf(buffer,"%*s %d %d %d",pos+i,nc+i,c+i);
		segment_id[i]=i;
		theta[i]=(float)nc[i]/(nc[i]+c[i]);
		el->from = i;
		el->to = i;
		el->loglik = 0 ;
		el->delta=0;
		el->segment_id = i;
		el->next = malloc(sizeof(node));
		el->next->prev = el;
		el = el->next;
	}
	el->prev->next=NULL;
	fclose(in);
}

void print_node(node* el){
		fprintf(stderr,"-------------\n");
		fprintf(stderr,"%d\t",el->segment_id);
		fprintf(stderr,"%d\t",el->from);
		fprintf(stderr,"%d\t",el->to);
		fprintf(stderr,"%.4f\t",el->loglik);
		fprintf(stderr,"%.4f\n",el->delta);
		fprintf(stderr,"-------------\n");
}


int main(int argc, char* argv[]){

	int i,j=0;
	int nlines=0;
	/*count lines*/
	FILE* in = fopen(argv[1],"r");
	const int sizeline=1000;
	char* buffer=malloc(sizeline*sizeof(char));
	while(!feof(in)){
		fgets(buffer,sizeline,in);
		nlines++;
	}
	nlines--;
	fprintf(stderr,"nlines:%d\n",nlines);
	fclose(in);
	/* 5 static data structures */
	float*  theta=malloc(nlines*sizeof(float));
	int*    pos=malloc(nlines*sizeof(int));
	int*    nc=malloc(nlines*sizeof(int));
	int*    c=malloc(nlines*sizeof(int));
	int*    segment_id=malloc(nlines*sizeof(int));
	float maxdelta; 
	float avgtheta;
	float lambda = 1;
	float delta;
	node* maxn;
	node* el, *head;
	head=malloc(sizeof(node));
	head->prev=NULL;
	el=head;
	//
	list_of_file(argv[1],head,nlines,theta,pos,nc,c,segment_id);
	// initialize delta, lik
	maxdelta = -FLT_MAX; 
	el=head;
	while(el!=NULL){
			if (el->next==NULL) {
				el->delta=-1000;
				break;
			}
			avgtheta=update_lik(el,theta,nc,c);
			delta=delta_lik(el,theta,nc,c);
			el = el->next;
	}	
	/* */
	print_list(head,0);	
	while(lambda<10){
		print_list(head,lambda);	
		for( el=head; el!=NULL; el=el->next ){
			if (el->delta > maxdelta) { maxdelta=el->delta; maxn=el; }
		}
		fprintf(stderr,"max delta:%f\n",maxn->delta);
		if ((maxn->delta + lambda) >0){
			/* merge */
			fprintf(stderr,"merging...\n");
			j++;
			node* el2=maxn->next;
			for( i=el2->from; i<=el2->to; i++ ){
				segment_id[i]=maxn->segment_id;
			}
			fprintf(stderr,"maxn before merging\n");
			print_node(maxn);
			maxn->next=el2->next;
			el2->next->prev=maxn;
			maxn->to = el2->to;
			update_lik(maxn,theta,nc,c);
			free(el2);
			/* change local deltas */
			delta_lik(maxn->prev,theta,nc,c);
			delta_lik(maxn,theta,nc,c);
			fprintf(stderr,"maxn after merging\n");
			print_node(maxn);
		} else {
			//increment lambda
			lambda++;
		}	
	}
	fprintf(stderr,"%d merging operation(s)\n",j);
}
