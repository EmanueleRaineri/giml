#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <errno.h>

#define DEBUG 0

typedef struct node{
	// from, to refer to indexes for the four arrays mentioned above, 
	// not to real positions along 
	// the chromosome
	int from;
	int to;
	float loglik;
	float delta;
	float sumtheta;
	int segment_id;
	struct node* prev;
	struct node* next;
} node;

typedef struct heap{
	unsigned int size;
	node** heap;
} heap;

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

int print_list(node* head, float lambda){
	
	int le = 0;	
	node* el = head;
	while(el!=NULL){
		printf("%d\t",el->segment_id);
		printf("%d\t",el->from);
		printf("%d\t",el->to);
		printf("%.4f\t",el->loglik);
		printf("%.4f\t",el->delta);
		printf("%.4f\n",lambda);
		le++;
		el = el->next;
	}
	return(le);
}

float update_lik(node* el, float* theta, int* nc, int* c){
	float sumtheta=0,avgtheta;
	int i;
	for( i = el->from; i <= el->to; i++ ){
		sumtheta+=theta[i];		
	}
	el->sumtheta = sumtheta;
	avgtheta = sumtheta/(el->to - el->from +1);
	//fprintf(stderr,"avgtheta:%g\n",avgtheta);
	el->loglik=0;
	for( i = el->from; i <= el->to; i++ ){
		el->loglik+=dbinom(nc[i],nc[i]+c[i],avgtheta);
	}
	return (avgtheta);
}

float delta_lik( node* el, int* pos, float* theta, int* nc, int* c ){
	int i;
	float sumtheta=0,avgtheta;
	if (el->next==NULL) {
		fprintf(stderr,"el->next==NULL in delta_lik\n");
		exit(1);
	}
	//for( i = el->from; i <= el->next->to; i++ ){
	//	sumtheta+=theta[i];		
	//}
	sumtheta=(el->sumtheta)+(el->next->sumtheta);
	avgtheta=sumtheta/(el->next->to - el->from +1);
	//fprintf(stderr,"%d->%d avg theta=%f\n",el->segment_id,el->next->segment_id,avgtheta);
	el->delta=0;
	for( i = el->from; i <= el->next->to; i++ ){
		el->delta+=dbinom(nc[i],nc[i]+c[i],avgtheta);
		//fprintf(stderr,"%d : %f\t",pos[i],el->delta);
	}
	//fprintf(stderr,"\n");
	//fprintf(stderr,"%d->%d total lik=%f\n",el->segment_id,el->next->segment_id,el->delta);
	//fprintf(stderr,"%d->%d  lik1=%f lik2=%f\n",el->segment_id,el->next->segment_id,el->loglik,el->next->loglik);
	el->delta=el->delta - el->loglik - el->next->loglik;
	//fprintf(stderr,"%d->%d delta lik=%f\n",el->segment_id,el->next->segment_id,el->delta);
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
		el->delta = 0 ;
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

int print_segmentation(node* head, float lambda, int* pos, int*nc, int* c, float* theta){
	node *el;
	int i,n,le=0;
	float mintheta,maxtheta,avgtheta,t;
	for( el=head; el!=NULL; el=el->next ){
		n = el->to-el->from+1;
		avgtheta=0;mintheta=FLT_MAX;maxtheta=-FLT_MAX;
		for ( i=el->from; i<=el->to; i++ ){
			t=theta[i];
			avgtheta+=t;
			if (t>maxtheta) maxtheta=t;
			if (t<mintheta) mintheta=t;
		}
		avgtheta=avgtheta/n;
		printf("%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4g\t%.4g\t%.4f\n",
		pos[el->from],pos[el->to],n,mintheta,avgtheta,maxtheta,el->loglik,el->delta,lambda);
		le++;
	}
	return(le);
}

node* find_max(node* head){
	float maxdelta; 
	node* maxn;
	node* el;
		maxdelta = -FLT_MAX; 
		for( el=head; el!=NULL; el=el->next ){
			if (el->delta > maxdelta) { maxdelta=el->delta; maxn=el; }
		}
		if (DEBUG) fprintf(stderr,"max delta:%f\n",maxn->delta);
	return(maxn);
}


int parent(int i){
	return ((i-1)/2);
}

int left(int i){
	return (2*i+1);
}

int right(int i){
	return (2*i+2);
}

void max_heapify(heap* h, int i){
	int l = left(i);
	int r = right(i);
	int largest;
	node* tmp;
	if (l < h->size && h->heap[l]->delta > h->heap[i]->delta ) largest = l;
	else largest = i;
	if (r < h->size && h->heap[r]->delta  > h->heap[largest]->delta ) largest=r;
	if ( largest != i ){
		tmp = h->heap[largest];
		h->heap[largest] = h->heap[i];
		h->heap[i] = tmp;
		max_heapify(h,largest);
	}
}

void heap_increase_key(heap* h , int i , float key ){
	// in practice use this only for de novo insertion
	if (  key < h->heap[i]->delta ){
		fprintf(stderr,"new key smaller than current key\n");
		exit(1);
	}
	h->heap[i]->delta=key;
	node* tmp;
	while( i>0 && h->heap[parent(i)]->delta < h->heap[i]->delta   ){
		tmp=h->heap[i];
		h->heap[i]=h->heap[parent(i)];
		h->heap[parent(i)]=tmp;		
		i=parent(i);
	}
}

void max_heap_insert (  heap* h , node* n  ){
	//assumes a max heap to start with
	float key = n->delta;
	h->size++;
	n->delta=-FLT_MAX;
	h->heap[h->size] = n;
	heap_increase_key( h , h->size , key );
}

void max_heap_delete (heap* h, int i){
	h->heap[i]=h->heap[h->size-1];
	h->size--;
	max_heapify(h,i);
}


void print_heap(heap* h){
	int i;
	for (i=0;i<h->size-1;i++){
		fprintf(stderr,"%d:%.4f ",i,h->heap[i]->delta);
	}
	fprintf(stderr,"%d:%.4f\n",h->size-1,h->heap[h->size-1]->delta);
}

node* heap_extract_max (heap* h){
	node* maxn = h->heap[0];
	h->heap[0]=h->heap[h->size-1];
	h->size--;
	max_heapify(h,0);
	return(maxn);
}

int main(int argc, char* argv[]){
	int i,j=0;
	int nlines=0,le;
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
	float*  theta = malloc(nlines*sizeof(float));
	int*    pos = malloc(nlines*sizeof(int));
	int*    nc = malloc(nlines*sizeof(int));
	int*    c = malloc(nlines*sizeof(int));
	int*    segment_id = malloc(nlines*sizeof(int));
	int ilambda = 0 ;
	float avgtheta;
	float lambda[13] = {0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000};
	float delta;
	node* el, *head;
	head=malloc(sizeof(node));
	head->prev=NULL;
	heap* h=malloc(sizeof(heap));
	h->heap=malloc(nlines*sizeof(node*));
	h->size=0;
	el=head;
	//
	list_of_file(argv[1],head,nlines,theta,pos,nc,c,segment_id);
	// initialize  lik
	el=head;
	while(el!=NULL){
			//if (el->next==NULL) {
			//	el->delta=-1000;
			//	break;
			//}
			avgtheta=update_lik(el,theta,nc,c);
			//delta=delta_lik(el,pos,theta,nc,c);
			el = el->next;
	}
	//initialize delta
	el=head;
	while(1){
			if (el->next==NULL) {
				el->delta=-1000;
				break;
			}
			delta=delta_lik(el,pos,theta,nc,c);
			max_heap_insert(h,el);
			el = el->next;
	}
	/* */
	le=print_list(head,0);	
	node* maxn;
	node* maxn2;
	int loopc=0;
	print_heap(h);
	while(1){
		//maxdelta = -FLT_MAX; 
		//for( el=head; el!=NULL; el=el->next ){
		//	if (el->delta > maxdelta) { maxdelta=el->delta; maxn=el; }
		//}
		//if (DEBUG) fprintf(stderr,"max delta:%f\n",maxn->delta);
		maxn = find_max(head);
		if (loopc==0){
			maxn2 = heap_extract_max(h);
			fprintf(stderr,"%.4f\t%.4f\n",maxn->loglik,maxn2->loglik);
			print_heap(h);
		}
		if ((maxn->delta + lambda[ilambda]) >0){
			/* merge */
			fprintf(stderr,"merging...\n");
			j++;
			node* el2=maxn->next;
			for( i= el2->from; i <= el2->to; i++ ){
				segment_id[i]=maxn->segment_id;
			}
			//
			if (DEBUG) {
				fprintf(stderr,"maxn before merging\n");
				print_node(maxn);
			}
			//
			maxn->next=el2->next;
			if (el2->next!=NULL){
				el2->next->prev=maxn;
			}
			maxn->to = el2->to;
			free(el2);
			//
			update_lik( maxn , theta , nc , c );
			/* change local deltas */
			if (maxn->prev != NULL){
				delta_lik( maxn->prev , pos , theta , nc , c );
			}
			if (maxn->next!=NULL){
				delta_lik( maxn , pos , theta , nc , c );
			} else {
				maxn->delta=-1000;
			}
			//
			if (DEBUG) {
				fprintf( stderr , "maxn after merging\n" );
				print_node( maxn );
			}
		} else {
			//le=print_list(head,lambda);	
			le=print_segmentation(head,lambda[ilambda],pos,nc,c,theta);
			fprintf( stderr , "%d segment(s)\n" , le );
			if (ilambda<12)
				ilambda++;
			else break;
		}	
		loopc++;
	}
	//print_heap(h);
	fprintf( stderr , "%d loops %d merging operation(s)\n" , loopc, j );
	return(0);
}
