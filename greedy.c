#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <errno.h>
#define NDEBUG 1
#include <assert.h>

#define DEBUG 0

typedef struct table{
	char* chrom;
	float*  theta;
	float*  loglik;
	int*    pos;
	int*    nc;
	int*    c;
	int*    segment_id;
}table;

typedef struct node{
	// from, to refer to indexes for the four arrays mentioned above, 
	// not to real positions along 
	// the chromosome
	int from;
	int to;
	float loglik;
	float delta;
	float sumtheta;
	float mletheta;
	int sum_nc;
	int sum_c;
	int segment_id;
	int heapidx;
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
	return ( bcoeff + x*log(p) + (size-x)*log(1-p) );
}

int print_list(FILE* stream, node* head, float lambda){
	int le = 0;	
	node* el = head;
	while(el!=NULL){
		fprintf(stream,"%d\t",el->segment_id);
		fprintf(stream, "%d\t",el->from);
		fprintf(stream,"%d\t",el->to);
		fprintf(stream,"%.4f\t",el->loglik);
		fprintf(stream,"%.4f\t",el->delta);
		fprintf(stream,"%.4f\n",lambda);
		le++;
		el = el->next;
	}
	return(le);
}

void update_lik(node* el, table* data){
	float sumtheta=0;
	int sum_nc=0;
	int sum_c=0;
	int i;
	int* nc=data->nc;
	int* c = data->c;
	float* theta = data->theta;
	for( i = el->from; i <= el->to; i++ ){
		sum_nc+=nc[i];
		sum_c+=c[i];
		sumtheta+=theta[i];		
	}
	el->sumtheta = sumtheta;
	el->sum_nc=sum_nc;
	el->sum_c=sum_c;
	el->mletheta=(float)sum_nc/(sum_nc+sum_c);
	if (el->mletheta==0 && el->sum_nc>0) {
		fprintf(stderr,"update_lik:invalid mle %.4f\n",el->mletheta);
		exit(1);
	}
	el->loglik=0;
	for( i = el->from; i <= el->to; i++ ){
		el->loglik+=dbinom( nc[i] , nc[i]+c[i] , el->mletheta );
	}
}

void delta_lik( node* el, table* data ){
	int i,sum_nc,sum_c;
	float mletheta;
	int* nc=data->nc;
	int* c = data->c;
	if (el->next == NULL) {
		fprintf(stderr,"el->next==NULL in delta_lik (el->delta=%.4f)\n",el->delta);
		el->delta=-FLT_MAX;
		return;
	}
	sum_nc=el->sum_nc+el->next->sum_nc;
	sum_c = el->sum_c+el->next->sum_c;
	mletheta=(float)sum_nc/(sum_c+sum_nc);
	if (mletheta==0 && sum_nc>0) {
		fprintf(stderr,"delta_lik:invalid mle %.4f\n",mletheta);
		exit(1);
	}
	el->delta=0;
	for( i = el->from; i <= el->next->to; i++ ){
		el->delta+=dbinom(nc[i],nc[i]+c[i],mletheta);
	}
	el->delta = el->delta - el->loglik - el->next->loglik;
}

int list_of_file(char* fname,node* head, int nlines, table* data ){
	FILE* in=fopen(fname,"r");
	const int sizeline=1000;
	char* buffer=malloc(sizeline*sizeof(char));
	node* el = head;
	int i;
	for( i=0; i<nlines; i++ ) {
		buffer =	fgets(buffer,sizeline,in);
		sscanf(buffer,"%s %d %d %d",data->chrom,data->pos+i,data->nc+i,data->c+i);
		data->segment_id[i]=i;
		data->theta[i]=(float)data->nc[i]/(data->nc[i]+data->c[i]);
		data->loglik[i]=dbinom(data->nc[i],data->nc[i]+data->c[i],data->theta[i]);
		el-> from = i;
		el-> to = i;
		el-> loglik = data->loglik[i];
		el-> delta = 0 ;
		el-> segment_id = i;
		el-> next = malloc(sizeof(node));
		el-> next->prev = el;
		el->sum_nc+=data->nc[i];
		el->sum_c+=data->c[i];
		el->sumtheta+=data->theta[i];
		el->mletheta=(float)el->sum_nc/(el->sum_nc+el->sum_c);
		el = el->next;
	}
	el->prev->next=NULL;
	free(el);
	el=NULL;
	fclose(in);
	return 0;
}

void free_list(node* head){
	node* el;
	for(el=head->next;el!=NULL;el=el->next){
		free(el->prev);
	}
	free(el);
}

void print_node(node* el, int* pos){
		fprintf(stderr,"-------------\n");
		fprintf(stderr,"id:%d\t",el->segment_id);
		fprintf(stderr,"from:%d[%d]\t",el->from,pos[el->from]);
		fprintf(stderr,"to:%d[%d]\t",el->to,pos[el->to]);
		fprintf(stderr,"loglik:%.4f\t",el->loglik);
		fprintf(stderr,"delta:%.4f\n",el->delta);
		fprintf(stderr,"-------------\n");
}

int print_segmentation(FILE* stream, node* head, float lambda, table* data, float* totloglik ){
	node *el;
	int i,n,le=0;
	float mintheta,maxtheta,t,adjlik;
	*totloglik=0;
	for( el=head; el!=NULL; el=el->next ){
		n = el->to-el->from+1;
		mintheta=FLT_MAX;maxtheta=-FLT_MAX;adjlik=0;
		for ( i=el->from; i<=el->to; i++ ){
			t=data->theta[i];
			if (t>maxtheta) maxtheta=t;
			if (t<mintheta) mintheta=t;
			adjlik+=data->loglik[i];
		}
		*totloglik+=el->loglik;
		fprintf(stream,"%s\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4g\t%.4g\t%.4f\n",
		data->chrom,data->pos[el->from],data->pos[el->to],n,
		mintheta,el->mletheta,maxtheta,
		el->loglik-adjlik,el->delta,lambda);
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

unsigned int parent(unsigned int i){
	//return ((i-1)/2);
	return ((i-1)>>1);
}

unsigned int left(unsigned int i){
	return ((i<<1)+1);
}

unsigned int right(unsigned int i){
	return ((i<<1)+2);
}

int heap_wrong_index(heap* h){
	unsigned int i;
	for(i=0;i<h->size;i++){
		if (h->heap[i]->heapidx!=i) {
			fprintf(stderr,"wrong index at %d\n",i);
			return(1);
		}
	}
	return(0);
}

void max_heapify(heap* h, unsigned int i){
	unsigned int l = left(i);
	unsigned int r = right(i);
	unsigned int largest;
	node* tmp;
	if (l < h->size && h->heap[l]->delta > h->heap[i]->delta ) largest = l;
	else largest = i;
	if (r < h->size && h->heap[r]->delta  > h->heap[largest]->delta ) largest=r;
	if ( largest != i ){
		tmp = h->heap[largest];
		h->heap[largest] = h->heap[i];
		h->heap[i] = tmp;
		h->heap[i]->heapidx=i;
		h->heap[largest]->heapidx=largest;
		max_heapify(h,largest);
	}
}

void heap_increase_key(heap* h , unsigned int i , float key ){
	// in practice use this only for de novo insertion
	if (  key < h->heap[i]->delta ){
		fprintf(stderr,"new key smaller than current key\n");
		exit(1);
	}
	if (i<0 || i>=h->size){
		fprintf(stderr,"position %d does not exist\n",i);
		exit(1);
	}
	h->heap[i]->delta=key;
	node* tmp;
	while( i > 0 && h->heap[parent(i)]->delta < h->heap[i]->delta   ){
		tmp = h->heap[i];
		h->heap[i] = h->heap[parent(i)];
		h->heap[parent(i)] = tmp;	
		h->heap[i]->heapidx = i;
		h->heap[parent(i)]->heapidx = parent(i);	
		i=parent(i);
	}
}

int check_heap_integrity(heap* h){
	unsigned int i,p;
	for(i=h->size-1;i>0;i--){
		p=parent(i);
		if (h->heap[i]->delta>h->heap[p]->delta){
			fprintf(stderr,"node %d not good\n",p);
			fprintf(stderr,"parent(%d):%f,left(%d):%f,right(%d):%f\n",
				p,h->heap[p]->delta,left(p),
				h->heap[left(p)]->delta,
				right(p),h->heap[right(p)]->delta);
			return 1;
		}
	}
	return 0;
}

void heap_insert (  heap* h , node* n  ){
	//assumes a max heap to start with
	float key = n->delta;
	h->size++;
	n->delta=-FLT_MAX;
	h->heap[h->size-1] = n;
	h->heap[h->size-1]->heapidx = h->size-1;
	heap_increase_key( h , h->size-1 , key );
}

node* heap_extract_max (heap* h){
	node* maxn = h->heap[0];
	h->heap[0]=h->heap[h->size-1];
	h->heap[0]->heapidx=0;
	h->size--;
	max_heapify(h,0);
	return(maxn);
}

void heap_delete (heap* h, unsigned int i){
	node* del;
	heap_increase_key(h,i,FLT_MAX);
	del=heap_extract_max(h);
	if (del==NULL) {
		fprintf(stderr,"heap_delete:invalied del\n");
		exit(1);
	}
}

void print_heap(heap* h){
	if (h->size==0) {
		fprintf(stderr,"empty heap\n");
		return;	
	}
	unsigned int i;
	for (i=0; i< h->size-1; i++){
		fprintf(stderr,"%d:%d:%.4f ",i,h->heap[i]->heapidx,h->heap[i]->delta);
	}
	fprintf(stderr,"%d:%d:%.4f\n",h->size-1,h->heap[h->size-1]->heapidx,h->heap[h->size-1]->delta);
}

float find_heap_max(heap* h){
	unsigned int i;
	float m=-FLT_MAX;
	for(i=0;i<h->size;i++){
		if (h->heap[i]->delta>m) {
			m=h->heap[i]->delta;
		}
	}
	return(m);
}

void merge1(node* maxn, node* tmpnext, heap* h, table* data){
	if (DEBUG) fprintf(stderr,"merging 1\n");
	heap_delete(h,tmpnext->heapidx);
	maxn->to = tmpnext->to;
	update_lik(maxn, data);
	maxn->next=tmpnext->next;
	delta_lik(maxn,data);
	maxn->next->prev=maxn;
	heap_insert(h,maxn);
	free(tmpnext);
	tmpnext=NULL;
}

void merge2(node* maxn, node* tmpnext, heap* h, table* data){
	heap_delete( h , tmpnext->heapidx );
	heap_delete( h , maxn->prev->heapidx );
	maxn->to = tmpnext->to;
	update_lik( maxn, data );
	update_lik( maxn->prev, data );
	maxn->next = tmpnext->next;
	maxn->next->prev = maxn;
	delta_lik( maxn,data );
	delta_lik( maxn->prev,data );	
	heap_insert( h, maxn );
	heap_insert( h , maxn->prev );
	free( tmpnext );
	tmpnext = NULL;
}

void merge3(node* maxn, node* tmpnext, heap* h, table* data){
	heap_delete( h , tmpnext->heapidx );
	if (maxn->prev != NULL) heap_delete( h , maxn->prev->heapidx );
	maxn->to = tmpnext->to;
	update_lik( maxn, data );
	if (maxn->prev != NULL) update_lik( maxn->prev, data );
	maxn->next=NULL;
	maxn->delta=-FLT_MAX;
	if (maxn->prev != NULL) delta_lik( maxn->prev, data );	
	heap_insert( h, maxn );
	if (maxn->prev != NULL) heap_insert( h , maxn->prev );
	free( tmpnext );
	tmpnext = NULL;
}

int main(int argc, char* argv[]){
	#if defined(NDEBUG)
	fprintf(stderr,"Assert disabled.\n");
	#else
	fprintf(stderr,"Assert enabled.\n");
	#endif
	int i,mergec=0;
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
	table* data;
	data=malloc(sizeof(table));

	data->chrom = malloc(1000);
	data->theta = malloc(nlines*sizeof(float));
	data->loglik = malloc(nlines*sizeof(float));
	data->pos = malloc(nlines*sizeof(int));
	data->nc = malloc(nlines*sizeof(int));
	data->c = malloc(nlines*sizeof(int));
	data->segment_id = malloc(nlines*sizeof(int));
	
	int ilambda = 0 ;
	float lambda[13] = {0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000};
	node* el, *head;
	head = malloc(sizeof(node));
	head->prev = NULL;
	heap* h = malloc(sizeof(heap));
	h->heap = malloc(nlines*sizeof(node*));
	h->size = 0;
	el = head;
	list_of_file( argv[1] , head , nlines , data );
	//initialize  lik
	//fprintf(stderr,"initialize lik...");
	//el=head;
	//while(el!=NULL){
	//	update_lik(el,data);
	//	el = el->next;
	//}
	//fprintf(stderr,"done\n");
	//initialize delta
	fprintf(stderr,"initialize delta...");
	el=head;
	while(1){
			if ( el->next == NULL ) {
				el->delta = -FLT_MAX;
				heap_insert( h , el );
				break;
			}
			delta_lik( el , data );
			heap_insert( h , el );
			assert ( heap_wrong_index(h)==0 );
			el = el->next;
	}
	fprintf(stderr,"done\n");
	le=print_list(stderr,head,0);	
	node* maxn;
	node* tmpnext;
	int loopc=0;
	float totloglik;
	if (DEBUG) print_heap(h);
	while(1 && le>1){
		if (DEBUG) print_heap(h);
		maxn = heap_extract_max(h);
		assert (heap_wrong_index(h)==0);
		if (maxn->prev==NULL && maxn->next==NULL) break;
		if ( (maxn->delta + lambda[ilambda]) >0 ){
			/***** merge ****/
			fprintf( stderr , "merging %d\n", maxn->segment_id );
			if (DEBUG) {
				fprintf(stderr,"maxn before merging\n");
				print_node( maxn , data->pos );
			}
			mergec++;
			tmpnext=maxn->next;
			assert(tmpnext!=NULL);
			for( i= tmpnext->from ; i <= tmpnext->to ; i++ ){
					data->segment_id[i] = maxn->segment_id;
			}	
			/** 3 cases : (1)max==head **/
			if (maxn->prev==NULL && tmpnext!=NULL && tmpnext->next!=NULL ){
				assert(maxn==head);
				if (DEBUG) fprintf(stderr,"merging 1\n");
			 	merge1( maxn,  tmpnext,  h,  data);
			} 
			/**  (2)max==node in the middle of long list **/
			if (maxn->prev!=NULL && tmpnext!=NULL && tmpnext->next!=NULL){
				if (DEBUG) fprintf(stderr,"merging 2\n");
			 	merge2( maxn,  tmpnext,  h,  data);
				goto finish;	
			}
			/** (3)max==last but one node, merging with the last **/
			if (tmpnext->next == NULL){
				if (DEBUG) fprintf(stderr,"merging 3\n");
			 	merge3( maxn,  tmpnext,  h,  data);
			}
			finish:
				assert ( heap_wrong_index(h)==0 );
				if (DEBUG) {
					fprintf( stderr , "maxn after merging\n" );
					print_node( maxn ,data->pos );
				}
				if (DEBUG) print_segmentation(stderr,head,lambda[ilambda],data,&totloglik);
				assert ( heap_wrong_index(h)==0 );
		} else { // no possible mergings
			assert ( heap_wrong_index(h)==0 );
			heap_insert( h , maxn );
			assert ( heap_wrong_index(h)==0 );
			le = print_segmentation(stdout, head , lambda[ilambda] , data, &totloglik );
			fprintf( stderr , "lambda %.4f %d segment(s) total loglik=%.4f\n" , lambda[ilambda], le, totloglik );
			if (ilambda<12)
				ilambda++;
			else break;
		}	
		loopc++;
	}
	fprintf( stderr , "%d loop(s) %d merging operation(s)\n" , loopc, mergec );
	free_list(head);
	return(0);
}
