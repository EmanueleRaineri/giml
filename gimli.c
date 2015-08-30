#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <errno.h>
#include <string.h>
#define NDEBUG 1 //this must appear before the inclusion of assert.h
#include <assert.h>
//
#define DEBUG 0
#define D1 100.0
#define D0 10000.0
#define LINE 1000
#define MAXLINES 3000000
#define FILE_END 0
#define CHANGE_CHROM 1
#define MAXNLAMBDA 100

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
	/*	from, to refer to indexes for the four arrays mentioned above, 
		not to real positions along 
		the chromosome 
	*/
	int from;
	int to;
	float loglik;
	float delta;
	float mletheta;
	float mean;
	float var;
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

float rho(int d){
	float dpen;
	dpen = D1/(1+exp(-d/D0))-D1/2;
	/*if ( dpen < 1 ){
		fprintf(stderr,"dpen<1 : %f\n",dpen);
		exit(1);
	}*/
	if ( dpen > FLT_MAX ){
		fprintf(stderr,"dpen>FLT_MAX : %f with d:%d\n",dpen,d);
		exit(1);
	}
	return dpen;
}

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

void print_node(node* el, int* pos){
		fprintf(stderr,"id:%d\t",el->segment_id);
		fprintf(stderr,"from:%d[%d]\t",el->from,pos[el->from]);
		fprintf(stderr,"to:%d[%d]\t",el->to,pos[el->to]);
		fprintf(stderr,"sum_nc:%d\t",el->sum_nc);
		fprintf(stderr,"sum_c:%d\t",el->sum_c);
		fprintf(stderr,"loglik:%.4f\t",el->loglik);
		fprintf(stderr,"delta:%.4f\n",el->delta);
}

int print_list(FILE* stream, node* head, int* pos, float lambda){
	int le = 0;	
	node* el = head;
	while(el!=NULL){
   		print_node(el,pos);
		le++;
		el = el->next;
	}
	return(le);
}

void mean_var(node* el, table* data, float* mean, float* var){
	float mold,mnew,sold,snew,xk;
	int k;
	int n = el->to - el->from +1;
	if (n==1){
		*mean = data->theta[el->from];
		*var=0.0;
		return;
	}
	mold = data->theta[el->from];
	sold=0.0;
	for ( k=1 ; k< n ; k++ ){
		xk = data->theta[el->from+k];
		mnew = mold + ( xk - mold )/( k + 1 );
		snew = sold + ( xk - mold )* ( xk - mnew );
		mold = mnew;
		sold = snew;
	}
	*mean= mnew;
	*var = snew/(n-1);
}

void update_lik(node* el, table* data){
	int sum_nc=0;
	int sum_c=0;
	int i;
	int* nc=data->nc;
	int* c = data->c;
	for( i = el->from; i <= el->to; i++ ){
		sum_nc += nc[i];
		sum_c  += c[i];
	}
	el->sum_nc   = sum_nc;
	el->sum_c    = sum_c;
	el->mletheta = (float)sum_nc/(sum_nc+sum_c);
	if ( fabs(el->mletheta)<FLT_EPSILON && el->sum_nc>0 ) {
		fprintf( stderr, "update_lik:invalid mle %.4f\n", el->mletheta );
		exit( 1 );
	}
	el->loglik=0;
	for( i = el->from; i <= el->to; i++ ){
		el->loglik += dbinom( nc[i] , nc[i]+c[i] , el->mletheta );
	}
	mean_var(el, data, &(el->mean),&(el->var));
}

void delta_lik( node* el, table* data ){
	int i,sum_nc,sum_c;
	float mletheta;
	int* nc=data->nc;
	int* c = data->c;
	int d;
	float dpen;
	if (el->next == NULL) {
		fprintf(stderr,"el->next==NULL in delta_lik (el->delta=%.4f)\n",el->delta);
		el->delta=-FLT_MAX;
		return;
	}
	sum_nc=el->sum_nc+el->next->sum_nc;
	sum_c = el->sum_c+el->next->sum_c;
	mletheta=(float)sum_nc/(sum_c+sum_nc);
	if ( ( fabs(mletheta)<FLT_EPSILON && sum_nc>0) || 
		( fabs(mletheta - 1 )<FLT_EPSILON && sum_c > 0) ) {
		fprintf( stderr, "delta_lik:invalid mle %.4f\n" , mletheta );
		print_node(el,data->pos);
		print_node(el->next,data->pos);
		exit(1);
	}
	el->delta=0;
	for( i = el->from; i <= el->next->to; i++ ){
		el->delta+=dbinom(nc[i],nc[i]+c[i],mletheta);
	}
	el->delta = el->delta - el->loglik - el->next->loglik ;
	if (el->delta>FLT_EPSILON){
		fprintf( stderr, "warning: delta_lik:invalid el->delta %g\n" , 
			el->delta );
		print_node(el,data->pos);
		print_node(el->next,data->pos);
	}
	// here distance penalty 
	d =  data->pos[el->next->from] - data->pos[el->to]; 
	if ( d<=0 ){
		fprintf( stderr , "d<=0 : %d (%d -> %d )\n" , d , data->pos[el->to],
		data->pos[el->next->from] );
		exit(1);
	}
	dpen = rho(d);
	fprintf(stderr,"rho:%d\t%f\t%f\n",d,dpen,el->delta);
	el->delta = (el->delta)-dpen;
}

void fill_node(node* el, table* data, int i){
		/*
			this fills the node when it contains on locus only,
			during the initial upload of the input file
		*/
		el->from = i;
		el->to   = i;
		el->loglik = data->loglik[i];
		el->delta = 0 ;
		el->segment_id = i;	
		el->sum_nc = data->nc[i];
		el->sum_c  = data->c[i];
		el->mletheta = (float)el->sum_nc/(el->sum_nc+el->sum_c);
		el->mean = data->theta[i];
		el->var  = 0.0;
}

int list_of_file( FILE* in, char* buffer , node* head, table* data ,int* status ){
	char* refchr = malloc( LINE*sizeof(char) );
	char* tmp = malloc( LINE*sizeof(char) );
	sscanf( buffer , "%s %*d %*d %*d" , refchr );
	fprintf( stderr , "refchr:%s\n" , refchr );
	node* el = head;
	int i=0;
	while(1){
		if ( i >= MAXLINES){
			fprintf(stderr,"too many positions in chr\n");
			exit(1);
		}
		sscanf( buffer , "%s %*d %*d %*d" , tmp );	
		if( strcmp( tmp , refchr ) !=0 ) {
			*status = CHANGE_CHROM;
			break;
		}
		sscanf(buffer,"%s %d %d %d",data->chrom,data->pos+i,data->nc+i,data->c+i);
		data->segment_id[i] = i;
		data->theta[i]  = (float)data->nc[i]/(data->nc[i]+data->c[i]);
		data->loglik[i] = dbinom(data->nc[i],data->nc[i]+data->c[i],data->theta[i]);
		fill_node(el,data,i);
		el->next = malloc(sizeof(node));
		if (el->next==NULL){
			fprintf(stderr,"out of memory at line %d\n",__LINE__);
			exit(1);
		}
		el->next->prev = el;
		el = el->next;
		i++;
		buffer = fgets(buffer,LINE,in);
		if (feof(in)) {
			*status=FILE_END;
			break;
		}
		if (NULL == buffer){
			fprintf(stderr,"problem in reading input\n");
			exit(1);
		}
	}
	switch(*status){
		case FILE_END:
			if (el->prev==NULL){
				fprintf(stderr,"el->prov is NULL\n");
				exit(1);
			}
			el->prev->next=NULL;
			free(el); el=NULL;
			break;
		case CHANGE_CHROM:
			el->prev->next=NULL;	
			free(el);
			el=NULL;
			break;
		default:
			fprintf(stderr,"illegal status:%d at line %d\n",*status,__LINE__);
			exit(1);
	}
	free(refchr); free(tmp);
	refchr=NULL; tmp=NULL;
	return i;
}

void free_list(node* head){
	node* el;
	for(el=head->next;el!=NULL;el=el->next){
		free(el->prev);
	}
	free(el);
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
		fprintf( stream , "%s\t%d\t%d\t%d\t" ,
		data->chrom,data->pos[el->from],data->pos[el->to],n);
		fprintf(stream,"%.4f\t%.4f\t%.4f\t", 
			mintheta,el->mletheta,maxtheta);
		fprintf(stream,"%.4f\t%.4g\t", 
			el->mean,el->var);
		fprintf(stream,"%.4g\t%.4g\t%.4f\n",
			el->loglik-adjlik, el->delta, lambda );
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
		fprintf(stderr,"new key:%f current key:%f\n",key,h->heap[i]->delta);
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
		i = parent(i);
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
		fprintf(stderr,"heap_delete:invalid del\n");
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
	//(1)max==head
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
	/**  (2)max==node in the middle of long list **/
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
	/** (3)max==last but one node, merging with the last **/
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

void free_all(node* head, table* data, heap* h){
	free_list( head );
	free( data->chrom );      
	free( data->theta );      
	free( data->loglik );     
	free( data->pos );       
	free( data->nc );       
	free( data->c );
	free( data->segment_id );
	free( data );
	free( h->heap );
	h->heap=NULL;
	free( h );
	h=NULL;
}

int main(int argc, char* argv[]){
	fprintf(stderr,"GIMLI -- emanuele.raineri@gmail.com\n");
	fprintf(stderr,"source code timestamp: %s\n", __TIMESTAMP__ );
	fprintf(stderr,"compiled on %s %s\n",__DATE__ , __TIME__ );
	fprintf(stderr,"FLT_EPSILON=%g\n",FLT_EPSILON);
	#if defined(NDEBUG)
	fprintf(stderr,"Assert disabled.\n");
	#else
	fprintf(stderr,"Assert enabled.\n");
	#endif
	int i,mergec;
	int nlines,le;
	char* buffer=malloc(LINE*sizeof(char));
	FILE *in;
	int cutoff=1;

	if ( argc != 3 ){
		fprintf(stderr,"usage: gimli FILE LAMBDAS\n");
		exit(1);
	}

	if ( strcmp( argv[1] , "-" ) == 0 ){
		in=stdin;
	}else{
		in = fopen(argv[1],"r");
		if ( in == NULL ){
			fprintf(stderr,"can't open %s\n" , argv[1] );
			exit( 1 );
		}
	}

	int nlambda=0;

	float *lambda = malloc(MAXNLAMBDA*sizeof(float));
		
	char *c=argv[2];
	char *ref=c;

	while(1){
		if ( *c=='\0' ) {
			lambda[nlambda++]=atof( ref );	
			break;
		}
		if (*c++==':') {
			lambda[ nlambda++ ]=atof( ref );
			ref = c;
		}
	}

	fprintf( stderr , "%d lambda(s)\n", nlambda );
	for( i=0; i<nlambda; i++ ){
		fprintf(stderr, "lambda[%d]=%.4g\n",i,lambda[i]);				
	}

	table* data;
	
	int ilambda=0,status=CHANGE_CHROM ;
	node* el, *head;
	heap* h;	
	node* maxn;
	node* tmpnext;
	
	buffer = fgets(buffer,LINE,in);

	if (buffer==NULL || feof(in)){
		fprintf(stderr,"can't read from input file\n");
		fprintf(stderr,"feof status:%d\n",feof(in));
		exit(1);
	}
	
read:

	data             = malloc(sizeof(table));
	data->chrom      = malloc(LINE);
	data->theta      = malloc(MAXLINES*sizeof(float));
	data->loglik     = malloc(MAXLINES*sizeof(float));
	data->pos        = malloc(MAXLINES*sizeof(int));
	data->nc         = malloc(MAXLINES*sizeof(int));
	data->c          = malloc(MAXLINES*sizeof(int));
	data->segment_id = malloc(MAXLINES*sizeof(int));
	
	head = malloc(sizeof(node));
	head->prev = NULL;

	nlines = list_of_file( in, buffer,  head  , data, &status );
	fprintf(stderr,"nlines:%d\n",nlines);
	
	data->theta      =  realloc( data->theta , nlines*sizeof(float) );
	data->loglik     =  realloc( data->loglik ,  nlines*sizeof(float) );
	data->pos        =  realloc( data->pos , nlines*sizeof(int) );
	data->nc         =  realloc( data->nc , nlines*sizeof(int) );
	data->c          =  realloc( data->c ,  nlines*sizeof(int) );
	data->segment_id =  realloc( data->segment_id , nlines*sizeof(int) );
	
	/* HEAP */
	h = malloc(sizeof(heap));
	h->heap = malloc(nlines*sizeof(node*));
	h->size = 0;
	/* *** */
	
	fprintf(stderr,"initializing nodes...\n");
	el=head;
	/* initialize nodes in the list computing delta L to next node */	
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
	int loopc=0;
	float deltalik,totloglik=0;
	if (DEBUG) print_heap( h );
	mergec = 0;	
	if ( cutoff ){
		fprintf( stderr , "with distance penalty\n" );
	} else{
		fprintf( stderr , "no distance penalty\n" );
	}
	while(1){
		if (DEBUG) print_heap(h);
		maxn = heap_extract_max(h);
		assert (heap_wrong_index(h)==0);
		if (maxn->prev==NULL && maxn->next==NULL) {
			for (i=ilambda; i<nlambda;i++){
				le = print_segmentation(stdout, head , lambda[i] , data, &totloglik );
			}
			break;
		}
		deltalik=maxn->delta + lambda[ilambda];
		loopc++;
		if ( deltalik >=0 ){
			/* merge maxn with the following node */
			fprintf( stderr , "merging %d ( %d -> %d )@ lambda=%.4f maxn->delta=%.4f deltalik=%.4f\n", 
				maxn->segment_id,data->pos[maxn->from],
				data->pos[maxn->next->to], lambda[ilambda], maxn->delta, deltalik );
			if ( DEBUG ) {
				fprintf(stderr,"maxn before merging\n");
				print_node( maxn , data->pos );
			}
			mergec++;
			tmpnext=maxn->next;
			assert(tmpnext != NULL);
			for( i= tmpnext->from ; i <= tmpnext->to ; i++ ){
					data->segment_id[i] = maxn->segment_id;
			}	
			/** 3 cases : (1)max==head **/
			if ( maxn->prev==NULL && tmpnext!=NULL && tmpnext->next!=NULL ){
				assert(maxn==head);
				if (DEBUG) fprintf(stderr,"merging 1\n");
			 	merge1( maxn,  tmpnext,  h,  data);
			} 
			/**  (2)max==node in the middle of long list **/
			if ( maxn->prev!=NULL && tmpnext!=NULL && tmpnext->next!=NULL ){
				if (DEBUG) fprintf(stderr,"merging 2\n");
			 	merge2( maxn,  tmpnext,  h,  data);
				goto finish;	
			}
			/** (3)max==last but one node, merging with the last **/
			if ( tmpnext->next == NULL ){
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
		} else { // can't merge, deltalik<0
			fprintf( stderr , "not merging %d ( %d -> %d )@ lambda=%.4f maxn->delta=%.4f deltalik=%.4f\n", 
				maxn->segment_id,data->pos[maxn->from],
				data->pos[maxn->next->to],lambda[ilambda],maxn->delta, deltalik );
			assert ( heap_wrong_index(h)==0 );
			heap_insert( h , maxn );
			assert ( heap_wrong_index(h)==0 );
			le = print_segmentation(stdout, head , lambda[ilambda] , data, &totloglik );
			fprintf( stderr , "lambda %.4f %d segment(s) total loglik=%.4f\n" , lambda[ilambda], le, totloglik );
			if ( ilambda < ( nlambda - 1 ) ) ilambda++;
			else break;
		}	
	}
	fprintf( stderr , "%d loop(s) %d merging operation(s)\n" , loopc, mergec );
	free_all(head,data,h);
	
	switch(status){
		case FILE_END:
			break;
		case CHANGE_CHROM:
			ilambda = 0;
			goto read;
			break;
		default:
			fprintf(stderr, "illegal status:%d at line %d\n", status, __LINE__ );
			exit(1);
	}
	
	if ( in != stdin ) fclose( in );
	free( buffer );
	return( 0 );
}