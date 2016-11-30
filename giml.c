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
	float deltalik;
	float dpen;
	float psi;
	float mletheta;
	float mlese;
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
	if ( dpen > FLT_MAX ){
		fprintf(stderr,"dpen>FLT_MAX : %f with d:%d\n",dpen,d);
		exit(1);
	}
	return dpen;
}

float dist_pen(node* el, table* data){
	int x;
	if (el->next==NULL){
		return INFINITY;
	}
	x=data->pos[el->next->from] - data->pos[el->to]; 
	if ( x<=0 ){
		fprintf( stderr , "x<=0 : %d (%d -> %d )\n" , x , data->pos[el->to],
		data->pos[el->next->from] );
		exit(1);
	}
	return(rho(x));
}

float logbinom(int x, int size, float p){
	if (p==1){
		if (x==size) return 0;
		else {
			fprintf(stderr,"invalid logbinom call x:%d size:%d p:%f\n",
			x,size,p);
			exit(1);
		}
	}
	if (p==0){
		if (x==0) return 0;
		else {
			fprintf(stderr,"invalid logbinom call x:%d size:%d p:%f\n",
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
		fprintf(stderr,"delta:%.4f\n",el->psi);
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

/*void mean_var(node* el, table* data, float* mean, float* var){
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
}*/

float sum_log_lik(int from, int to, float theta, table* data){
	int i;
	float sum=0;
	for(i=from;i<=to;i++){
		sum += logbinom(data->nc[i],data->nc[i]+data->c[i],theta);
	}
	return sum;
}

float node_lik(node* el, table* data){
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
	/*el->mlese = sqrt(el->mletheta*(1-el->mletheta)/(sum_nc+sum_c));*/
	if ( fabs(el->mletheta)<FLT_EPSILON && el->sum_nc>0 ) {
		fprintf( stderr, "node_lik:invalid mle %.4f\n", el->mletheta );
		exit( 1 );
	}
	//mean_var(el, data, &(el->mean),&(el->var));
	return (sum_log_lik(el->from,el->to,el->mletheta,data));
}

float lik_merge( node* el, table* data ){
	/* likelihood of 2 merged nodes */
	if (el->next == NULL) {
		fprintf(stderr,
			"el->next==NULL in delta_lik (el->psi=%.4f)\n",
			el->psi);
		return -INFINITY;
	}
	float mletheta;
	float lm=0;
	int sum_nc,sum_c;
	sum_nc=el->sum_nc+el->next->sum_nc;
	sum_c = el->sum_c+el->next->sum_c;
	mletheta=(float)sum_nc/(sum_c+sum_nc);
	lm=sum_log_lik(el->from,el->next->to,mletheta,data);
	return lm;
}

float delta_lik( node* el, table* data ){
	int sum_nc,sum_c;
	float mletheta;
	float dlik;
	if (el->next == NULL) {
		fprintf(stderr,"el->next==NULL in delta_lik (el->psi=%.4f)\n",el->psi);
		return -INFINITY;
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
	dlik = lik_merge( el , data ) - el->loglik - el->next->loglik ;
	if (dlik > FLT_EPSILON){
		fprintf( stderr, 
		"warning: delta_lik:invalid dlik (%g)-(%g)-(%g)=%g\n" , 
		 lik_merge(el,data), el->loglik, el->next->loglik, dlik );
		print_node(el,data->pos);
		print_node(el->next,data->pos);
	}
	return dlik;
}

void set_psi(node* el, table * data){
	if (el->next==NULL) {
		el->deltalik=-INFINITY;
		el->dpen=INFINITY;
		el->psi=-INFINITY;
		return;
	}
	el->dpen=dist_pen( el,  data);
	el->deltalik=delta_lik( el,  data);
	el->psi=( el->deltalik ) - ( el->dpen );
}

node* init_node(node* el, table* data, int i){
		/*
			this fills the node when it contains on locus only,
			during the initial upload of the input file
		*/
		el->from = i;
		el->to   = i;
		el->loglik= logbinom(data->nc[i],data->nc[i]+data->c[i],data->theta[i]);
		el->psi = 0 ;
		el->segment_id = i;	
		el->sum_nc = data->nc[i];
		el->sum_c  = data->c[i];
		el->mletheta = (float)el->sum_nc/(el->sum_nc+el->sum_c);
		el->mlese=0;
		el->mean = data->theta[i];
		el->var  = 0.0;
		return(el);
}

int table_of_file(FILE* in, char* buffer ,table* data ,int* status){
	char* refchr = malloc( LINE*sizeof(char) );
	char* tmp = malloc( LINE*sizeof(char) );
	sscanf( buffer , "%s %*d %*d %*d" , refchr );
	fprintf( stderr , "refchr:%s\n" , refchr );
	int i=0;
	while(1){
		if ( i >= MAXLINES){
			fprintf(stderr,"too many positions in chr\n");
			exit(1);
		}
		sscanf( buffer , "%s %*d %*d %*d" , tmp );	
		if( strcmp( tmp , refchr ) !=0 ) { *status = CHANGE_CHROM; break; }
		sscanf(buffer,"%s %d %d %d",data->chrom,data->pos+i,data->nc+i,data->c+i);
		data->segment_id[i] = i;
		data->theta[i]  = (float)data->nc[i]/(data->nc[i]+data->c[i]);
		if (data->nc[i]+data->c[i] > 0) {i++;}
		else {fprintf(stderr,"skipping (%s\t%d)\n",data->chrom,data->pos[i]);}
		buffer = fgets(buffer,LINE,in);
		if (feof(in)) { *status=FILE_END; break; }
		if (NULL == buffer){
			fprintf(stderr,"problem in reading input\n");
			exit(1);
		}
	}
	free(refchr); free(tmp);
	refchr=NULL; tmp=NULL;
	return i;
}

node* list_of_table( node * head, table* data, int le ){
	node* el = head;
	int i;
	for (i=0;i<le;i++){
		el = init_node(el,data,i);
		el->next = malloc(sizeof(node));
		if (el->next==NULL){
			fprintf(stderr,"out of memory at line %d\n",__LINE__);
			exit(1);
		}
		el->next->prev = el;
		el = el->next;
	}
	el->prev->next=NULL;
	free(el);
	return head;
}

void free_list(node* head){
	node* el;
	for(el=head->next;el!=NULL;el=el->next){
		free(el->prev);
	}
	free(el);
}

void prettyprint(FILE* stream, int* v, int start, int end, char sep){
	int i;
	for ( i=start; i<end ; i++ ){
		fprintf(stream,"%d%c",v[i],sep);			
	}
	fprintf(stream,"%d",v[end]);			
}

int print_segmentation(FILE* stream, node* head, float lambda, table* data ){
	node *el;
	int le=0;
	for( el=head; el!=NULL; el=el->next ){
		le++;
		fprintf( stream , "%s\t%d\t%d\t%g\t%d\t%d\t%f" ,
		data->chrom, data->pos[el->from], data->pos[el->to], 
		lambda, el->from, el->to, el->mletheta );
		fprintf( stream, "\n" );
	}
	return(le);
}

node* find_max(node* head){
	float maxdelta; 
	node* maxn;
	node* el;
	maxdelta = -INFINITY; 
	for( el=head; el!=NULL; el=el->next ){
		if (el->psi > maxdelta) { maxdelta=el->psi; maxn=el; }
	}
	if (DEBUG) fprintf(stderr,"max delta:%f\n",maxn->psi);
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
	if (l < h->size && h->heap[l]->psi > h->heap[i]->psi ) largest = l;
	else largest = i;
	if (r < h->size && h->heap[r]->psi  > h->heap[largest]->psi ) largest=r;
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
	if (  key < h->heap[i]->psi ){
		fprintf(stderr,"new key smaller than current key\n");
		fprintf(stderr,"new key:%f current key:%f\n",key,h->heap[i]->psi);
		exit(1);
	}
	if (i<0 || i>=h->size){
		fprintf(stderr,"position %d does not exist\n",i);
		exit(1);
	}
	h->heap[i]->psi=key;
	node* tmp;
	while( i > 0 && h->heap[parent(i)]->psi < h->heap[i]->psi   ){
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
		if (h->heap[i]->psi>h->heap[p]->psi){
			fprintf(stderr,"node %d not good\n",p);
			fprintf(stderr,"parent(%d):%f,left(%d):%f,right(%d):%f\n",
				p,h->heap[p]->psi,left(p),
				h->heap[left(p)]->psi,
				right(p),h->heap[right(p)]->psi);
			return 1;
		}
	}
	return 0;
}

void heap_insert (  heap* h , node* n  ){
	//assumes a max heap to start with
	float key = n->psi;
	h->size++;
	n->psi=-INFINITY;
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
	heap_increase_key(h,i,INFINITY);
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
		fprintf(stderr,"%d:%d:%.4f ",i,h->heap[i]->heapidx,h->heap[i]->psi);
	}
	fprintf(stderr,"%d:%d:%.4f\n",h->size-1,h->heap[h->size-1]->heapidx,
		h->heap[h->size-1]->psi);
}

float find_heap_max(heap* h){
	unsigned int i;
	float m=-INFINITY;
	for(i=0;i<h->size;i++){
		if (h->heap[i]->psi>m) {
			m=h->heap[i]->psi;
		}
	}
	return(m);
}

void merge(node* maxn, node* tmpnext, heap* h, table* data){
	heap_delete(h,tmpnext->heapidx);
	if ( maxn->prev != NULL ) heap_delete( h , maxn->prev->heapidx );
	maxn->to = tmpnext->to;
	maxn->loglik=node_lik(maxn, data);
	maxn->next=tmpnext->next;
	if (maxn->next!=NULL) maxn->next->prev = maxn;
	set_psi(maxn,data);
	heap_insert( h, maxn );
	if ( maxn->prev != NULL )  {
		set_psi(maxn->prev,data);
		heap_insert( h , maxn->prev );
	}
	free( tmpnext );
}

void free_all(node* head, table* data, heap* h){
	free_list( head );
	free( data->chrom );      
	free( data->theta );      
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
	data->pos        = malloc(MAXLINES*sizeof(int));
	data->nc         = malloc(MAXLINES*sizeof(int));
	data->c          = malloc(MAXLINES*sizeof(int));
	data->segment_id = malloc(MAXLINES*sizeof(int));
	head = malloc(sizeof(node));
	head->prev = NULL;
	nlines = table_of_file(in, buffer,  data, &status );
	fprintf(stderr,"nlines:%d\n",nlines);
	head = list_of_table( head, data, nlines);	
	data->theta      =  realloc( data->theta , nlines*sizeof(float) );
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
		set_psi(el ,data);
		heap_insert( h , el );
		if ( el->next == NULL ) break;
		el = el->next;
	}
	fprintf(stderr,"done\n");
	int loopc=0;
	float curlambda,gain;
	mergec = 0;	
	while(1){
		loopc++;
		maxn = heap_extract_max(h);
		curlambda = lambda[ilambda];
		gain = maxn->psi + curlambda;
		fprintf( stderr,
		"%s\t%d:%d->%d\t%f\t%f\t%f\t%g\t%f\t%s\n",
			data->chrom,
			maxn->segment_id, maxn->from, maxn->to,
			maxn->deltalik, maxn->dpen, maxn->psi,
			curlambda, gain, gain>0?"m":"o" ); 
		if ( gain >=0 ){
			/* merge maxn with the following node */
			mergec++;
			tmpnext = maxn->next;
			assert( tmpnext != NULL );
			for( i= tmpnext->from ; i <= tmpnext->to ; i++ ){
					data->segment_id[i] = maxn->segment_id;
			}
			merge( maxn, tmpnext, h, data );
		} else { // can't merge, deltalik<0
			heap_insert( h , maxn );
			le = print_segmentation( stdout, head , curlambda , data );
			fprintf( stderr , "@lambda=%f, %d segment(s)\n", curlambda, le );
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
			fprintf(stderr,"change of chromosome\n");
			ilambda = 0;
			goto read;
			break;
		default:
			fprintf( stderr, 
			"illegal status:%d at line %d\n", status, __LINE__ );
			exit( 1 );
	}
	if ( in != stdin ) fclose( in );
	free( buffer );
	return( 0 );
}
