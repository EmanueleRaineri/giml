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
	if (el->next==NULL) {
		fprintf(stderr,"el->next==NULL in delta_lik\n");
		exit(1);
	}
	for( i = el->from; i <= el->next->to; i++ ){
		sumtheta+=theta[i];		
	}
	avgtheta=sumtheta/(el->next->to - el->from +1);
	fprintf(stderr,"%d->%d avg theta=%f\n",el->segment_id,el->next->segment_id,avgtheta);
	el->delta=0;
	for( i = el->from; i <= el->next->to; i++ ){
		el->delta+=dbinom(nc[i],nc[i]+c[i],avgtheta);
	}
	fprintf(stderr,"%d->%d total lik=%f\n",el->segment_id,el->next->segment_id,el->delta);
	el->delta=el->delta - el->loglik - el->next->loglik;
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
	float*  theta=malloc(nlines*sizeof(float));
	int*    pos=malloc(nlines*sizeof(int));
	int*    nc=malloc(nlines*sizeof(int));
	int*    c=malloc(nlines*sizeof(int));
	int*    segment_id=malloc(nlines*sizeof(int));
	int ilambda=0;
	float maxdelta; 
	float avgtheta;
	float lambda[13] = {0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000};
	float delta;
	node* maxn;
	node* el, *head;
	head=malloc(sizeof(node));
	head->prev=NULL;
	el=head;
	//
	list_of_file(argv[1],head,nlines,theta,pos,nc,c,segment_id);
	// initialize delta, lik
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
	le=print_list(head,0);	
	while(1){
		maxdelta = -FLT_MAX; 
		for( el=head; el!=NULL; el=el->next ){
			if (el->delta > maxdelta) { maxdelta=el->delta; maxn=el; }
		}
		fprintf(stderr,"max delta:%f\n",maxn->delta);
		if ((maxn->delta + lambda[ilambda]) >0){
			/* merge */
			fprintf(stderr,"merging...\n");
			j++;
			node* el2=maxn->next;
			for( i=el2->from; i<=el2->to; i++ ){
				segment_id[i]=maxn->segment_id;
			}
			//
			fprintf(stderr,"maxn before merging\n");
			print_node(maxn);
			//
			maxn->next=el2->next;
			if (el2->next!=NULL){
				el2->next->prev=maxn;
			}
			maxn->to = el2->to;
			update_lik(maxn,theta,nc,c);
			free(el2);
			/* change local deltas */
			if (maxn->prev != NULL){
				delta_lik(maxn->prev,theta,nc,c);
			}
			if (maxn->next!=NULL){
				delta_lik(maxn,theta,nc,c);
			} else {
				maxn->delta=-1000;
			}
			//
			fprintf( stderr , "maxn after merging\n" );
			print_node( maxn );
		} else {
			//le=print_list(head,lambda);	
			le=print_segmentation(head,lambda[ilambda],pos,nc,c,theta);
			fprintf( stderr , "%d segment(s)\n" , le );
			if (ilambda<12)
				ilambda++;
			else break;
		}	
	}
	fprintf( stderr , "%d merging operation(s)\n" , j );
}
