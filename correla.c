#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DELTA 100
#define SIZE 2000000
#define HSIZE 10000
#define LSIZE 1000
#define ULIM 1000000

int bin(float x){
	if ( x == ULIM ) return HSIZE;
	return (int)x/DELTA;
}

int main(int argc, char* argv[]){
	fprintf(stderr,"RAND_MAX=%d\n",RAND_MAX);
	unsigned int seed = atoi(argv[1]);
	srand(seed);
	int le=0;
	int i;
	int* pos = malloc(SIZE*sizeof(int));	
	float* met = malloc(SIZE*sizeof(float));
	int* counter= malloc((HSIZE+1)*sizeof(int));	
	float* x = malloc((HSIZE+1)*sizeof(float));
	float* y = malloc((HSIZE+1)*sizeof(float));
	float* xy= malloc((HSIZE+1)*sizeof(float));
	float* x2 = malloc((HSIZE+1)*sizeof(float));
	float* y2 = malloc((HSIZE+1)*sizeof(float));

	for(i=0;i<HSIZE+1;i++){
		x[i]=0;y[i]=0;xy[i]=0;
		x2[i]=0;y2[i]=0;
		counter[i]=0;
	}

	char* buffer=malloc(LSIZE*sizeof(char));
	int lc=0;
	while(1){
		buffer=fgets(buffer,LSIZE,stdin);
		if (feof(stdin)) break;
		if (buffer==NULL) break;
		sscanf(buffer,"%d %f",pos+lc,met+lc);
		lc++;
	}
	fprintf(stderr,"read %d lines\n",lc);
	fprintf(stderr,"first pos %d %f\n",pos[0],met[0]);
	fprintf(stderr,"last pos %d %f\n",pos[lc-1],met[lc-1]);
	int ssize=10000000;
	int b,diff,pos1,pos2,rnd1,rnd2;
	i = 0;
	while( i < ssize ){
		//fprintf(stdout,"%d\n",i);
		rnd1=((float)rand()/RAND_MAX)*(lc-1);
		rnd2=((float)rand()/RAND_MAX)*(lc-1);
		if (rnd1==rnd2){
			fprintf(stderr,"clash!\n");
			continue;
		}
		pos1=pos[rnd1];
		pos2=pos[rnd2];
		if ( pos1 > 250000000 || pos2> 250000000 ) {
			fprintf(stderr,"oops %d %d %d %d %f %f\n",
			rnd1,rnd2,pos1,pos2,met[rnd1],met[rnd2]);
			exit(EXIT_FAILURE);
		}
		diff = abs(pos1-pos2);
		if (diff>ULIM) continue;
		if (diff>pos[lc-1]) {
			fprintf(stderr,"diff too big, %d->%d %d->%d %f->%f\n",
			rnd1,rnd2,pos1,pos2,met[rnd1],met[rnd2]);
			exit(EXIT_FAILURE);
		}
		
		b = bin(diff);
		fprintf(stderr,"%d\t%d\n",diff,b);
		x[b]+=met[rnd1];
		y[b]+=met[rnd2];
		xy[b]+=met[rnd1]*met[rnd2];
		x2[b]+=met[rnd1]*met[rnd1];
		y2[b]+=met[rnd2]*met[rnd2];
		counter[b]+=1;
		i++;
	}
	float sx,sy;
	for( i=0 ; i<20; i++){
		float mux = x[i]/counter[i];
		float muy = y[i]/counter[i];
		float cov = 	xy[i]/counter[i]-mux*muy;
		
		fprintf(stdout,"%d %d\t",i*DELTA,(i+1)*DELTA);
		fprintf(stdout,"%d\t",counter[i]);
		fprintf(stdout,"%f\t",xy[i]/counter[i]);
		fprintf(stdout,"%f\t",mux);
		fprintf(stdout,"%f\t",muy);
		fprintf(stdout,"%f\t",cov);
		sx = sqrt(1.0/(counter[i]-1)*(x2[i]-counter[i]*mux*mux));
		sy = sqrt(1.0/(counter[i]-1)*(y2[i]-counter[i]*muy*muy));
		fprintf(stdout,"%f\n",cov/(sx*sy));
	}
}
