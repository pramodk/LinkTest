/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2011                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include <sion.h>

#define MAXN  300000

int getrandomize_vector_search (int n, int *vec);
int getrandomize_vector_swap (int n, int *vec);
int getpartner (int penr, int n, int *vec, int *pepartner);
int getpartneropt (int penr, int n, int *vec, int *pepartner);

int main(int argc, char **argv) {

  int *vec1, *vec2, *pepartner, *pepartneropt;
  int r,i,n,countequal,countdiff, part;
  double starttime,endtime;

  vec1 = (int *) malloc(sizeof(int)*MAXN);
  if (vec1 == NULL) {
    fprintf(stderr, "cannot allocate vec[%lu] ...\n", (unsigned long) sizeof(int)*MAXN);
    exit(0);
  }
  vec2 = (int *) malloc(sizeof(int)*MAXN);
  if (vec1 == NULL) {
    fprintf(stderr, "cannot allocate vec[%lu] ...\n", (unsigned long) sizeof(int)*MAXN);
    exit(0);
  }

  pepartner = (int *) malloc(sizeof(int)*MAXN);
  if (pepartner == NULL) {
    fprintf(stderr, "cannot allocate pepartner[%lu] ...\n", (unsigned long) sizeof(int)*MAXN);
    exit(0);
  }

  pepartneropt = (int *) malloc(sizeof(int)*MAXN);
  if (pepartneropt == NULL) {
    fprintf(stderr, "cannot allocate pepartneropt[%lu] ...\n", (unsigned long) sizeof(int)*MAXN);
    exit(0);
  }

  n=MAXN;

  for(r=1;r<=72;r++) {
    printf("\n");
    n=r*4*1024; 
    /* n=r*16;  */
    starttime=_sion_get_time();
    getrandomize_vector_swap(n,vec1);
    endtime=_sion_get_time();
    countequal=0;  for(i=0;i<n;i++) { if(vec1[i]==i) countequal++;};
    printf("swap:       %10d %10.4fs (%d equals)\n",n,endtime-starttime,countequal);


    /* for(i=0;i<n;i++) vec1[i]=i;  */

    for(part=8;part<=8;part++) {
      for(i=0;i<n;i++) vec2[i]=vec1[i];
      starttime=_sion_get_time();
      getpartner(part,n,vec1,pepartner);
      endtime=_sion_get_time();
      printf("partner:    %10d %10.4fs\n",n,endtime-starttime);

      starttime=_sion_get_time();
      getpartneropt(part,n,vec2,pepartneropt);
      endtime=_sion_get_time();
      countdiff=0;  for(i=0;i<n;i++) { if(pepartner[i]!=pepartneropt[i]) countdiff++;};
      printf("partneropt: %10d %10.4fs  (%d diffs)\n",n,endtime-starttime,countdiff);
    }
  }

  if(0) {
    for(r=1;r<=72;r++) {
      n=r*4*1024;
      starttime=_sion_get_time();
      getrandomize_vector_search(n,vec1);
      endtime=_sion_get_time();
      countequal=0;  for(i=0;i<n;i++) { if(vec1[i]==i) countequal++;};
      printf("search: %10d %10.4fs (%d equals)\n",n,endtime-starttime,countequal);
    
    }
  }

#ifdef OUTPUTRANDOMDIST
  if(my_rank==0) {
    printf("Random distribution: "); 
    for(i=0;i<n;i++) { printf("%d,",vec[i]); };
    printf("\n");
  }
#endif


}


int getrandomize_vector_search (int n, int *vec) {
  int *helpline=NULL;
  int i;
  int pos;

  helpline = (int *) malloc(sizeof(int)*n);
  if (helpline == NULL) {
    fprintf(stderr, "cannot allocate helpline[%lu] ...\n", (unsigned long) sizeof(int)*n);
    exit(0);
  }
  for(i=0;i<n;i++) {helpline[i]=i;}

  srand(1);
  for(i=0;i<n;i++) {
    do {
      pos= (int) ((double) rand() * (double) (n) / (double) RAND_MAX );
      /*           printf("%03d: pos=%3d -->%3d\n",i,pos,helpline[pos]);  */
    } while((helpline[pos]<0 || (pos>=n)));
    vec[i]=pos;
    helpline[pos]=-1;
    /*         printf("line[%03d]=%03d (%d)\n",i,pos,helpline[pos]);  */
  }
  free(helpline);
  return(1);
}

int getrandomize_vector_swap (int n, int *vec) {

  int help,pos;
  int i;
  for(i=0;i<n;i++) {vec[i]=i;}
  for(i=0;i<n;i++) {
    pos= (int) ((double) rand() * (double) (n) / (double) RAND_MAX );
    /* 	printf(" exch i=%3d <-> pos=%3d\n",i,pos);  */
    help=vec[i];
    vec[i]=vec[pos];
    vec[pos]=help;
  }
  return(1);
}

int getpartner (int penr, int n, int *vec, int *pepartner) {

  int help,pos;
  int a0,a1,b0,b1,maxp,from,to,p,l,b,i;

  a0=1;maxp=a1=(int) (n/2);
  b0=a1+1;b1=n-1;
    
  for(l=0;l<n;l++) {
    /*   get partners */
    for(p=0;p<(maxp-1);p++) {
      from=vec[a1-p];  to=vec[b0+p];
      if(from==penr) { pepartner[l]=to+1;	    } 
      if(to==penr)   { pepartner[l]=-from-1;    } 
    }
    from=vec[0];	to=vec[a0];
    if(from==penr) { pepartner[l]=to+1;  } 
    if(to==penr)   { pepartner[l]=-from-1; } 


    if(0) {
      printf("norm l=%2d: ",l);
      printf(" %2d -",vec[0]);
      for(i=0;i<n-1;i++) {
	printf(" %2d",vec[i+1]); 
	  if(i==0) printf("  ");
	  if(i==maxp-1) printf("  - ");
      }
      printf(" penr=%2d p=%2d partner[%2d]=%2d\n",penr,p,l,pepartner[l]); 
    }


    /*  rotate */
    help=vec[b1];
    for(b=b1;b>0;b--) {
      vec[b]=vec[b-1];
    }
    vec[a0]=help;
  }

  return(1);
}


/*
             a1
    a0 2  3  4               1  2  3
 0 - 1 |  |  |     -->	 0-7 |  |  |
       7  6  5 	             6  5  4
       b1    b0  
*/

#define v(i) vecn1[(i-l+(n-1))%(n-1)]

int getpartneropt (int penr, int n, int *vec, int *pepartner) {

  int help,pos;
  int a1,a2,an,b1,bn;
  int maxp,from,to,p,l,b,i;
  int *vecn1, pnext=-1;


  vecn1=vec+1;

  a1=0;
  a2=1;
  an=(int) (n/2) -1;
  b1=an+1;
  bn=n-2;

  printf("n=%d: a1=%d, a2=%d,an=%d b1=%d bn=%d\n",n,a1,a2,an,b1,bn);

  if(vec[0]==penr) {
    for(l=0;l<n;l++) {
      pepartner[l]=v(0)+1;          /* p -> x */
    }
  } else {
    p=-1;
    for(pos=a1;(pos<=bn && (p==-1));pos++) {
      if(vecn1[pos]==penr) p=pos;
    }
    
    for(l=0;l<n;l++) {
      if((p>=a2) && (p<=an)) {
	pepartner[l]=v(bn-(p-a2))+1; /* p -> x */
	p=p+1;
      }else if((p>=b1) && (p<=bn)) {
	pepartner[l]=-v(a2+(bn-p))-1; /* x <- p */
	p=p+1;
	if(p>bn) p=a1;
      } else if(p==a1) {
	pepartner[l]=-vec[0]-1;       /* x <- p */
	p=p+1;
      }

      if(0) {
	printf("opt  l=%2d: ",l);
	printf(" %2d -",vec[0]);
	for(i=0;i<n-1;i++) {
	  printf(" %2d",v(i)); 
	  if(i==0) printf("  ");
	  if(i==an) printf("  - ");
	}
	printf(" penr=%2d p=%2d partner[%2d]=%2d\n",penr,p,l,pepartner[l]); 
      }

    }
    
  }
    

  
  return(1);
}
