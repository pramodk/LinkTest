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
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "mpi.h"

#include "mpilinktest_util.h"

#define WAKEUP   102
#define COLPRINT 103


double time2mbs( double time, _linktest_spec *linktest_spec) {
  if(time>0) {
    return(linktest_spec->length_of_message*sizeof(char)/time/1024/1024);
  } else {
    return(0.0);
  }
}

int _linktest_errorprint(int rc, int level, const char *format, ...) {
  int       rank=-1;
  va_list ap;
  
  switch (level) {
    case _LINKTEST_ERROR_RETURN: 
      fprintf(stderr,"LINKTEST_ERROR_RETURN on rank %d, rc=%d: ",rank,rc);
      va_start(ap, format);
      vfprintf(stderr, format, ap);
      va_end(ap);
      fprintf(stderr,"\n");
      return (rc);
      break;
    case _LINKTEST_ERROR_WARN: 
      fprintf(stderr,"LINKTEST_ERROR_WARN on rank %d, rc=%d: ",rank,rc);
      va_start(ap, format);
      vfprintf(stderr, format, ap);
      va_end(ap);
      fprintf(stderr,"\n");
      return (rc);
      break;
    case _LINKTEST_ERROR_ABORT: 
      fprintf(stderr,"LINKTEST_ERROR_ABORT on rank %d, rc=%d: \n",rank,rc);
      va_start(ap, format);
      vfprintf(stderr, format, ap);
      va_end(ap);
      fprintf(stderr,"\n");
      exit (rc);
      break;
  default:
      fprintf(stderr,"LINKTEST_ERROR_UNKNOWN on rank %d, rc=%d: ",rank,rc);
      va_start(ap, format);
      vfprintf(stderr, format, ap);
      va_end(ap);
      fprintf(stderr,"\n");
      return (rc);
  }
  return (rc);
}

int collective_print_gather(char *cbuffer)
{
  int       rank, size, p;
  char     *lbuffer;


  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(size*MAXCHARLEN > 64*1024*1024) {
    return(collective_print(cbuffer));
  }

  if(rank==0) {
    lbuffer = (char *) malloc(MAXCHARLEN * size);
    if(!lbuffer) {
      fprintf(stderr,"could allocate buffer of size %d\n",MAXCHARLEN * size);
      MPI_Abort(MPI_COMM_WORLD,1);
    }
  }
  else  lbuffer = NULL;


   
  MPI_Gather(cbuffer, MAXCHARLEN, MPI_CHAR, lbuffer, MAXCHARLEN, MPI_CHAR, 0, MPI_COMM_WORLD);

  if (rank == 0) {

    for (p = 0; p < size; p++) {
      /*       fprintf(stderr,"%06d:%s\n",p,lbuffer+p*MAXCHARLEN);  */
      fprintf(stderr, lbuffer + p * MAXCHARLEN);
    }
  }

  if(rank==0) free(lbuffer);

  return (1);
}

int collective_print(char *cbuffer)
{
  int       rank, size, p;
  int       dummy = 0;
  char      lbuffer[MAXCHARLEN];
  MPI_Status status;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    /*     fprintf(stderr,"%06d:%s",0,cbuffer); */
    fprintf(stderr, cbuffer);

    for (p = 1; p < size; p++) {
      if (p > 0) {
        MPI_Send(&dummy, 1, MPI_INT, p, WAKEUP, MPI_COMM_WORLD);
        MPI_Recv(lbuffer, MAXCHARLEN, MPI_CHAR, p, COLPRINT, MPI_COMM_WORLD, &status);
	/*      fprintf(stderr,"%06d:%s",p,lbuffer); */
        if (strlen(lbuffer) > 0)
          fprintf(stderr, lbuffer);
      }
    }

  }
  else {

    MPI_Recv(&dummy, 1, MPI_INT, 0, WAKEUP, MPI_COMM_WORLD, &status);
    MPI_Send(cbuffer, MAXCHARLEN, MPI_CHAR, 0, COLPRINT, MPI_COMM_WORLD);

  }
  return (1);
}
 

int perm (int step, int n, int *line) {
    int i,b,help;
    int numblocks=n/step;
    for(b=0;b<numblocks;b+=2) {
	for(i=0;i<step;i++) {
	    help=line[b*step+i];
	    line[b*step+i]=line[(b+1)*step+i];
	    line[(b+1)*step+i]=help;
	}
    }
  return(1);
}

/*
     a0      a1
       2  3  4               1  2  3
   0-1 |  |  |     -->	 0-7 |  |  |
       7  6  5 	             6  5  4
       b1    b0  
*/

int getpartnerplan_pe (int penr, int n, int *pepartner, int do_mix, int rank) {
  int *line=NULL;
#ifdef OLDRANDOMGEN      
  int *helpline=NULL;
#endif
  int l,i,b,help,sum;
  int a0,a1,b0,b1,maxp,from,to,p;
  double starttime, timedelay;

  starttime = MPI_Wtime();
  
  line = (int *) malloc(sizeof(int)*n);
  if (line == NULL) {
    fprintf(stderr, "cannot allocate line[%lu] ...\n", (unsigned long) sizeof(int)*n);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
#ifdef OLDRANDOMGEN      
  helpline = (int *) malloc(sizeof(int)*n);
  if (helpline == NULL) {
    fprintf(stderr, "cannot allocate helpline[%lu] ...\n", (unsigned long) sizeof(int)*n);
    MPI_Abort(MPI_COMM_WORLD,1);
  }
#endif

  if(do_mix) {
    if(rank==0) {
      /* generate randomized vector */
#ifdef OLDRANDOMGEN      
      int pos;
      for(i=0;i<n;i++) {helpline[i]=i;}
      srand(1);
      for(i=0;i<n;i++) {
        do {
          pos= (int) ((double) rand() * (double) (n) / (double) RAND_MAX );
/*           printf("%03d: pos=%3d -->%3d\n",i,pos,helpline[pos]);  */
        } while((helpline[pos]<0 || (pos>=n)));
        line[i]=pos;
        helpline[pos]=-1;
/*         printf("line[%03d]=%03d (%d)\n",i,pos,helpline[pos]);  */
      }
#else
      int pos;
      for(i=0;i<n;i++) {line[i]=i;}
      for(i=0;i<n;i++) {
	pos= (int) ((double) rand() * (double) (n) / (double) RAND_MAX );
/* 	printf(" exch i=%3d <-> pos=%3d\n",i,pos);  */
	help=line[i];
        line[i]=line[pos];
	line[pos]=help;
      }
#endif      
    }

    timedelay = MPI_Wtime()-starttime;
    if(rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",rank,"[mixpe]",timedelay);
    starttime = MPI_Wtime();
    
    MPI_Bcast(&line[0],n,MPI_INT,0,MPI_COMM_WORLD);

    timedelay = MPI_Wtime()-starttime;
    if(rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",rank,"[bcast]",timedelay);

    sum=0;
    for(i=0;i<n;i++) {sum+=line[i];}
/*     printf("PE[%06d] do_mix=1 sum=%08d\n",penr,sum);   */
    
  } else {
    for(i=0;i<n;i++) {line[i]=i;}

    sum=0;
    for(i=0;i<n;i++) {sum+=line[i];}
/*     printf("PE[%06d] domix=0 sum=%08d\n",penr,sum);   */
  }
  
#ifdef OUTPUTRANDOMDIST
  if(rank==0) {
    printf("Random distribution: "); 
    for(i=0;i<n;i++) { printf("%d,",line[i]); };
    printf("\n");
  }
#endif

  a0=1;maxp=a1=(int) (n/2);
  b0=a1+1;b1=n-1;
    
  for(l=0;l<n;l++) {
    /*   get partners */
    for(p=0;p<(maxp-1);p++) {
      from=line[a1-p];  to=line[b0+p];
      if(from==penr) { pepartner[l]=to+1;	    } 
      if(to==penr)   { pepartner[l]=-from-1;    } 
    }
    from=line[0];	to=line[a0];
    if(from==penr) { pepartner[l]=to+1;  } 
    if(to==penr)   { pepartner[l]=-from-1; } 

    /*  rotate */
    help=line[b1];
    for(b=b1;b>0;b--) {
      line[b]=line[b-1];
    }
    line[a0]=help;
  }

  timedelay = MPI_Wtime()-starttime;
  if(rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",rank,"[pelist]",timedelay);

  return(1);
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

/*
             a1
    a0 2  3  4               1  2  3
 0 - 1 |  |  |     -->	 0-7 |  |  |
       7  6  5 	             6  5  4
       b1    b0  
*/

#define v(i) vecn1[(i-l+(n-1))%(n-1)]

int getpartneropt (int penr, int n, int *vec, int *pepartner) {

  int pos;
  int a1,a2,an,b1,bn;
  int p,l,i;
  int *vecn1;


  vecn1=vec+1;

  a1=0;
  a2=1;
  an=(int) (n/2) -1;
  b1=an+1;
  bn=n-2;

  /* printf("n=%d: a1=%d, a2=%d,an=%d b1=%d bn=%d\n",n,a1,a2,an,b1,bn);  */

  if(vec[0]==penr) {
    for(l=0;l<n;l++) {
      pepartner[l]=v(0)+1;          /* p -> x */
      /* printf(" penr=%2d p=%2d partner[%2d]=%2d\n",penr,p,l,pepartner[l]);  */
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
