/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2011                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include <sys/times.h>
#include <sys/resource.h>

#define VERSION  3

#define FNAMELEN 255

#define LEVEL_us  1000 
#define LEVEL_mbs 100 
#define proc_A 0
#define proc_B 1
#define PING 100
#define PONG 101
#define WAKEUP 102
#define PTIMINGS 103
#define LOCATION 104
#define MAPPING  105
#define ACCESS   106
#define ATIMINGS 107
#define number_of_messages_default 3
#define length_of_message_max   (1*1024*1024) 
#define length_of_message_a2a   (72*4*1024) 
#define length_of_message_default   (128*1024) 

int number_of_messages = 3;
int length_of_message = (128*1024);

double avg_time_all=0.0, min_time_all=0.0, max_time_all=0.0;
double a2a_avg_time_all=0.0, a2a_min_time_all=0.0, a2a_max_time_all=0.0;
double *ptimings=NULL;
double *stimings=NULL;
double *atimings=NULL;
int        *pepartner=NULL;


struct rusage rusagestr;
int print_memusage(int i);

double time2mbs( double time) {
  if(time>0) {
    return(length_of_message*sizeof(char)/time/1024/1024);
  } else {
    return(0.0);
  }
}


double alltoall (	      int my_rank,
			      int tasks,
			      int iterations
			      ) {
  double timeused=0.0;
  char in[length_of_message_a2a];
  char out[length_of_message_a2a];
  int i,maxlen;   
  double start, finish, time;
  maxlen=1;
  for(i=0;i<tasks*maxlen;i++)  out[i] = my_rank%256;
  start = MPI_Wtime();
  MPI_Alltoall(out,maxlen,MPI_BYTE,in,maxlen,MPI_BYTE,MPI_COMM_WORLD);
  finish = MPI_Wtime();
  
  time = finish - start;
  timeused=(float)(time);
  /* '!!! ckeck results !!! */
  return(timeused);
}

double gather (    	      int my_rank,
			      int size,
			      int datasize,
			      int iterations
			      ) {
  double timeused=0.0;
  char* in=NULL;;
  char* out=NULL;;
  int i;   
  double start, finish, time;

  in = (char *) malloc(sizeof(char)*size*datasize);
  if (in == NULL) {
    fprintf(stderr, "cannot allocate in[%lu] ...\n", (unsigned long) sizeof(char)*size*datasize);
    MPI_Abort(MPI_COMM_WORLD,1);
  }  
  out = (char *) malloc(sizeof(char)*datasize);
  if (out == NULL) {
    fprintf(stderr, "cannot allocate out[%lu] ...\n", (unsigned long) sizeof(char)*datasize);
    MPI_Abort(MPI_COMM_WORLD,1);
  }  

  for(i=0;i<datasize;i++)  out[i] = my_rank%256;

  start = MPI_Wtime();
  MPI_Gather(out,datasize,MPI_CHAR,in,datasize,MPI_CHAR,0,MPI_COMM_WORLD);
  finish = MPI_Wtime();
  
  time = finish - start;
  timeused=(float)(time);

  /* '!!! ckeck results !!! */
  free(out);out=NULL;
  free(in);in=NULL;
  return(timeused);
}


void usage(char *name) {
  fprintf(stderr, "Usage: %s options\n\nwith the following optional options (default values in parathesis):\n\n",name);

  fprintf(stderr, "  [-a 0|1]                  do alltoall mode (0 or 1)     (0) \n");
  fprintf(stderr, "  [-i <iterations>]         number of pingpong iterations (3) \n");
  fprintf(stderr, "  [-s <size>]               message size in Bytes    (131072) \n");
  fprintf(stderr, "  [-k <size>]               message size in KBytes     (128k) \n");
  fprintf(stderr, "  [-T <min>]                minmal runtime in min        (-1) \n");
  fprintf(stderr, "  [-M 0|1]                  randomized processor numbers  (0) \n");
  fprintf(stderr, "  [-S 0|1]                  do searialized test          (no) \n");
  fprintf(stderr, "  [-W 0|1]                  write full protocol           (1) \n");
  exit(1);
}

int main(int argc, char *argv[])
{
  int my_rank,size;
  double starttime,timedelay,globalstarttime;
  int i;
  double     timeused;
  int        dsize;
  int        do_alltoall;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  do_alltoall=   0;

  if(my_rank<=3) {  
    print_memusage(my_rank);
  }

  if(my_rank==0) {
    /* parse command line */
    i=1;
    while( i < argc ) {
      if( argv[i][0] == '-' ) {
        switch( argv[i][1] ) {
        case 'a':
          do_alltoall = atoi(argv[++i]);
          break;
        default:
          usage(argv[0]);
        }
      } else {
        usage(argv[0]);
      }
      i++;
    }
  }

  MPI_Bcast(&do_alltoall,1, MPI_INT, 0, MPI_COMM_WORLD);

  /*                                    */  globalstarttime = starttime = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  /*                                    */  timedelay = MPI_Wtime()-starttime;
  if(my_rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",my_rank,"[first sync]",timedelay);

  if(my_rank==0) { 
    printf("\n\n");
    printf("----------------------------------------------------------\n");
    printf("mpilinktest: Number of MPI-Task:        %10d\n",size);
    printf("mpilinktest: Number of Iteration:       %10d\n",number_of_messages);
    printf("mpilinktest: alltoall:                  %d\n",do_alltoall);
    printf("----------------------------------------------------------\n");
  } 

  
  for(dsize=1;dsize<136000;dsize*=2) {
    timeused=gather(my_rank,size,dsize,1);
    if(my_rank==0) printf("timings[%03d] %-30s datasize=%6d t=%10.6fs\n",my_rank,"[gather]",dsize,timeused);
  }

  /*                                    */  timedelay = MPI_Wtime()-globalstarttime;
  if(my_rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",my_rank,"[all]",timedelay);
  MPI_Finalize();
  return(0);
}




struct rusage rusagestr;

int print_memusage(int i)
{
  if (getrusage(RUSAGE_SELF, &rusagestr) != 0)
    return -1;

  printf("PE: %04d: MEMUSAGE: %12.8f\n",i, rusagestr.ru_maxrss / 1024.0);

  return 1;
} 


