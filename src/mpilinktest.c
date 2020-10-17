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

#include <mpi.h>

#include "mpilinktest_datastructures.h"
#include "mpilinktest_datastructures_par.h"

#include "mpilinktest_util.h"
#include "mpilinktest_mapping.h"
#include "mpilinktest_output_sion.h"
#include "mpilinktest_output_sion_par.h"
#include "mpilinktest_util_machine.h"
#include "mpilinktest_kernel.h"
#include "limprioqueue.h"


static int work ( int step, _linktest_spec *linktest_spec, _linktest_stat *linktest_stat, lim_queue *lqueue);

static void usage(char *name) {
  fprintf(stderr, "Usage: %s options\n\nwith the following optional options (default values in parathesis):\n\n",name);

  fprintf(stderr, "  [-a 0|1]                  do alltoall mode (0 or 1)                   (0) \n");
  fprintf(stderr, "  [-i <iterations>]         number of pingpong iterations               (3) \n");
  fprintf(stderr, "  [-I <iterations>]         number of warmup pingpong iterations        (0) \n");
  fprintf(stderr, "  [-s <size>]               message size in Bytes                  (131072) \n");
  fprintf(stderr, "  [-k <size>]               message size in KBytes                   (128k) \n");
  fprintf(stderr, "  [-T <min>]                minimal runtime in min                     (-1) \n");
  /* fprintf(stderr, "  [-B <niter>]              number of cycles of the whole meas.        (-1) \n"); */
  fprintf(stderr, "  [-M 0|1]                  randomized processor numbers                (0) \n");
  fprintf(stderr, "  [-S 0|1]                  do serialized test                         (no) \n");
  fprintf(stderr, "  [-W 0|1]                  write full protocol (SION)                  (1) \n");
  fprintf(stderr, "  [-Y 0|1]                  parallel(1) or serial(0) I/O for protocol   (1) \n");
  fprintf(stderr, "  [-Z <num>]                test top <num> slowest pairs again     (ntasks) \n");
  fprintf(stderr, "  [-V ]                     show version                                    \n");
  exit(1);
}

int main(int argc, char *argv[])
{
  int my_rank, my_size, t, i;
  double starttime,timedelay,globalstarttime;
  double starttimewrk,timedelaywrk,wrktime;
  double starttimestep,timedelaystep,steptime;
  double     time;
  int        dummy, step;
  char      *location;
  int       *mapping;
  int        restart_loop; 
  char      cbuffer[MAXCHARLEN];

  _linktest_spec *linktest_spec;
  _linktest_stat *linktest_stat;

  /* for searching slowest pairs */
  lim_queue    lqueue;
  lim_queue    lqueue_help;
  lim_queue   *sendq, *recvq;
  MPI_Datatype node_dt;
  node         testnode[2];
  MPI_Datatype type[4] = {MPI_DOUBLE,MPI_INT,MPI_INT,MPI_UB};
  MPI_Aint     disp[4];
  int          blocklen[4]={1,1,1,1};
  MPI_Aint     base;
  
  /* MPI Init */
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &my_size);

  if(my_size%2==1) {
    if(my_rank==0) {
      fprintf(stderr, "odd number of tasks (%d), please start mpilinktest with even number of tasks   ... aborting\n", my_size);
    }
    MPI_Abort(MPI_COMM_WORLD,1);
  }
  
  /* Allocate and init internal data structures */
  linktest_spec = _linktest_alloc_spec();
  _linktest_init_spec(linktest_spec);
  linktest_spec->rank=my_rank; 
  linktest_spec->size=my_size; 

  if(linktest_spec->rank==0) {
    /* parse command line */
    i=1;
    while( i < argc ) {
      if( argv[i][0] == '-' ) {
        switch( argv[i][1] ) {
        case 'a':
          linktest_spec->do_alltoall = atoi(argv[++i]);
          break;
        case 'M':
          linktest_spec->do_mix = atoi(argv[++i]);
          break;
	  /*         case 'm': */
	  /*           min_bandwidth = (float) atof(argv[++i]); */
	  /*           break; */
        case 'i':
          linktest_spec->number_of_messages = atoi(argv[++i]);
          break;
        case 'I':
          linktest_spec->number_of_warmup_messages = atoi(argv[++i]);
          break;
        case 'D':
          linktest_spec->sleeptime = atoi(argv[++i]);
          break;
        case 'S':
          linktest_spec->do_serial = atoi(argv[++i]);
          break;
        case 'T':
          linktest_spec->control_wtime= atof(argv[++i]) * 60.0;
          break; 
        case 'N':
          linktest_spec->control_niter= atoi(argv[++i]);
          break; 
        case 'C':
          linktest_spec->collectpnum = atoi(argv[++i]);
          break; 
        case 'W':
          linktest_spec->do_sion = atoi(argv[++i]);
          break;
        case 'Y':
          linktest_spec->do_sion_par = atoi(argv[++i]);
          break;
	case 's':
          linktest_spec->length_of_message = atoi(argv[++i]);
          if(linktest_spec->length_of_message>LENGTH_OF_MESSAGE_MAX) {
	    _linktest_errorprint(-1,_LINKTEST_ERROR_WARN,
				 "Message %d size reduced to %d\n",linktest_spec->length_of_message,LENGTH_OF_MESSAGE_MAX);
            linktest_spec->length_of_message=LENGTH_OF_MESSAGE_MAX;
          }
          break;
        case 'k':
          linktest_spec->length_of_message  = atoi(argv[++i]);
          linktest_spec->length_of_message *= 1024;
          if(linktest_spec->length_of_message>LENGTH_OF_MESSAGE_MAX) {
	    _linktest_errorprint(-1,_LINKTEST_ERROR_WARN,
				 "Message %d size reduced to %d\n",linktest_spec->length_of_message,LENGTH_OF_MESSAGE_MAX);
            linktest_spec->length_of_message=LENGTH_OF_MESSAGE_MAX;
          }
          break;
        case 'V':
	  fprintf(stderr,"FZJ Linktest version: %d.%dp%d\n",VERSION,VERSIONSUB,VERSIONPATCHLEVEL);
	  exit(1);
        case 'Z':
          linktest_spec->max_stest = atoi(argv[++i]);
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

  if( (linktest_spec->max_stest<0) || (linktest_spec->max_stest>linktest_spec->size) ) {
    linktest_spec->max_stest=linktest_spec->size;
  }


  _linktest_distribute_spec_mpi(linktest_spec);
  
  /* Allocate and init internal data structures */
  linktest_stat = _linktest_alloc_stat();
  _linktest_init_stat(linktest_stat,linktest_spec);

  /*                                    */  globalstarttime = starttime = MPI_Wtime();
  MPI_Barrier(MPI_COMM_WORLD);
  /*                                    */  timedelay = MPI_Wtime()-starttime;
  if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[first sync]",timedelay);

  if(linktest_spec->rank==0) { 
    _linktest_print_spec(linktest_spec);

  } 

  /***********************************************************************************/
  /* get Mapping                                                                     */
  /***********************************************************************************/
  /*                                    */  starttime = MPI_Wtime();
  get_mapping(&location,&mapping);

#if MAPPING_DIM==6
  sprintf(cbuffer, "task[%06d] on %s (%d,%d,%d,%d,%d,%d) mem=%10.4f MB \n", linktest_spec->rank,location,
	  mapping[0],mapping[1],mapping[2],mapping[3],mapping[4],mapping[5],(double) (get_memusage()/1024.0/1024.0));
#else
  sprintf(cbuffer, "task[%06d] on %s (%d,%d,%d,%d) mem=%10.4f MB \n", linktest_spec->rank,location,
	  mapping[0],mapping[1],mapping[2],mapping[3],get_memusage()/1024.0/1024.0);
#endif
  collective_print_gather(cbuffer);

  strncpy(linktest_stat->location,location,MAX_LOCATION_LENGTH-1); 
  for(i=0;i<MAPPING_DIM;i++)   linktest_stat->mapping[i]=mapping[i];

#ifdef LINKTEST_BGQ
  if(linktest_spec->rank==0) { 
    print_memusage(linktest_spec->rank);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(linktest_spec->rank==1) { 
    print_memusage(linktest_spec->rank);
  }
#endif  

  /*                                    */  timedelay = MPI_Wtime()-starttime;
  MPI_Barrier(MPI_COMM_WORLD);
  if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[mapping]",timedelay);


  /***********************************************************************************/
  /* PRETEST                                                                         */
  /***********************************************************************************/
#ifdef PRETEST
  /*                                    */  starttime = MPI_Wtime();
  if(linktest_spec->rank==0) {
    printf("\n\nPRE-test if all links are up:\n");
    printf("-----------------------------\n");
  }  

  for(t=1;t<size;t++) {
    if(linktest_spec->rank==0) {
      printf(" Task %3d to Task %3d:",0,t);
      fflush(stdout);
    }
    time=pingpong(0,t,linktest_spec->rank,3);
    if(linktest_spec->rank==0) {
      printf(" %10.3fus (%7.2f MB/s)\n",time*1000*1000,time2mbs(time));
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  /*                                    */  timedelay = MPI_Wtime()-starttime;
  if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[pretest]",timedelay);
#endif

  if(linktest_spec->rank==0) {
    printf("\n\nStarting Test of all connections:\n");
    printf("--------------------------------\n");
  }  

  _linktest_reinit_stat(linktest_stat,linktest_spec);

  /***********************************************************************************/
  /* Init partnerplan                                                                */
  /***********************************************************************************/

  /*                                    */  starttime = MPI_Wtime();
  if(linktest_spec->do_mix) {
    getrandomize_vector_swap(linktest_spec->size,linktest_stat->taskvec);
  } else {
    for(i=0;i<linktest_spec->size;i++) {linktest_stat->taskvec[i]=i;}
  }
  /*                                    */  timedelay = MPI_Wtime()-starttime;
  if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[randvec]",timedelay);

  /*                                    */  starttime = MPI_Wtime();
  getpartneropt(linktest_spec->rank,linktest_spec->size,linktest_stat->taskvec,linktest_stat->pepartner);
  /*                                    */  timedelay = MPI_Wtime()-starttime;
  if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[getpart]",timedelay);

  /* getpartnerplan_pe(linktest_spec->rank,size,pepartner,do_mix,linktest_spec->rank);  */
  /* pepartner contains for each step number of communication (from -> to) partner:
     pepartner[step]<0   rank is 'to'   part, 'from'-rank is  -pepartner[step]-1 
     pepartner[step]>0   rank is 'from' part, 'to'-rank   is   pepartner[step]-1 
  */

  /* check partner list */
  MPI_Barrier(MPI_COMM_WORLD);
  {
    long psum=0,pasum=0;
    for(t=0;t<linktest_spec->size;t++) {
      psum+=linktest_stat->pepartner[t];
      if(linktest_stat->pepartner[t]<0)  pasum+= -linktest_stat->pepartner[t];
      else                               pasum+=  linktest_stat->pepartner[t];
    }
    if(linktest_spec->rank%1000==0) {
      printf("PE%05d: psum=%d pasum=%ld do_mix=%ld\n",linktest_spec->rank,psum,pasum,linktest_spec->do_mix);
    }
  }
  if(linktest_spec->rank==0)  printf(" ... after getpartnerplan_pe\n"); 


  /* print partner list */
  if(0) {
    for(t=0;t<linktest_spec->size;t++) { 
      MPI_Barrier(MPI_COMM_WORLD);
      if(t==linktest_spec->rank) {
	printf("partners of %2d: ",linktest_spec->rank);
	for(i=0;i<linktest_spec->size;i++) {
	  if(linktest_stat->pepartner[i]<0)  printf("(%3d->   ) ",-linktest_stat->pepartner[i]-1);
	  else                               printf("[   ->%3d] ", linktest_stat->pepartner[i]-1);
	}
	printf("\n");
      }
    }

  }

  /* preparing MPI datatype for lim_queue */
  {
    MPI_Address(testnode,disp);
    MPI_Address(&testnode[0].from,disp+1);
    MPI_Address(&testnode[0].to,  disp+2);
    MPI_Address(testnode+1,       disp+3);
    base=disp[0];
    for(i=0;i<4;i++) disp[i]-=base;
    if(linktest_spec->rank==0) {
      printf("PE[%d] disp %d %d %d %d base=%10ld sizeof(node)=%d sizeof(MPI_Aint)=%d\n",linktest_spec->rank,
	     (int) disp[0],(int) disp[1],(int) disp[2],(int) disp[3],(long) base,(int) sizeof(node),(int) sizeof(MPI_Aint));
    }
    fflush(stdout);
    MPI_Type_struct(4,blocklen,disp,type,&node_dt);
    MPI_Type_commit(&node_dt);
    MPI_Barrier(MPI_COMM_WORLD);
  }  


  /***********************************************************************************/
  /* Work                                                                            */
  /***********************************************************************************/
  /*                                    */  starttime = MPI_Wtime();

  limqueue_initalize(&lqueue,      linktest_spec->size);
  limqueue_initalize(&lqueue_help, linktest_spec->size);

  restart_loop=1;
  wrktime=steptime=0.0;

  while(restart_loop) {

    /* reinit data structures */
    _linktest_reinit_stat(linktest_stat,linktest_spec);
    limqueue_reinitalize(&lqueue);
    limqueue_reinitalize(&lqueue_help);
    linktest_spec->control_current_iter++;
  
    /*                                    */  starttimestep = MPI_Wtime();
    for(step=1;step<linktest_spec->size;step++) {
      /*                                    */  starttimewrk = MPI_Wtime();
      work(step, linktest_spec, linktest_stat, &lqueue_help);
      /*                                    */  timedelaywrk = MPI_Wtime()-starttimewrk;
      wrktime+=timedelaywrk;

      if((linktest_spec->rank==0) && (step%linktest_spec->collectpnum==0)) {
	
	printf("analyse summary: min. %11.6fus (%7.2f MB/s) max. %11.6fus (%7.2f MB/s) avg. %11.6fus (%7.2f MB/s)\n",
	       linktest_stat->min_time_all*1000*1000,time2mbs(linktest_stat->min_time_all,linktest_spec),
	       linktest_stat->max_time_all*1000*1000,time2mbs(linktest_stat->max_time_all,linktest_spec),
	       linktest_stat->avg_time_all/step*1000*1000,time2mbs(linktest_stat->avg_time_all/step,linktest_spec));
	printf("timing summary:  %d steps used %10.6fs (%10.6fs/step),        to do: %d steps,                   estimate time: %10.6fs\n",
	       step,wrktime,wrktime/step,
	       linktest_spec->size-step,wrktime/step*(linktest_spec->size-step)
	       );
	printf("\n");
	
	
      }
    }
    /*                                    */  timedelaystep = MPI_Wtime()-starttimestep;
    steptime+=timedelaystep;
    MPI_Barrier(MPI_COMM_WORLD);
    /*                                    */  timedelay = MPI_Wtime()-starttime;
    if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[work]",timedelay);
    
    /***********************************************************************************/
    /* Average results                                                                 */
    /***********************************************************************************/
    if(linktest_spec->rank==0) {
      linktest_stat->avg_time_all/=(linktest_spec->size-1);
      linktest_stat->a2a_avg_time_all/=(linktest_spec->size-1);
      printf("\n\n");
      printf("RESULT: Min Time:            %20.8fus (%7.2f MB/s) [iter %2d ]\n",linktest_stat->min_time_all*1000*1000,
	     time2mbs(linktest_stat->min_time_all,linktest_spec),linktest_spec->control_current_iter);
      printf("RESULT: Max Time:            %20.8fus (%7.2f MB/s) [iter %2d ]\n",linktest_stat->max_time_all*1000*1000,
	     time2mbs(linktest_stat->max_time_all,linktest_spec),linktest_spec->control_current_iter);
      printf("RESULT: Avg Time:            %20.8fus (%7.2f MB/s) [iter %2d ]\n",linktest_stat->avg_time_all*1000*1000,
	     time2mbs(linktest_stat->avg_time_all,linktest_spec),linktest_spec->control_current_iter);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    /***********************************************************************************/
    /* searching slowest pairs                                                         */
    /***********************************************************************************/
    /*                                    */  starttime = MPI_Wtime();
    limqueue_sort(&lqueue_help);
    limqueue_merge_in(&lqueue,&lqueue_help);
    limqueue_reinitalize(&lqueue_help);



    /* merging timings for seach of slowest pairs  */
    {
      int stride1=1;
      int stride2=2;
      int sender,receiver;
      MPI_Status status;

      sendq=&lqueue;
      recvq=&lqueue_help;

      while(stride1<linktest_spec->size) {
	if(linktest_spec->rank==0) printf(" %03d: starting collect for stride  %6d\n",linktest_spec->rank,stride1);
	if( (linktest_spec->rank%stride2) == 0) { /* receive data */
	  sender=linktest_spec->rank+stride1;
	  if(sender<linktest_spec->size) {
	    /* 	  printf(" %03d: RECV from %03d %d el\n",linktest_spec->rank,sender,recvq->max_elems);fflush(stdout); */
	    MPI_Recv(recvq->elements,recvq->max_elems,node_dt,sender,1001,MPI_COMM_WORLD,&status);
	    MPI_Get_count(&status,node_dt,&recvq->qsize);
	    /* 	  printf(" %03d: RECV from %03d ready %d el\n",linktest_spec->rank,sender,recvq->qsize);fflush(stdout); */
	    limqueue_merge_in(sendq,recvq);
	  }
	}
	if( (linktest_spec->rank%stride2) == stride1) { /* send data */
	  receiver=linktest_spec->rank-stride1;
	  /* 	printf(" %03d: SEND to   %03d %d el\n",linktest_spec->rank,receiver,sendq->qsize);fflush(stdout);  */
	  MPI_Send(sendq->elements,sendq->qsize,node_dt,receiver,1001,MPI_COMM_WORLD);
	}
	stride1*=2;
	stride2*=2;
      }
    }


    MPI_Bcast(&sendq->qsize,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Bcast(sendq->elements,sendq->qsize,node_dt,0,MPI_COMM_WORLD);

    /*                                    */  timedelay = MPI_Wtime()-starttime;
    if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[search_slow]",timedelay);

    /* test again the slowest pairs */
    /*                                    */  starttime = MPI_Wtime();
    {
      int i;
      MPI_Status status;

      for(i=0;i<linktest_spec->max_stest;i++) {
	if(linktest_spec->rank==0) {
	  printf(" % 6d: PINGPONG % 6d  <->  % 6d: ",i,sendq->elements[i].from,sendq->elements[i].to);fflush(stdout);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(sendq->elements[i].from == linktest_spec->rank) {
	  /* 	printf(" %03d: PINGPONG %03d   ->  %03d\n",linktest_spec->rank,linktest_spec->rank,sendq->elements[i].to);fflush(stdout); */
	  time =pingpong(linktest_spec->rank,                 sendq->elements[i].to, linktest_spec); 
	} 
	if(sendq->elements[i].to == linktest_spec->rank) {
	  /* 	printf(" %03d: PINGPONG %03d  <-   %03d\n",linktest_spec->rank,sendq->elements[i].from,linktest_spec->rank);fflush(stdout); */
	  dummy =pingpong(sendq->elements[i].from, linktest_spec->rank, linktest_spec); 
	}
	MPI_Barrier(MPI_COMM_WORLD);
      
	if(sendq->elements[i].from != 0) {
	  /* send result to task zero */
	  if(sendq->elements[i].from == linktest_spec->rank) {
	    MPI_Send(&time,1,MPI_DOUBLE,0,1002,MPI_COMM_WORLD);
	  }
	  if(linktest_spec->rank==0) {
	    MPI_Recv(&time,1,MPI_DOUBLE,sendq->elements[i].from,1002,MPI_COMM_WORLD,&status);
	  }
	}
	if(linktest_spec->rank==0) {
	  linktest_stat->stimings_new[i]=time;
	  linktest_stat->stimings_old[i]=sendq->elements[i].time;
	  linktest_stat->sfrom[i]=sendq->elements[i].from;
	  linktest_stat->sto[i]=sendq->elements[i].to;

	  printf(" 1st: %10.4fus (%7.2f MB/s)    2nd: %10.4fus (%7.2f MB/s)\n",
		 sendq->elements[i].time*1000*1000,time2mbs(sendq->elements[i].time,linktest_spec),
		 linktest_stat->stimings_new[i]*1000*1000,time2mbs(linktest_stat->stimings_new[i],linktest_spec));
	}
      }
      if(linktest_spec->rank==0) {
	double min_time_stest=1e20;
	double max_time_stest=-1e20;
	double avg_time_stest=0.0;
	for(i=0;i<linktest_spec->max_stest;i++) {
	  if(linktest_stat->stimings_new[i]<min_time_stest) {
	    min_time_stest=linktest_stat->stimings_new[i];
	  }
	  if(linktest_stat->stimings_new[i]>max_time_stest) {
	    max_time_stest=linktest_stat->stimings_new[i];
	  }
	  avg_time_stest+=linktest_stat->stimings_new[i];
	}
	avg_time_stest/=linktest_spec->size;

	printf("RESULT: Min Time Stest:      %20.8fus (%7.2f MB/s) [iter %2d ]\n",
	       min_time_stest*1000*1000,time2mbs(min_time_stest,linktest_spec),linktest_spec->control_current_iter);
	printf("RESULT: Max Time Stest:      %20.8fus (%7.2f MB/s) [iter %2d ]\n",
	       max_time_stest*1000*1000,time2mbs(max_time_stest,linktest_spec),linktest_spec->control_current_iter);
	printf("RESULT: Avg Time Stest:      %20.8fus (%7.2f MB/s) [iter %2d ]\n",
	       avg_time_stest*1000*1000,time2mbs(avg_time_stest,linktest_spec),linktest_spec->control_current_iter);
      }
    }
    /*                                    */  timedelay = MPI_Wtime()-starttime;
    if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[test_slow]",timedelay);

    /* write data to SION file */
    if(linktest_spec->do_sion) {

#ifdef LINKTEST_BGQ
  if(linktest_spec->rank==0) { 
    print_memusage(linktest_spec->rank);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(linktest_spec->rank==1) { 
    print_memusage(linktest_spec->rank);
  }
#endif  

      if(linktest_spec->do_sion_par) {
	linktest_output_sion_parallel_mpi(linktest_spec, linktest_stat);
      } else {
	linktest_output_sion_2stage_mpi(linktest_spec, linktest_stat);
      }

#ifdef LINKTEST_BGQ
  if(linktest_spec->rank==0) { 
    print_memusage(linktest_spec->rank);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if(linktest_spec->rank==1) { 
    print_memusage(linktest_spec->rank);
  }
#endif  

    }


    /* check if test should run again */
    if(linktest_spec->rank==0) {
      if(linktest_spec->control_wtime>0) {
	if( (steptime + timedelaystep) < linktest_spec->control_wtime) restart_loop=1;
	else                                                           restart_loop=0;
      } else if(linktest_spec->control_niter>0) {
	if( linktest_spec->control_current_iter < linktest_spec->control_niter)                  restart_loop=1;
	else                                                           restart_loop=0;
      } else {
	restart_loop=0;
      }
      if(restart_loop==1) {
	fprintf(stdout,"LINKTEST[%2d]: end test not reached: control_wtime=%10.2fs, elapsed time=%10.4fs, estimate step time=%10.4fs niter=%d<=%d\n",
		linktest_spec->control_current_iter,linktest_spec->control_wtime,steptime,timedelaystep, linktest_spec->control_current_iter, linktest_spec->control_niter);fflush(stdout);
	fprintf(stderr,"LINKTEST[%2d]: end test not reached: control_wtime=%10.2fs, elapsed time=%10.4fs, estimate step time=%10.4fs niter=%d<=%d\n",
		linktest_spec->control_current_iter,linktest_spec->control_wtime,steptime,timedelaystep, linktest_spec->control_current_iter, linktest_spec->control_niter);fflush(stderr);

      }
    }
    MPI_Bcast(&restart_loop,1,MPI_INT,0,MPI_COMM_WORLD);
    
  }
  


  /*                                    */  timedelay = MPI_Wtime()-globalstarttime;
  if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[all]",timedelay);



  MPI_Finalize();
  return(0);
}


int work ( int step, _linktest_spec *linktest_spec, _linktest_stat *linktest_stat, lim_queue *lqueue) {

  double timep_max,timep_min,timep,timea_max,timea_min, timea;
  double time, dummy, timep_avg, timea_avg;
  int    partner, partner_rank,l;
  
  
  timea_min=0.0;
  timea_max=0.0;

  if(linktest_spec->rank==0) {
      if(linktest_spec->do_serial)  printf(" Serial   PingPong for step %5d: \n",step);
      else                          printf(" Parallel PingPong for step %5d: ",step);
  }
  if(linktest_spec->do_alltoall) {
    MPI_Barrier(MPI_COMM_WORLD);
    timea=alltoall(linktest_spec);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&timea,&timea_min,1,MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
    MPI_Reduce(&timea,&timea_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
    MPI_Reduce(&timea,&timea_avg,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    timea_avg/=linktest_spec->size;
    linktest_stat->atimings[step]+=timea;
  } else {
    MPI_Barrier(MPI_COMM_WORLD);
  }

  partner=linktest_stat->pepartner[step];
  if(partner<0) partner_rank=-partner-1;
  else          partner_rank= partner-1;
  
  if(!linktest_spec->do_serial) {
    if(partner<0) {
      MPI_Barrier(MPI_COMM_WORLD);
      dummy=pingpong(partner_rank,linktest_spec->rank     ,linktest_spec); 
      MPI_Barrier(MPI_COMM_WORLD);
      time =pingpong(linktest_spec->rank     ,partner_rank,linktest_spec); 
      MPI_Barrier(MPI_COMM_WORLD);
    } else {
      MPI_Barrier(MPI_COMM_WORLD);
      time= pingpong(linktest_spec->rank     ,partner_rank,linktest_spec); 
      MPI_Barrier(MPI_COMM_WORLD);
      dummy=pingpong(partner_rank,linktest_spec->rank     ,linktest_spec); 
      MPI_Barrier(MPI_COMM_WORLD);
    }
  } else {
    
    for(l=0;l<linktest_spec->size;l++) {
     
      MPI_Barrier(MPI_COMM_WORLD);
      if (l==linktest_spec->rank) {
	  time= pingpong(linktest_spec->rank     ,partner_rank,linktest_spec); 
	  printf(" %6d->%6d:%12.5fus (%7.2f MB/s) (l=%d)\n",linktest_spec->rank, partner_rank, 	   time*1000*1000,time2mbs(time,linktest_spec),l);
      }
      if (l==partner_rank) {
	  dummy=pingpong(partner_rank,linktest_spec->rank     ,linktest_spec); 
      }
      MPI_Barrier(MPI_COMM_WORLD);
   }
  }
  linktest_stat->ptimings[partner_rank]+=time;

  /* printf("WF[%4d] set ptimings[%d]=%e (%d)\n",linktest_spec->rank,partner_rank,ptimings[partner_rank],partner);  */
  linktest_stat->accesspattern[partner_rank]=step;
  limqueue_unsort_insert(lqueue, time, linktest_spec->rank, partner_rank);


  /* get results of this step (min,max) */
  timep=time; 
  MPI_Reduce(&timep,&timep_min, 1,  MPI_DOUBLE,MPI_MIN,0,MPI_COMM_WORLD);
  MPI_Reduce(&timep,&timep_max, 1,  MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  MPI_Reduce(&timep,&timep_avg, 1,  MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
  timep_avg/=linktest_spec->size;

  
  if(linktest_spec->rank==0) {

    if(timep_min<linktest_stat->min_time_all) {
      linktest_stat->min_time_all=timep_min;
    }
    if(timep_max>linktest_stat->max_time_all) {
      linktest_stat->max_time_all=timep_max;
    }
    linktest_stat->avg_time_all+=timep_avg;

    if(linktest_spec->do_alltoall) {
	if(timea_min<linktest_stat->a2a_min_time_all) {
	    linktest_stat->a2a_min_time_all=timea_min;
	}
	if(timea_max>linktest_stat->a2a_max_time_all) {
	    linktest_stat->a2a_max_time_all=timea_max;
	}
	linktest_stat->a2a_avg_time_all+=timea_avg;
    }
    printf("  avg=%10.3fus (%7.2f MB/s) min=%10.3fus (%7.2f MB/s) max=%10.3fus (%7.2f MB/s)",
	   timep_avg*1000*1000,time2mbs(timep_avg,linktest_spec),
	   timep_min*1000*1000,time2mbs(timep_min,linktest_spec),
	   timep_max*1000*1000,time2mbs(timep_max,linktest_spec)	   );
    if(linktest_spec->do_alltoall) {
      printf(" [ alltoall min=%10.3fus max=%10.3fus ]\n",
	     timea_min*1000*1000,
	     timea_max*1000*1000                   );
    } else {
      printf("\n");
    }
  }
  
  return(1);
}






