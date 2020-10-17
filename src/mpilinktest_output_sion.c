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

#include "mpilinktest_util.h"
#include "mpilinktest_datastructures.h"

#include "sion.h"

int linktest_output_sion_parallel ( _linktest_spec *linktest_spec ) {

#ifdef USE_SIONser
  MPI_Comm   iocomm;
  sion_int64 sion_chunksize;
  sion_int32 sion_fsblksize;
  sion_int32 mappingout[4];
  sion_int32 statdata_int[10];
  int        sion_numfiles;
  int        inioset;
  double     statdata_dp[10];
  char *ptr;
  int        left,bwrote;
  char       filename_bin[256]="pingpong_results_bin.sion";
  char      *newfname;
  MPI_Comm   comm_local;
  MPI_Status status;
  FILE      *fp;
  int        p,sid,grank;

  /***********************************************************************************/
  /* output to sion file                                                             */
  /***********************************************************************************/
  if(linktest_spec->do_sion) {
    sion_numfiles = 1;
    
    if((size%collectpnum)!=0) {
      if(my_rank==0) {
	printf("WARNING: reduce collectpnum to 1 size=%d, collectpnum=%d \n",size,collectpnum);
      }
      collectpnum=1;
    }
    
    sion_chunksize=   MAX_LOCATION_LENGTH;        /* location */
    sion_chunksize+=  4 * sizeof(sion_int32);     /* mapping */
    sion_chunksize+=  size*sizeof(double);        /* timings */
    sion_chunksize+=  size*sizeof(sion_int32);    /* accesspattern */
    if(do_alltoall)   sion_chunksize+=  size*sizeof(double);                   /* alltoall timings, if option selected */
    
    sion_chunksize*=collectpnum;
    
    if(my_rank==0)    sion_chunksize+=10*sizeof(double)+10*sizeof(sion_int32);        /* statistics */
    if(my_rank==0)    sion_chunksize+=size * (2*sizeof(double)+2*sizeof(sion_int32)); /* single tests results */
    
    sion_fsblksize=2*1024*1024;
    
    inioset=((my_rank%collectpnum)==0);
    /*   printf("PE%04d: .. inioset %d / %d\n",my_rank,inioset,collectpnum); */
    MPI_Comm_split(MPI_COMM_WORLD,inioset,my_rank,&iocomm);
    
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    /* open sion file */
    if(inioset) {
      /*                                    */  starttime = MPI_Wtime();
      sid=sion_paropen_mpi(filename_bin,"bw",
			   &sion_numfiles,
			   iocomm,
			   &comm_local,
			   &sion_chunksize,
			   &sion_fsblksize,
			   &grank,
			   &fp,
			   &newfname);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /*                                    */  timedelay = MPI_Wtime()-starttime;
    if(my_rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",my_rank,"[sionopen]",timedelay);
    
    /*                                    */  starttime = MPI_Wtime();
    /*   printf("PE%04d: inioset %d / %d\n",my_rank,inioset,collectpnum); */
    
    if(inioset) {
      for(p=my_rank;p<(my_rank+collectpnum);p++) {
	/*       printf("PE%04d: start loop for %d\n",my_rank,p); */
	
	if(p>my_rank) {
	  /*         printf("PE%04d: send wakeup to %d\n",my_rank,p);fflush(stdout); */
	  MPI_Send(&dummy, 1, MPI_INT, p, WAKEUP, MPI_COMM_WORLD);
	  /*         printf("PE%04d: receive data from %d\n",my_rank,p); */
	  MPI_Recv(location, MAX_LOCATION_LENGTH, MPI_CHAR, p, LOCATION, MPI_COMM_WORLD, &status);
	  MPI_Recv(mapping, 4, MPI_INT, p, MAPPING, MPI_COMM_WORLD, &status);
	  MPI_Recv(ptimings, size, MPI_DOUBLE, p, PTIMINGS, MPI_COMM_WORLD, &status);
	  MPI_Recv(accesspattern, size, MPI_INTEGER4, p, ACCESS, MPI_COMM_WORLD, &status);
	  if(do_alltoall) {
	    MPI_Recv(atimings, size, MPI_DOUBLE, p, ATIMINGS, MPI_COMM_WORLD, &status);
	  }
	  /*         printf("PE%04d: received data from %d\n",my_rank,p); */
	}
	
	/* write location */
	left=MAX_LOCATION_LENGTH;
	ptr=location;
	while(left>0)  {
	  bwrote=fwrite(ptr, 1, left, fp);
	  left-=bwrote;
	  ptr+=bwrote;
	}
	
	/* write mapping */
	mappingout[0]=mapping[0]; mappingout[1]=mapping[1]; mappingout[2]=mapping[2]; mappingout[3]=mapping[3];
	/*       printf("MAPPING[%04d] = %d,%d,%d,%d (%s)\n",p,mapping[0],mapping[1],mapping[2],mapping[3],location); */
	left=4*sizeof(sion_int32);
	ptr=(char *) &mappingout[0];
	while(left>0)  {
	  bwrote=fwrite(ptr, 1, left, fp);
	  left-=bwrote;
	  ptr+=bwrote;
	}
	
	if((my_rank==0) && (p==0))  {
	  /* write statistics on rank 0 */
	  for(i=0;i<10;i++)  statdata_int[i]=0;
	  for(i=0;i<10;i++)  statdata_dp[i]=0.0;
	  statdata_int[0]=number_of_messages;
	  statdata_int[1]=length_of_message;
	  statdata_int[2]=do_alltoall;
	  statdata_int[3]=collectpnum;
	  statdata_int[4]=do_serial;
	  statdata_int[5]=do_mix;
	  statdata_int[9]=VERSION;
	  statdata_dp[0]=min_time_all;
	  statdata_dp[1]=max_time_all;
	  statdata_dp[2]=avg_time_all;    
	  if(do_alltoall) {
	    statdata_dp[3]=a2a_min_time_all;
	    statdata_dp[4]=a2a_max_time_all;
	    statdata_dp[5]=a2a_avg_time_all;
	  }
	  left=10*sizeof(sion_int32);
	  ptr=(char *) &statdata_int[0];
	  while(left>0)  {
	    bwrote=fwrite(ptr, 1, left, fp);
	    left-=bwrote;
	    ptr+=bwrote;
	  }
	  left=10*sizeof(double);
	  ptr=(char *) &statdata_dp[0];
	  while(left>0)  {
	    bwrote=fwrite(ptr, 1, left, fp);
	    left-=bwrote;
	    ptr+=bwrote;
	  }
	}


	/*       for(t=0;t<size;t++) printf("ptimings[%5d]=%e\n",t,ptimings[t]);   */
    
	/* write timings */
	left=size*sizeof(double);
	ptr=(char *) ptimings;
	while(left>0)  {
	  bwrote=fwrite(ptr, 1, left, fp);
	  left-=bwrote;
	  ptr+=bwrote;
	}
	/* write accesspattern */
	left=size*sizeof(sion_int32);
	ptr=(char *) accesspattern;
	while(left>0)  {
	  bwrote=fwrite(ptr, 1, left, fp);
	  left-=bwrote;
	  ptr+=bwrote;
	}


	/* write a2a timings */
	if(do_alltoall) {
	  left=size*sizeof(double);
	  ptr=(char *) atimings;
	  while(left>0)  {
	    bwrote=fwrite(ptr, 1, left, fp);
	    left-=bwrote;
	    ptr+=bwrote;
	  }
	}

	/* write single test results */
	if((my_rank==0) && (p==0))  {

	  /* new timings */
	  left=size*sizeof(double);
	  ptr=(char *) stimings;
	  while(left>0)  {
	    bwrote=fwrite(ptr, 1, left, fp);
	    left-=bwrote;
	    ptr+=bwrote;
	  }
	  
	  /* old timings */
	  for(i=0;i<size;i++)   stimings[i]=sendq->elements[i].time;
	  left=size*sizeof(double);
	  ptr=(char *) stimings;
	  while(left>0)  {
	    bwrote=fwrite(ptr, 1, left, fp);
	    left-=bwrote;
	    ptr+=bwrote;
	  }

	  /* from */
	  for(i=0;i<size;i++)   accesspattern[i]=sendq->elements[i].from;
	  left=size*sizeof(sion_int32);
	  ptr=(char *) accesspattern;
	  while(left>0)  {
	    bwrote=fwrite(ptr, 1, left, fp);
	    left-=bwrote;
	    ptr+=bwrote;
	  }

	  /* to */
	  for(i=0;i<size;i++)   accesspattern[i]=sendq->elements[i].to;
	  left=size*sizeof(sion_int32);
	  ptr=(char *) accesspattern;
	  while(left>0)  {
	    bwrote=fwrite(ptr, 1, left, fp);
	    left-=bwrote;
	    ptr+=bwrote;
	  }
	}
	/*       printf("PE%04d: wrote data for %d\n",my_rank,p); */

      } /* for p */

    } else { /* ! inioset */
      int sendp= ((int) (my_rank/collectpnum)) * collectpnum;
      /*     printf("PE%04d: waiting for wakeup from %d\n",my_rank,sendp); */
      MPI_Recv(&dummy, 1, MPI_INT, sendp, WAKEUP, MPI_COMM_WORLD,&status);
      /*     printf("PE%04d: received wakeup from %d\n",my_rank,sendp); */
      MPI_Send(location, MAX_LOCATION_LENGTH, MPI_CHAR, sendp, LOCATION, MPI_COMM_WORLD);
      MPI_Send(mapping, 4, MPI_INT, sendp, MAPPING, MPI_COMM_WORLD);
      MPI_Send(ptimings, size, MPI_DOUBLE, sendp, PTIMINGS, MPI_COMM_WORLD);
      MPI_Send(accesspattern, size, MPI_INTEGER4, sendp, ACCESS, MPI_COMM_WORLD);
      if(do_alltoall) {
	MPI_Send(atimings, size, MPI_DOUBLE, sendp, ATIMINGS, MPI_COMM_WORLD);
      }
      /*     printf("PE%04d: ready with send to %d\n",my_rank,sendp); */
    } /* inioset */

    MPI_Barrier(MPI_COMM_WORLD);
    /*                                    */  timedelay = MPI_Wtime()-starttime;
    if(my_rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",my_rank,"[sionwrite]",timedelay);

    /* flush sion file */
    /*                                    */  starttime = MPI_Wtime();
    if(inioset) {
      sion_flush(sid);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /*                                    */  timedelay = MPI_Wtime()-starttime;
    if(my_rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",my_rank,"[sionflush]",timedelay);

    /* close sion file */
    /*                                    */  starttime = MPI_Wtime();
  
    if(inioset) {
      sion_parclose_mpi(sid);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /*                                    */  timedelay = MPI_Wtime()-starttime;
    if(my_rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",my_rank,"[sionclose]",timedelay);

  } /* (do_sion) */

#endif
  return(1);
}


int linktest_output_sion_collect_local_data ( _linktest_spec *linktest_spec, _linktest_stat *linktest_stat,
					      char **buffer, long  *buffer_size) {
  int        i,j,rc = 1;
  char      *my_buffer;
  size_t     my_buffer_size;      
  char       *ptr;
  sion_int32 mappingout[MAPPING_DIM];
  sion_int32 statdata_int[10];
  double     statdata_dp[10];
  sion_int32 version_int[3];

  /* calculate buffer size */
  my_buffer_size    = strlen(LINKTEST_ID);                                                         /* id */
  my_buffer_size   += 3 * sizeof(sion_int32);                                                      /* version information */
  my_buffer_size   += MAX_LOCATION_LENGTH;                                                         /* location */
  my_buffer_size   += MAPPING_DIM * sizeof(sion_int32);                                            /* mapping */
  if(linktest_spec->rank==0)       my_buffer_size+=10*sizeof(double)+10*sizeof(sion_int32);        /* stat data */
  my_buffer_size   += linktest_spec->size*sizeof(double);                                          /* timings */
  my_buffer_size   += linktest_spec->size*sizeof(int32_t);                                         /* accesspattern */
  if(linktest_spec->do_alltoall)   {
    my_buffer_size += linktest_spec->size*sizeof(double);                                          /* alltoall timings, if option selected */
  }
  if(linktest_spec->rank==0)    {
    my_buffer_size += linktest_spec->size * (2*sizeof(double)+2*sizeof(sion_int32));               /* single tests results */
  }
       
  /* allocate buffer */
  
  my_buffer = (char *) malloc((size_t) my_buffer_size);
  if (my_buffer == NULL) {
    return(_linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
				"cannot allocate my_buffer[%lu] ...\n", (unsigned long) sizeof(char) * my_buffer_size));
  }  

  /* gather local data into my_buffer */
  ptr=my_buffer;

  /* Id */
  memcpy(ptr,LINKTEST_ID,strlen(LINKTEST_ID));    ptr+=strlen(LINKTEST_ID);

  /* version information */
  version_int[0]=VERSION;
  version_int[1]=VERSIONSUB;
  version_int[2]=VERSIONPATCHLEVEL;
  memcpy(ptr,version_int,3 * sizeof(sion_int32));    ptr+=3*sizeof(sion_int32);

  /* location */
  memcpy(ptr,linktest_stat->location,MAX_LOCATION_LENGTH);    ptr+=MAX_LOCATION_LENGTH;

  /* mapping */
  for(i=0;i<MAPPING_DIM;i++)   mappingout[i]=linktest_stat->mapping[i];
  memcpy(ptr,mappingout,MAPPING_DIM*sizeof(sion_int32)); ptr+=MAPPING_DIM*sizeof(sion_int32);
  /* stat data */
  if(linktest_spec->rank==0) {
    for(i=0;i<10;i++)  statdata_int[i]=0;
    for(i=0;i<10;i++)  statdata_dp[i]=0.0;
    statdata_int[0]=linktest_spec->number_of_messages;
    statdata_int[1]=linktest_spec->length_of_message;
    statdata_int[2]=linktest_spec->do_alltoall;
    statdata_int[3]=1;  /* linktest_spec->collectpnum; */
    statdata_int[4]=linktest_spec->do_serial;
    statdata_int[5]=linktest_spec->do_mix;
    statdata_int[9]=VERSION;
    memcpy(ptr,statdata_int,10*sizeof(sion_int32)); ptr+=10*sizeof(sion_int32);

    statdata_dp[0]=linktest_stat->min_time_all;
    statdata_dp[1]=linktest_stat->max_time_all;
    statdata_dp[2]=linktest_stat->avg_time_all;    
    if(linktest_spec->do_alltoall) {
      statdata_dp[3]=linktest_stat->a2a_min_time_all;
      statdata_dp[4]=linktest_stat->a2a_max_time_all;
      statdata_dp[5]=linktest_stat->a2a_avg_time_all;
    }
    memcpy(ptr,statdata_dp,10*sizeof(double)); ptr+=10*sizeof(double);

  }
  /* timings */
  /*
  if(linktest_spec->rank <8) {
    for(j=0;j<linktest_spec->size;j++) {
      printf("WF: %d   %2d->%10.8f\n",linktest_spec->rank,j,linktest_stat->ptimings[j]);
    }
  }
  */
  memcpy(ptr,linktest_stat->ptimings,linktest_spec->size*sizeof(double)); ptr+=linktest_spec->size*sizeof(double);

  /* accesspattern */
  memcpy(ptr,linktest_stat->accesspattern,linktest_spec->size*sizeof(int32_t)); ptr+=linktest_spec->size*sizeof(int32_t);

  /* alltoall */
  if(linktest_spec->do_alltoall)   {
    memcpy(ptr,linktest_stat->atimings,linktest_spec->size*sizeof(double)); ptr+=linktest_spec->size*sizeof(double);
  }
  /* single test results */
  if(linktest_spec->rank==0) {
    /* new timings */
    memcpy(ptr,linktest_stat->stimings_new,linktest_spec->size*sizeof(double)); ptr+=linktest_spec->size*sizeof(double);
    /* old timings */
    memcpy(ptr,linktest_stat->stimings_old,linktest_spec->size*sizeof(double)); ptr+=linktest_spec->size*sizeof(double);
    /* from timings */
    memcpy(ptr,linktest_stat->sfrom,linktest_spec->size*sizeof(int32_t)); ptr+=linktest_spec->size*sizeof(int32_t);
    /* to timings */
    memcpy(ptr,linktest_stat->sto,linktest_spec->size*sizeof(int32_t)); ptr+=linktest_spec->size*sizeof(int32_t);
  }

  *buffer=my_buffer;
  *buffer_size=my_buffer_size;

  _linktest_free_stat(  linktest_stat, linktest_spec ); 

  if(linktest_spec->rank==0) {
    printf("linktest_output_sion_collect_local_data[%d] alloc+init local buffer of size %d bytes for %d tasks\n",linktest_spec->rank,my_buffer_size,linktest_spec->size);
  }

  return(rc);

}

