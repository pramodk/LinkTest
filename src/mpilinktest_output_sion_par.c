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

#include "mpilinktest_util.h"
#include "mpilinktest_datastructures.h"
#include "mpilinktest_output_sion.h"
#include "mpilinktest_output_sion_par.h"

#include "sion.h"

#define WAKEUP 102
#define BUFFERID 102

int linktest_output_sion_parallel_mpi ( _linktest_spec *linktest_spec, _linktest_stat *linktest_stat ) {

#ifdef USE_SION
  double     starttime, timedelay;
  MPI_Comm   iocomm;
  sion_int64 sion_chunksize,sion_chunksize_sum;
  sion_int32 sion_fsblksize;
  int        sion_numfiles;
  char       filename_bin[256];
  char      *newfname;
  MPI_Comm   comm_local;
  FILE      *fp;
  int        bwrote,sid,grank;
  char      *my_buffer;
  long       my_buffer_size;      

  /***********************************************************************************/
  /* collect local data                                                              */
  /***********************************************************************************/
  
  linktest_output_sion_collect_local_data ( linktest_spec, linktest_stat, &my_buffer, &my_buffer_size);


  /***********************************************************************************/
  /* output to sion file                                                             */
  /***********************************************************************************/

  /* Open SION file */
  MPI_Barrier(MPI_COMM_WORLD);
  /*                                    */  starttime = MPI_Wtime();
  sion_numfiles = 1;
  sion_fsblksize=-1;
  grank=linktest_spec->rank;
  sion_chunksize=my_buffer_size;
  iocomm=MPI_COMM_WORLD;
  newfname=NULL;
  if ( (linktest_spec->control_wtime>0) || (linktest_spec->control_niter>0) ) {
    sprintf(filename_bin, "%s_i%03d.sion", "pingpong_results_bin",linktest_spec->control_current_iter);
  } else {
    sprintf(filename_bin, "%s.sion", "pingpong_results_bin");
  }
  sid=sion_paropen_mpi(filename_bin,"bw,ansi,collective",
		       &sion_numfiles,
		       iocomm,
		       &comm_local,
		       &sion_chunksize,
		       &sion_fsblksize,
		       &grank,
		       &fp,
		       &newfname);

  MPI_Barrier(MPI_COMM_WORLD);
  /*                                    */  timedelay = MPI_Wtime()-starttime;
  if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[sionopen]",timedelay);
    
  /*                                    */  starttime = MPI_Wtime();
    
  /* write buffer */
  bwrote=sion_coll_fwrite_mpi(my_buffer,1,sion_chunksize,sid);
  if (bwrote != sion_chunksize) {
    _linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
			 "problems writing data to sionfile on rank %d, my_buffer[%lu] ...\n", 
			 linktest_spec->rank, (unsigned long) sizeof(char) * sion_chunksize);
  }  
  /*                                    */  timedelay = MPI_Wtime()-starttime;
  MPI_Reduce(&sion_chunksize,&sion_chunksize_sum, 1, SION_MPI_INT64 ,MPI_SUM,0,MPI_COMM_WORLD);

  if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[sionwrite]",timedelay);
  if(linktest_spec->rank==0) printf("OUTPUT: writing data to file %-30s (%10.4f MB)\n",filename_bin,sion_chunksize_sum/1024.0/1024.0);

  /*                                    */  starttime = MPI_Wtime();
  sion_parclose_mpi(sid);
  /*                                    */  timedelay = MPI_Wtime()-starttime;
  if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[sionclose]",timedelay);

#endif
  return(1);
}


int linktest_output_sion_2stage_mpi ( _linktest_spec *linktest_spec, _linktest_stat *linktest_stat ) {

#ifdef USE_SION
  double     starttime, timedelay;
  char       filename_bin[256]="pingpong_results_bin.sion";
  char      *my_buffer;
  long       my_buffer_size;      

  sion_int64 *sion_chunksizes = NULL;
  int        *sion_globalranks = NULL;
  sion_int64 sion_chunksize, sion_chunksize_sum;
  sion_int32 sion_fsblksize;
  int         ntasks = linktest_spec->size;
  int         nfiles = 1;
  int         outsid;
  FILE       *outfp;
  int         t,dummy;
  MPI_Status status;

  /*                                    */  starttime = MPI_Wtime();
  /***********************************************************************************/
  /* collect local data                                                              */
  /***********************************************************************************/
  linktest_output_sion_collect_local_data ( linktest_spec, linktest_stat, &my_buffer, &my_buffer_size);
  /* printf("linktest_output_sion_collect_local_data[%d], got local data of size %ld\n",linktest_spec->rank,(long) my_buffer_size); */

  /* collect chunk_sizes */
  sion_chunksize=my_buffer_size;
  if(linktest_spec->rank==0) {
    sion_chunksizes = (sion_int64 *) malloc(ntasks * sizeof(sion_int64));
    if (sion_chunksizes == NULL) {
      fprintf(stderr,"linktest_output_sion_2stage_mpi: cannot allocate chunksizes of size %lu, aborting ...\n",
	      (unsigned long) ntasks * sizeof(sion_int64));
      return(0);
    }
  }
  MPI_Gather(&sion_chunksize, 1, SION_MPI_INT64,  sion_chunksizes, 1, SION_MPI_INT64, 0, MPI_COMM_WORLD);

  MPI_Barrier(MPI_COMM_WORLD);
  /*                                    */  timedelay = MPI_Wtime()-starttime;
  if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[sioncollect]",timedelay);
    
  /*                                    */  starttime = MPI_Wtime();
  /***********************************************************************************/
  /* serial output to sion file                                                      */
  /***********************************************************************************/
  sion_globalranks = (int *) malloc(ntasks * sizeof(int));
  if (sion_globalranks== NULL) {
    fprintf(stderr,"linktest_output_sion_2stage_mpi: cannot allocate sion_globalranks size %lu, aborting ...\n",
	    (unsigned long) ntasks * sizeof(int));
    return(0);
  }
  for(t=0;t<ntasks;t++) {
    sion_globalranks[t]=t;
  }

  /* open the file on rank 0 and do serial collect and write */
  if(linktest_spec->rank==0) {
    sion_fsblksize=-1;
    nfiles=16;
    outsid = sion_open(filename_bin, "wb", &ntasks, &nfiles, &sion_chunksizes, &sion_fsblksize, &sion_globalranks, &outfp);
    if(outsid<0) {
      fprintf(stderr,"linktest_output_sion_2stage_mpi: sion_open returned %d\n",outsid);
    }
    /* printf("linktest_output_sion_collect_local_data[%d], sion file opened on rank 0\n",linktest_spec->rank); */
    
    for(t=0;t<linktest_spec->size;t++) {
 	if(t>0) {
	  /* printf("linktest_output_sion_collect_local_data[%d], send dummy to rank %d\n",linktest_spec->rank,t); */
	  MPI_Send(&dummy, 1, MPI_INT, t, WAKEUP, MPI_COMM_WORLD);
	  MPI_Recv(my_buffer, my_buffer_size, MPI_CHAR, t, BUFFERID, MPI_COMM_WORLD, &status);
	  /* printf("linktest_output_sion_collect_local_data[%d], got buffer from rank %d\n",linktest_spec->rank,t); */
	  sion_chunksize=sion_chunksizes[t];
	}
	/* printf("linktest_output_sion_collect_local_data[%d], seek pos for rank %d\n",linktest_spec->rank,t); */
	sion_seek_fp(outsid,t,SION_CURRENT_BLK,SION_CURRENT_POS,&outfp);
	/* printf("linktest_output_sion_collect_local_data[%d], write data of rank %d\n",linktest_spec->rank,t); */
	sion_fwrite(my_buffer,1,sion_chunksize,outsid);

    }
    /*                                    */  timedelay = MPI_Wtime()-starttime;
    if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[sioncollwr]",timedelay);

    /*                                    */  starttime = MPI_Wtime();
    sion_close(outsid);
    /*                                    */  timedelay = MPI_Wtime()-starttime;
    if(linktest_spec->rank==0) printf("timings[%03d] %-30s t=%10.6fs\n",linktest_spec->rank,"[sionclose]",timedelay);
    
  } else {
    /* printf("linktest_output_sion_collect_local_data[%d], wait for dummy from rank %d\n",linktest_spec->rank,0); */
    MPI_Recv(&dummy, 1, MPI_INT, 0, WAKEUP, MPI_COMM_WORLD,&status);
    /* printf("linktest_output_sion_collect_local_data[%d], set data to rank %d\n",linktest_spec->rank,0); */
    MPI_Send(my_buffer, sion_chunksize, MPI_CHAR, 0, BUFFERID, MPI_COMM_WORLD);
  }

  MPI_Reduce(&sion_chunksize,&sion_chunksize_sum, 1, SION_MPI_INT64 ,MPI_SUM,0,MPI_COMM_WORLD);

  if(linktest_spec->rank==0) printf("             %-30s %lld bytes\n","[sionwrite]",sion_chunksize_sum);


#endif
  return(1);
}
