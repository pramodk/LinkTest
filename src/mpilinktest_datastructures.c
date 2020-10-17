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

_linktest_spec * _linktest_alloc_spec() {
  _linktest_spec *linktest_spec;
  
  linktest_spec = (_linktest_spec *) malloc(sizeof(_linktest_spec));
  if (linktest_spec == NULL) {
    _linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
			 "cannot allocate spec structure of size %lu (linktest_spec), aborting ...\n", (unsigned long) sizeof(_linktest_spec));
    return(NULL);
  }
  
  return(linktest_spec);
}

int _linktest_init_spec( _linktest_spec *linktest_spec ) {

  int       rc = 1;

  if (linktest_spec == NULL) {
    return(_linktest_errorprint(0,_LINKTEST_ERROR_RETURN,"_linktest_init_spec: cannot initalized, data structure is not allocated, aborting ...\n"));
  }
  linktest_spec->rank               = -1;
  linktest_spec->size               = -1;

  linktest_spec->number_of_messages        = number_of_messages_default;
  linktest_spec->number_of_warmup_messages = number_of_warmup_messages_default;
  linktest_spec->length_of_message  = (128*1024);
  linktest_spec->do_alltoall        = 0;
  linktest_spec->sleeptime          = 0;
  linktest_spec->do_serial          = 0;
  linktest_spec->do_mix             = 0;
  linktest_spec->do_sion            = 1;
  linktest_spec->do_sion_par        = 1;
  linktest_spec->collectpnum        = 32;
  linktest_spec->control_niter      = -1;
  linktest_spec->control_wtime      = -1.0;
  linktest_spec->control_current_iter = 0;
  linktest_spec->max_stest          = -1;

  return (rc);
}

int _linktest_print_spec( _linktest_spec *linktest_spec ) {

  int       rc = 1;
  
  if (linktest_spec == NULL) {
    return(_linktest_errorprint(0,_LINKTEST_ERROR_RETURN,"_linktest_print_spec: cannot distribute, data structure is not allocated, aborting ...\n"));
  }
  printf("\n\n");
  printf("----------------------------------------------------------\n");
  printf("linktest: Number of MPI-Task:        %10d\n",linktest_spec->size);
  printf("linktest: Current Rank MPI-Task:     %10d\n",linktest_spec->rank);
  printf("linktest: Message length:            %10d Bytes\n",(int) (linktest_spec->length_of_message*sizeof(char)));
  printf("linktest: Message length:            %10.4f KBytes\n",(float) linktest_spec->length_of_message*sizeof(char)/1024);
  printf("linktest: Number of Iteration:       %10d\n",linktest_spec->number_of_messages);
  printf("linktest: Number of Iter. (Warmup):  %10d\n",linktest_spec->number_of_warmup_messages);
  printf("linktest: alltoall:                  %d\n",linktest_spec->do_alltoall);
  printf("linktest: sleeptime:                 %d us\n",linktest_spec->sleeptime);
  printf("linktest: collectpnum:               %d\n",linktest_spec->collectpnum);
  printf("linktest: mixing pe order:           %d\n",linktest_spec->do_mix);
  printf("linktest: serial test:               %d\n",linktest_spec->do_serial);
  printf("linktest: control wtime:             %f\n",linktest_spec->control_wtime);
  printf("linktest: write protocol (SION):     %d\n",linktest_spec->do_sion);
  printf("linktest:       par I/O:             %d\n",linktest_spec->do_sion_par);
  printf("linktest: max stest:                 %d\n",linktest_spec->max_stest);
  printf("----------------------------------------------------------\n");

  return(1);
}


_linktest_stat * _linktest_alloc_stat() {
  _linktest_stat *linktest_stat;
  
  linktest_stat = (_linktest_stat *) malloc(sizeof(_linktest_stat));
  if (linktest_stat == NULL) {
    _linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
			 "cannot allocate stat structure of size %lu (linktest_stat), aborting ...\n", (unsigned long) sizeof(_linktest_stat));
    return(NULL);
  }
 
   linktest_stat->ptimings=NULL;
   linktest_stat->atimings=NULL;
   linktest_stat->pepartner=NULL;
   linktest_stat->taskvec=NULL;
   char location[MAX_LOCATION_LENGTH];
   linktest_stat->accesspattern=NULL;

   linktest_stat->stimings_new=NULL;
   linktest_stat->stimings_old=NULL;
   linktest_stat->sfrom=NULL;
   linktest_stat->sto=NULL;

  return(linktest_stat);
}

int _linktest_init_stat( _linktest_stat *linktest_stat, _linktest_spec *linktest_spec ) {

  int rc = 1;
  int i;
  int numtasks;

  if (linktest_stat == NULL) {
    return(_linktest_errorprint(0,_LINKTEST_ERROR_RETURN,"_linktest_init_stat: cannot initalized, data structure is not allocated, aborting ...\n"));
  }
  linktest_stat->avg_time_all=0.0; 
  linktest_stat->min_time_all=1e20;
  linktest_stat->max_time_all=-1e20;
  linktest_stat->a2a_avg_time_all=0.0;
  linktest_stat->a2a_min_time_all=1e20;
  linktest_stat->a2a_max_time_all=-1e20;
  
  strncpy(linktest_stat->location,"unknown",MAX_LOCATION_LENGTH-1); 

  /***********************************************************************************/
  /* Allocate and initialize buffers                                                 */
  /***********************************************************************************/
  numtasks=linktest_spec->size;
#ifdef LINKTEST_SETMAXNP
  numtasks=LINKTEST_SETMAXNP;
#endif 

  linktest_stat->ptimings = (double *) malloc(sizeof(double) * numtasks);
  if (linktest_stat->ptimings == NULL) {
    return(_linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
				"cannot allocate ptimings[%lu] ...\n", (unsigned long) sizeof(double) * numtasks));
  }  
  for(i=0;i<numtasks;i++) linktest_stat->ptimings[i]=0.0;

  if(linktest_spec->do_alltoall) {
    linktest_stat->atimings = (double *) malloc(sizeof(double) * numtasks);
    if (linktest_stat->atimings == NULL) {
      return(_linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
				  "cannot allocate atimings[%lu] ...\n", (unsigned long) sizeof(double) * numtasks));
    }  
    for(i=0;i<numtasks;i++) linktest_stat->atimings[i]=0.0;
  }

  linktest_stat->accesspattern = (int32_t *) malloc(sizeof(int32_t) * numtasks);
  if (linktest_stat->accesspattern == NULL) {
    return(_linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
				"cannot allocate accesspattern[%lu] ...\n", (unsigned long) sizeof(int32_t) * numtasks));
  }  
  for(i=0;i<numtasks;i++) linktest_stat->accesspattern[i]=0;

  linktest_stat->taskvec = (int *) malloc(sizeof(int) * numtasks);
  if (linktest_stat->taskvec == NULL) {
    return(_linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
				"cannot allocate taskvec[%lu] ...\n", (unsigned long) sizeof(int) * numtasks));
  }  
  for(i=0;i<numtasks;i++) linktest_stat->taskvec[i]=0;

  linktest_stat->pepartner = (int *) malloc(sizeof(int) * numtasks);
  if (linktest_stat->pepartner == NULL) {
    return(_linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
				"cannot allocate pepartner[%lu] ...\n", (unsigned long) sizeof(int) * numtasks));
  }  
  for(i=0;i<numtasks;i++) linktest_stat->pepartner[i]=0;

  if(linktest_spec->rank==0) {
    linktest_stat->stimings_new = (double *) malloc(sizeof(double) * numtasks);
    if (linktest_stat->stimings_new == NULL) {
      return(_linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
				  "cannot allocate stimings_new[%lu] ...\n", (unsigned long) sizeof(double) * numtasks));
    }  
    for(i=0;i<numtasks;i++) linktest_stat->stimings_new[i]=0.0;
    
    linktest_stat->stimings_old = (double *) malloc(sizeof(double) * numtasks);
    if (linktest_stat->stimings_old == NULL) {
      return(_linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
				  "cannot allocate stimings_old[%lu] ...\n", (unsigned long) sizeof(double) * numtasks));
    }  
    for(i=0;i<numtasks;i++) linktest_stat->stimings_old[i]=0.0;
    linktest_stat->sfrom = (int32_t *) malloc(sizeof(int32_t) * numtasks);
    if (linktest_stat->sfrom == NULL) {
      return(_linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
				  "cannot allocate sfrom[%lu] ...\n", (unsigned long) sizeof(int32_t) * numtasks));
    }  
    for(i=0;i<numtasks;i++) linktest_stat->sfrom[i]=0;

    linktest_stat->sto = (int32_t *) malloc(sizeof(int32_t) * numtasks);
    if (linktest_stat->sto == NULL) {
      return(_linktest_errorprint(-1,_LINKTEST_ERROR_RETURN,
				  "cannot allocate sto[%lu] ...\n", (unsigned long) sizeof(int32_t) * numtasks));
    }  
    for(i=0;i<numtasks;i++) linktest_stat->sto[i]=0;
    
  }


  return (rc);
}

int _linktest_reinit_stat( _linktest_stat *linktest_stat, _linktest_spec *linktest_spec ) {

  int rc = 1;
  int i;

  if (linktest_stat == NULL) {
    return(_linktest_errorprint(0,_LINKTEST_ERROR_RETURN,"_linktest_init_stat: cannot initalized, data structure is not allocated, aborting ...\n"));
  }
  linktest_stat->avg_time_all=0.0; 
  linktest_stat->min_time_all=1e20;
  linktest_stat->max_time_all=-1e20;
  linktest_stat->a2a_avg_time_all=0.0;
  linktest_stat->a2a_min_time_all=1e20;
  linktest_stat->a2a_max_time_all=-1e20;

  /***********************************************************************************/
  /* Allocate and initialize buffers                                                 */
  /***********************************************************************************/
  if (linktest_stat->ptimings != NULL) {
    for(i=0;i<linktest_spec->size;i++) linktest_stat->ptimings[i]=0.0;
  }  

  if (linktest_stat->stimings_new != NULL) {
    for(i=0;i<linktest_spec->size;i++) linktest_stat->stimings_new[i]=0.0;
  }  

  if (linktest_stat->stimings_old != NULL) {
    for(i=0;i<linktest_spec->size;i++) linktest_stat->stimings_old[i]=0.0;
  }  

  if (linktest_stat->atimings != NULL) {
    for(i=0;i<linktest_spec->size;i++) linktest_stat->atimings[i]=0.0;
  }  

  if (linktest_stat->accesspattern != NULL) {
  for(i=0;i<linktest_spec->size;i++) linktest_stat->accesspattern[i]=0;
  }  

  /*
    pre-generated data should not be discarded (from getpartneropt)

  if (linktest_stat->taskvec != NULL) {
    for(i=0;i<linktest_spec->size;i++) linktest_stat->taskvec[i]=0;
  }  

  if (linktest_stat->pepartner != NULL) {
    for(i=0;i<linktest_spec->size;i++) linktest_stat->pepartner[i]=0;
  } 
  */ 
  return (rc);
}


int _linktest_free_stat( _linktest_stat *linktest_stat, _linktest_spec *linktest_spec ) {

  int rc = 1;

  if (linktest_stat == NULL) {
    return(_linktest_errorprint(0,_LINKTEST_ERROR_RETURN,
				"_linktest_free_stat: cannot initalized, data structure is not allocated, aborting ...\n"));
  }

  free(linktest_stat->ptimings); linktest_stat->ptimings=NULL;
  if(linktest_spec->do_alltoall) {
    free(linktest_stat->atimings); linktest_stat->atimings=NULL;
  }
  free(linktest_stat->accesspattern); linktest_stat->accesspattern=NULL;
  free(linktest_stat->taskvec); linktest_stat->taskvec=NULL;
  free(linktest_stat->pepartner); linktest_stat->pepartner=NULL;
  if(linktest_spec->rank==0) {
    free(linktest_stat->stimings_new); linktest_stat->stimings_new=NULL;
    free(linktest_stat->stimings_old); linktest_stat->stimings_old=NULL;
    free(linktest_stat->sfrom); linktest_stat->sfrom=NULL;
    free(linktest_stat->sto); linktest_stat->sto=NULL;
  }


  return (rc);
}



int _linktest_print_stat( _linktest_stat *linktest_stat ) {

  int       rc = 1;
  
  if (linktest_stat == NULL) {
    return(_linktest_errorprint(0,_LINKTEST_ERROR_RETURN,"_linktest_print_stat: cannot distribute, data structure is not allocated, aborting ...\n"));
  }
  printf("\n\n");
  printf("----------------------------------------------------------\n");
  printf("----------------------------------------------------------\n");

  return(1);
}
