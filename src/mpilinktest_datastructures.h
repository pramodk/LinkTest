/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2011                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/
#ifndef MPILINKTEST_DATASTRUCTURES_H
#define MPILINKTEST_DATASTRUCTURES_H

#define VERSION  1
#define VERSIONSUB  2
#define VERSIONPATCHLEVEL  1

#define FNAMELEN 255
#define MAXCHARLEN 256

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
#define number_of_warmup_messages_default 0
#define LENGTH_OF_MESSAGE_MAX       (1*1024*1024) 
#define length_of_message_a2a       (72*4*1024) 
#define length_of_message_default   (128*1024) 
#define MAX_LOCATION_LENGTH 256
#define LINKTEST_ID "LKTST"  /*!< Linktest identification string (offset: 0x00) */

#if VERSIONSUB>=4
#define MAPPING_DIM 6
#else
#define MAPPING_DIM 4
#endif


#ifdef USESION
  #include "sion.h"
#endif

typedef struct _linktest_spec_struct _linktest_spec;

struct _linktest_spec_struct {
  int rank;
  int size;

  int number_of_messages;
  int number_of_warmup_messages;
  int length_of_message;
  int do_alltoall;
  int sleeptime;
  int do_mix;
  int do_serial;
  int collectpnum;
  int do_sion;
  int do_sion_par;
  int control_current_iter;
  int control_niter;
  double control_wtime;
  int max_stest;
}; 

_linktest_spec * _linktest_alloc_spec();
int              _linktest_init_spec ( _linktest_spec *linktest_spec );
int              _linktest_distribute_spec( _linktest_spec *linktest_spec );
int              _linktest_print_spec( _linktest_spec *linktest_spec );


typedef struct _linktest_stat_struct _linktest_stat;

struct _linktest_stat_struct {
  double avg_time_all; 
  double min_time_all;
  double max_time_all;
  double a2a_avg_time_all;
  double a2a_min_time_all;
  double a2a_max_time_all;
  double *ptimings;
  double *stimings_new;
  double *stimings_old;
  int32_t *sfrom;
  int32_t *sto;
  double *atimings;
  int    *pepartner;
  int    *taskvec;
  char   location[MAX_LOCATION_LENGTH];
  int    mapping[MAPPING_DIM];
  int32_t *accesspattern;
}; 

_linktest_stat * _linktest_alloc_stat();
int              _linktest_init_stat ( _linktest_stat *linktest_stat, _linktest_spec *linktest_spec );
int              _linktest_reinit_stat ( _linktest_stat *linktest_stat, _linktest_spec *linktest_spec );
int              _linktest_distribute_stat( _linktest_stat *linktest_stat );
int              _linktest_print_stat( _linktest_stat *linktest_stat );
int              _linktest_free_stat( _linktest_stat *linktest_stat, _linktest_spec *linktest_spec );

#endif
