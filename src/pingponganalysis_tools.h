/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2011                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/
#ifndef PINGPONGANALYSIS_TOOLS
#define PINGPONGANALYSIS_TOOLS

#include "sion.h"
#include "mpilinktest_mapping.h"

#define MAXNULLCNT 1000

struct analysis_fields_struct {
  sion_int32 *accesspattern;
  double     *accesspatterndp;
  sion_int32 *fromlist;
  sion_int32 *tolist;
  double     *ptimings;
  double     *atimings;
  double     *stimings1;
  double     *stimings2;
  int        *distance;
  int        *mapping;
  long       *ptimings_stepcnt;
  double     *ptimings_step;
  char       *alllocation;
  double      timings_min;
  double      timings_max;
  char       *mask;
  double    *nodeslist_fromavg;
  double    *nodeslist_toavg;
  double    *nodeslist_frommax;
  double    *nodeslist_tomax;
  int       *nodeslist_fromcntmax;
  int       *nodeslist_tocntmax;
};

typedef struct analysis_fields_struct analysis_fields_t;

extern int length_of_message;

void usage(char *name);
double time2mbs( double time);
double mbs2time( double mbs);
int freaddata(void *ptr, long size, FILE *fp);
int allocate_fields ( analysis_fields_t *f, int numtasks, int opt_steps);
int read_maskfile(analysis_fields_t *fields, int numtasks, char *filename);

#endif
