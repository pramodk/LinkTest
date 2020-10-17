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

#include "pingponganalysis_tools.h"
#include "mpilinktest_datastructures.h"


int length_of_message = (128*1024);

void usage(char *name) {
  fprintf(stderr, "Usage: %s options <insionfn> \n\nwith the following optional options (default values in parathesis):\n\n",name);

  fprintf(stderr, "  [-a]                   generate accesspattern file (PPM) \n");
  fprintf(stderr, "  [-A]                   generate accesspattern file (ASCII) \n");
  fprintf(stderr, "  [-b]                   generate bandwidthpattern file (PPM) \n");
  fprintf(stderr, "  [-B]                   generate badlink list (ASCII) \n");
  fprintf(stderr, "  [-d]                   generate distancepattern file (dat) \n");
  fprintf(stderr, "  [-D]                   ... Unknown Option ... \n");
  fprintf(stderr, "  [-f]                   ... Unknown Option ... \n");
  fprintf(stderr, "  [-g]                   generate gnuplot 2d input file \n");
  fprintf(stderr, "  [-l] <minbw>           min. bandwidth, conection below will be reported  (def. 1 MB/s) \n");
  fprintf(stderr, "  [-L] <maxbw>           max. bandwidth, conection above will be reported  (def. 10000 MB/s) \n");
  fprintf(stderr, "  [-m]                   ... Unknown Option ... \n");
  fprintf(stderr, "  [-O]                   ... Unknown Option ... \n");
  fprintf(stderr, "  [-p]                   generate postscript report\n");
  fprintf(stderr, "  [-q]                   ... Unknown Option ... \n");
  fprintf(stderr, "  [-Q]                   ... Unknown Option ... \n");
  fprintf(stderr, "  [-s]                   Specify the number of steps. \n");
  fprintf(stderr, "  [-t]                   min. time \n");
  fprintf(stderr, "  [-T]                   max. time \n");
  fprintf(stderr, "  [-v]                   verbose mode \n");
  fprintf(stderr, "  [-V]                   print version \n");
  exit(1);
}


double time2mbs( double time) {
  if(time>0) {
    return((length_of_message * sizeof(char))/1024.0/1024.0/time );
  } else {
    return(0.0);
  }
}

double mbs2time( double mbs) {
  if(mbs>0) {
    return( (length_of_message*sizeof(char))/(mbs*1024*1024));
  } else {
    return(0.0);
  }
}


int freaddata(void *ptr, long size, FILE *fp) {
  int nullcount;
  size_t left,bread,bsumread,bwrote,bsumwrote;
  char *p=ptr;


  nullcount=0;
  left=size;
  bsumread=0;
  p=ptr;
  while(left>0) {
      bread=fread(p+bsumread,1,left,fp);
      if(bread==0) nullcount++; else nullcount=0;
      left-=bread;
      bsumread+=bread;
      p+=bread;
      if (nullcount>MAXNULLCNT) { fprintf(stderr, "timeout on read data , aborting ...\n"); exit(0);}
  }
  return(bsumread);
}

int allocate_fields ( analysis_fields_t *f, int numtasks, int opt_steps) {
  int t;
  

  f->accesspattern = (sion_int32 *)malloc(numtasks*sizeof(sion_int32));
  if (f->accesspattern==NULL) {
    fprintf(stderr, "cannot allocate accesspattern of size %lld , aborting ...\n", numtasks*sizeof(sion_int32));
    return(1);
  }
  f->accesspatterndp = (double *)malloc(numtasks*sizeof(double));
  if (f->accesspatterndp==NULL) {
    fprintf(stderr, "cannot allocate accesspatterndp of size %lld , aborting ...\n", numtasks*sizeof(double));
    return(1);
  }
  f->fromlist = (sion_int32 *)malloc(numtasks*sizeof(sion_int32));
  if (f->fromlist==NULL) {
    fprintf(stderr, "cannot allocate fromlist of size %lld , aborting ...\n", numtasks*sizeof(sion_int32));
    return(0);
  }
  f->tolist = (sion_int32 *)malloc(numtasks*sizeof(sion_int32));
  if (f->tolist==NULL) {
    fprintf(stderr, "cannot allocate tolist of size %lld , aborting ...\n", numtasks*sizeof(sion_int32));
    return(0);
  }
  f->ptimings = (double *)malloc(numtasks*sizeof(double));
  if (f->ptimings==NULL) {
    fprintf(stderr, "cannot allocate ptimings of size %lld , aborting ...\n", numtasks*sizeof(double));
    return(0);
  } 
  f->atimings = (double *)malloc(numtasks*sizeof(double));
  if (f->atimings==NULL) {
    fprintf(stderr, "cannot allocate atimings of size %lld , aborting ...\n", numtasks*sizeof(double));
    return(0);
  }
  f->stimings1 = (double *)malloc(numtasks*sizeof(double));
  if (f->stimings1==NULL) {
    fprintf(stderr, "cannot allocate stimings1 of size %lld , aborting ...\n", numtasks*sizeof(double));
    return(0);
  }
  f->stimings2 = (double *)malloc(numtasks*sizeof(double));
  if (f->stimings2==NULL) {
    fprintf(stderr, "cannot allocate stimings2 of size %lld , aborting ...\n", numtasks*sizeof(double));
    return(0);
  }
  f->distance = (int *)malloc(numtasks*sizeof(int));
  if (f->distance==NULL) {
    fprintf(stderr, "cannot allocate distance of size %lld , aborting ...\n", numtasks*sizeof(int));
    return(0);
  }
  f->mapping = (int *)malloc(MAPPING_DIM * numtasks*sizeof(int));
  if (f->mapping==NULL) {
    fprintf(stderr, "cannot allocate mapping_x of size %lld , aborting ...\n", MAPPING_DIM * numtasks*sizeof(int));
    return(0);
  }
  f->ptimings_stepcnt = (long *)malloc(opt_steps*sizeof(long));
  if (f->ptimings_stepcnt==NULL) {
    fprintf(stderr, "cannot allocate ptimings_stepcnt of size %lld , aborting ...\n", opt_steps*sizeof(long));
    return(0);
  }
  f->ptimings_step = (double *)malloc(opt_steps*sizeof(double));
  if (f->ptimings_step==NULL) {
    fprintf(stderr, "cannot allocate ptimings_step of size %lld , aborting ...\n", opt_steps*sizeof(double));
    return(0);
  }
  for(t=0;t<opt_steps;t++) {f->ptimings_step[t]=0.0;f->ptimings_stepcnt[t]=0;}
  f->timings_min=1e20;f->timings_max=-1e20;
  f->alllocation = (char *)malloc(numtasks*MAX_LOCATION_LENGTH);
  if (f->alllocation==NULL) {
    fprintf(stderr, "cannot allocate alllocation of size %d , aborting ...\n", numtasks*MAX_LOCATION_LENGTH);
    return(0);
  }
  for(t=0;t<numtasks*MAX_LOCATION_LENGTH;t++) {f->alllocation[t]='\0';}
#ifdef USEMASK
  f->mask = (char *)malloc(numtasks*numtasks*sizeof(char));
  if (f->mask==NULL) {
    fprintf(stderr, "cannot allocate mask of size %d , aborting ...\n", numtasks*numtasks);
    return(0);
  }
  for(t=0;t<numtasks*numtasks;t++) {f->mask[t]=0;}
#endif

  f->nodeslist_fromavg = (double *)malloc(numtasks*sizeof(double));
  if (f->nodeslist_fromavg==NULL) {
    fprintf(stderr, "cannot allocate nodeslist_fromavg of size %lld , aborting ...\n", numtasks*sizeof(double));
    return(1);
  } 

  f->nodeslist_toavg = (double *)malloc(numtasks*sizeof(double));
  if (f->nodeslist_toavg==NULL) {
    fprintf(stderr, "cannot allocate nodeslist_toavg of size %lld , aborting ...\n", numtasks*sizeof(double));
    return(1);
  } 

  f->nodeslist_frommax = (double *)malloc(numtasks*sizeof(double));
  if (f->nodeslist_frommax==NULL) {
    fprintf(stderr, "cannot allocate nodeslist_frommax of size %lld , aborting ...\n", numtasks*sizeof(double));
    return(1);
  } 

  f->nodeslist_tomax = (double *)malloc(numtasks*sizeof(double));
  if (f->nodeslist_tomax==NULL) {
    fprintf(stderr, "cannot allocate nodeslist_tomax of size %lld , aborting ...\n", numtasks*sizeof(double));
    return(1);
  } 

  f->nodeslist_fromcntmax = (int *)malloc(numtasks*sizeof(int));
  if (f->nodeslist_fromcntmax==NULL) {
    fprintf(stderr, "cannot allocate nodeslist_fromcntmax of size %lld , aborting ...\n", numtasks*sizeof(int));
    return(1);
  } 

  f->nodeslist_tocntmax = (int *)malloc(numtasks*sizeof(int));
  if (f->nodeslist_tocntmax==NULL) {
    fprintf(stderr, "cannot allocate nodeslist_tocntmax of size %lld , aborting ...\n", numtasks*sizeof(int));
    return(1);
  } 

  for(t=0;t<numtasks;t++) {f->nodeslist_fromavg[t]=0;f->nodeslist_toavg[t]=0.0;}
  for(t=0;t<numtasks;t++) {f->nodeslist_frommax[t]=0;f->nodeslist_tomax[t]=0.0;}
  for(t=0;t<numtasks;t++) {f->nodeslist_fromcntmax[t]=0;f->nodeslist_tocntmax[t]=0.0;}

  return(1);
}

int read_maskfile(analysis_fields_t *fields, int numtasks, char *filename) {
  FILE *infile;
  long from,to;
  long count;

#ifdef USEMASK
  printf("pingponganalysis: reading mask file %s ...\n",filename);
  infile=fopen(filename,"r");
  if (!infile) {
    fprintf(stderr, "cannot open file %s , aborting ...\n", filename);
    return(0);
  }

  count=0;
  while(fscanf(infile,"%ld %ld",&from,&to)==2) {
    if (
	(from>=0) && (from<numtasks) && 
	(to>=0) && (to<numtasks)
	) {
      fields->mask[from*numtasks+to]=1;
      count++;
    } else {
      printf("MASK %4ld %4ld out of range (%d)\n",from,to,numtasks);
    }

  }
  fclose(infile);
  printf("pingponganalysis: reading mask file %s ... %ld pairs found\n",filename,count);

#endif  
  return(1);
}
