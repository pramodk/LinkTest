/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2011                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/
#ifndef LINKTEST_UTIL_H
#define LINKTEST_UTIL_H
 /* Error types */
#define _LINKTEST_ERROR_RETURN         -10001
#define _LINKTEST_ERROR_ABORT          -10002
#define _LINKTEST_ERROR_WARN           -10003
#define _LINKTEST_ERROR_UNKNOWN        -10020

#include "mpilinktest_datastructures.h"

double time2mbs( double time, _linktest_spec *linktest_spec);
int _linktest_errorprint(int rc, int level, const char *format, ...);
int collective_print_gather(char *cbuffer);
int collective_print(char *cbuffer);
int getpartnerplan_pe (int penr, int n, int *pepartner, int mix, int my_rank);
int getrandomize_vector_search (int n, int *vec);
int getrandomize_vector_swap (int n, int *vec);
int getpartneropt (int penr, int n, int *vec, int *pepartner);
int perm (int step, int n, int *line);

#endif
