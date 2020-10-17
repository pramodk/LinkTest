/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2011                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/
#ifndef LIMPRIOQUEUE_H
#define LIMPRIOQUEUE_H
typedef struct node {
  double time;
  int from;
  int to;
} node ;
typedef node * node_ptr;
typedef struct limqueue {
  int qsize;
  int max_elems;
  int sorted;
  node_ptr elements;
} lim_queue ; 
void limqueue_initalize(lim_queue *a,int nodes);
void limqueue_reinitalize(lim_queue *a);
int limqueue_merge_in(lim_queue *a, lim_queue *b);
int  limqueue_insert(lim_queue *a, double time, int from, int to);
int  limqueue_print(lim_queue *a);
int  limqueue_unsort_insert(lim_queue *a, double time, int from, int to);
int limqueue_sort(lim_queue *a);
#endif
