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
#include "limprioqueue.h"

int main(int argc, char **argv) {
  lim_queue lqueue_a;
  lim_queue lqueue_b;
  int i;

  
  printf("\n-------- QUEUE A ---------- \n");
  limqueue_initalize(&lqueue_a, 5);
  
  limqueue_insert(&lqueue_a, 0.5,  0,1); limqueue_print(&lqueue_a);
  limqueue_insert(&lqueue_a, 0.7,  0,2); limqueue_print(&lqueue_a);
  limqueue_insert(&lqueue_a, 0.35, 0,3); limqueue_print(&lqueue_a);
  limqueue_insert(&lqueue_a, 0.1,  0,4); limqueue_print(&lqueue_a); 
  limqueue_insert(&lqueue_a, 0.9,  0,5); limqueue_print(&lqueue_a); 
  limqueue_insert(&lqueue_a, 0.3,  0,6); limqueue_print(&lqueue_a);  
  
  printf("\n-------- QUEUE B ---------- \n");
  limqueue_initalize(&lqueue_b, 5);
  limqueue_unsort_insert(&lqueue_b, 0.8, 0,1); limqueue_print(&lqueue_b);
  limqueue_unsort_insert(&lqueue_b, 0.6, 0,2); limqueue_print(&lqueue_b);
  limqueue_unsort_insert(&lqueue_b, 0.4, 0,3); limqueue_print(&lqueue_b); 
  limqueue_unsort_insert(&lqueue_b, 0.2, 0,4); limqueue_print(&lqueue_b); 
  limqueue_unsort_insert(&lqueue_b, 0.0, 0,4); limqueue_print(&lqueue_b); 
  limqueue_sort(&lqueue_b);limqueue_print(&lqueue_b);
  
  printf("\n-------- QUEUE MERGE ---------- \n");
  limqueue_print(&lqueue_a); 
  limqueue_print(&lqueue_b); 
  limqueue_merge_in(&lqueue_a, &lqueue_b);
  limqueue_print(&lqueue_a);
}
