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
#include "limprioqueue.h"

void  limqueue_initalize(lim_queue *a,int nodes) {
  a->qsize = 0; 
  a->sorted = 1;
  a->max_elems = nodes;
  a->elements = (node_ptr)malloc(sizeof(node)*((a->max_elems)));
  if (a->elements == NULL) {
    fprintf(stderr, "cannot allocate lqueue of size %ld, aborting ...\n", (unsigned long) sizeof(node)*((a->max_elems)));
    exit(1);
  }
}

void  limqueue_reinitalize(lim_queue *a) {
  a->qsize = 0; 
  a->sorted = 1;
}

int limqueue_insert(lim_queue *a, double time, int from, int to) {
  int i, inspos,lpos=-1,last;

  if(!a->sorted) {
    return(limqueue_unsort_insert(a,time,from,to));
  }

  /* check if new element will be stored in queue */
  if(a->qsize==a->max_elems) {
    if(time<a->elements[a->qsize-1].time) return(0);
  }

  /* find first element which is smaller */
  inspos=a->qsize-1;
  last=0;
  for(i=0;((i<a->qsize) && (!last)); i++) {
    if(a->elements[i].time<time) {
      inspos=i;
      last=1;
    }
  }
  /* move data behind insert pos one forward */
  if(inspos<a->qsize-1) {
    if(a->qsize<a->max_elems) lpos=a->qsize-1;
    else                      lpos=a->max_elems-1;
    for(i=lpos;i>=inspos; i--) {
      a->elements[i+1].time=a->elements[i].time;
      a->elements[i+1].from=a->elements[i].from;
      a->elements[i+1].to  =a->elements[i].to;
    }
  } 
  printf("Insert at %d  time=%10.6f lpos=%d\n",inspos,time,lpos);

  /* set new values */ 
  if(a->qsize<a->max_elems) a->qsize++; 
  a->elements[inspos].time=time;
  a->elements[inspos].from=from;
  a->elements[inspos].to  =to;
  return(1);
}

int limqueue_unsort_insert(lim_queue *a, double time, int from, int to) {
  int i, inspos,lpos=-1,last;

  /* check if queue is not full */
  if(a->qsize==a->max_elems) {
    fprintf(stderr,"cannot add element %f,%d,%d to queue, reason: queue full, please use limqueue_insert\n",time,from,to);
    return(0);
  }
  a->sorted=0;

  /* set new values */ 
  inspos=a->qsize;
  a->elements[inspos].time=time;
  a->elements[inspos].from=from;
  a->elements[inspos].to  =to;
  a->qsize++; 
  return(1);
}

int _limqueue_sort(const void *aptr, const void *bptr) {
  node_ptr a;
  node_ptr b;
  int rc;
  a=(node_ptr) aptr;
  b=(node_ptr) bptr;
  if(a->time>b->time) rc=-1;
  else {
    if(a->time == b->time) rc=0;
    else                   rc=1;
  }
/*   printf("   _limqueue_sort(%f,%f) -> rc=%d\n",a->time,b->time,rc); */
  return(rc);
}

int limqueue_sort(lim_queue *a) {
  int i, inspos,lpos=-1,last;

  if(a->sorted) {
    return(1);
  }

/*   printf("Sorting queue, %d elements\n",a->qsize); */
  qsort(a->elements,a->qsize,sizeof(node),_limqueue_sort);
  a->sorted=1;
  return(1);
}

int limqueue_merge_in(lim_queue *a, lim_queue *b) {
  int num_remove_nodes,nr ,last_a, last_b, ins_a, copy_from_pos, copy_to_pos;
  lim_queue *from_q;

/*   printf("Merging queues, queue a:%d elements, queue b:%d elements\n",a->qsize,b->qsize); */

  nr=num_remove_nodes=a->qsize+b->qsize - a->max_elems;
/*   printf("  -> nodes to remove: %d elements\n",num_remove_nodes); */

  /* find position of last staying element in each queue */
  last_a=a->qsize-1;
  last_b=b->qsize-1;
  while(nr>0) {
    if ( (last_a>=0) && (last_b>=0)) {
      if(a->elements[last_a].time>b->elements[last_b].time) {
	last_b--;
      } else {
	last_a--;
      }
    } else {
      /* one queue has no more elements to remove */
      if (last_a>=0) last_a--;
      if (last_b>=0) last_b--;
    }
    nr--;
  }
/*   printf("  -> last positions: queue a: %d , queue b: %d\n",last_a, last_b); */
 
  if(last_b>=0) {
    /* there are nodes to merge in */
    if(num_remove_nodes>0) {
      /* more than max_elems */
      ins_a=a->max_elems-1;
      a->qsize=a->max_elems;
    } else {
      ins_a=a->qsize+b->qsize-1;
      a->qsize=a->qsize+b->qsize;
    }
/*     printf("  -> starting insert  position in queue a: %d\n",ins_a); */
    
    /* merging */
    while((last_a>=0) || (last_b>=0)) {
/*       printf("  scan A[%2d] <-> B[%2d] ",last_a,last_b); */
      if ( (last_a>=0) && (last_b>=0)) {
	if(a->elements[last_a].time > b->elements[last_b].time) {
	  copy_to_pos  =ins_a;	  copy_from_pos=last_b; from_q=b;
/* 	  printf("   1 copy B[%d] -> A[%d]\n",last_b,ins_a); */
	  ins_a--;last_b--;
	} else {
	  copy_to_pos  =ins_a;	  copy_from_pos=last_a; from_q=a;
/* 	  printf("   1 copy A[%d] -> A[%d]\n",last_a,ins_a); */
	  ins_a--;last_a--;
	}
      } else {
	/* one queue has no more elements */
	if (last_a>=0) {
	  copy_to_pos  =ins_a;	  copy_from_pos=last_a; from_q=a;
/* 	  printf("   2 copy A[%d] -> A[%d]\n",last_a,ins_a); */
	  ins_a--;last_a--;
	}
	if (last_b>=0) {
	  copy_to_pos  =ins_a;	  copy_from_pos=last_b; from_q=b;
/* 	  printf("   2 copy B[%d] -> A[%d]\n",last_b,ins_a); */
	  ins_a--;last_b--;
	}
      }
      a->elements[copy_to_pos].time=from_q->elements[copy_from_pos].time;
      a->elements[copy_to_pos].from=from_q->elements[copy_from_pos].from;
      a->elements[copy_to_pos].to  =from_q->elements[copy_from_pos].to;
    }
  }
  return(1);
}

int  limqueue_print(lim_queue *a) {
  int i;
  printf("LQ: qs=%3d/%3d:",a->qsize,a->max_elems);
  for(i=0; i<a->qsize; i++) printf("%3d: [%10.4fus %3d %3d] ",i,a->elements[i].time*1000*1000,a->elements[i].from,a->elements[i].to);
  printf("\n");
  return(1);
}
