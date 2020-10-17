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

#ifdef LINKTEST_BGQ
#include <spi/include/kernel/memory.h>
#endif

#ifdef LINKTEST_LINUX
#include <sys/times.h>
#include <sys/resource.h>

struct rusage rusagestr;
long get_memusage() {
  if (getrusage(RUSAGE_SELF, &rusagestr) != 0) {
    return -1;
  }
  
  /* printf("get_memusage: rusage.ru_maxrss=%10ld\n",rusagestr.ru_maxrss);  */
  return rusagestr.ru_maxrss;
} 

int print_memusage(int i)
{
  if (getrusage(RUSAGE_SELF, &rusagestr) != 0)
    return -1;

  printf("PE: %04d: MEMUSAGE: %12.8f\n",i, rusagestr.ru_maxrss / 1024.0);

  return 1;
} 

#elif defined( LINKTEST_BGQ)
long get_memusage() {
  uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;

  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);
  return (long) heap;
} 

int print_memusage(int pe)
{
  uint64_t shared, persist, heapavail, stackavail, stack, heap, guard, mmap;

  Kernel_GetMemorySize(KERNEL_MEMSIZE_SHARED, &shared);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_PERSIST, &persist);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAPAVAIL, &heapavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACKAVAIL, &stackavail);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_STACK, &stack);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &heap);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_GUARD, &guard);
  Kernel_GetMemorySize(KERNEL_MEMSIZE_MMAP, &mmap);
  
  printf("PE: %04d: Allocated heap,  %10.5f MB, avail. heap:  %10.5f MB\n", pe, 
	 (double)heap/(1024*1024), (double)heapavail/(1024*1024));
  printf("PE: %04d: Allocated stack: %10.5f MB, avail. stack: %10.5f MB\n", pe, 
	 (double)stack/(1024*1024), (double)stackavail/(1024*1024));
  printf("PE: %04d: Memory: shared:  %10.5f MB, persist:      %10.5f MB, guard: %.2f MB, mmap: %.2f MB\n", pe, 
	 (double)shared/(1024*1024), (double)persist/(1024*1024), (double)guard/(1024*1024), (double)mmap/(1024*1024));
  return 1;
} 

#else

long get_memusage() {
  return 0;
} 

int print_memusage(int i)
{
  printf("PE: %04d: MEMUSAGE: ???\n",i);
  return 1;
} 


#endif
