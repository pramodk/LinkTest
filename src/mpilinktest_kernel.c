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

#include <mpi.h>

#include "mpilinktest_datastructures.h"
#include "mpilinktest_kernel.h"


double pingpong ( int from,  int to, _linktest_spec *linktest_spec ) {

  double timeused=0.0;

  /* warm up */
  _pingpong(linktest_spec->rank, from, to, linktest_spec->length_of_message, linktest_spec->number_of_warmup_messages);

  /* run */
  timeused=_pingpong(linktest_spec->rank, from, to, linktest_spec->length_of_message, linktest_spec->number_of_messages);

  return(timeused);
} 


double _pingpong ( int rank, int from,  int to, int length_of_message, int number_of_messages ) {

  double timeused=0.0;
  char buffer[LENGTH_OF_MESSAGE_MAX];
  int i;   
  double start, finish, time;
  MPI_Status status;

  /* printf("pingpong on %d:  from: %4d To: %4d iterations=%d\n",my_rank,from,to,iterations); */

  for (i = 1; i <= length_of_message; i++) buffer[i]=(char) (i%256);
  
  if (from==to) return(0.0);

/*   MPI_Barrier(MPI_COMM_WORLD);  */
  if (rank == from) 
  {

    start = MPI_Wtime();
     
    for (i = 1; i <= number_of_messages; i++)
    {
      MPI_Send(buffer, length_of_message, MPI_BYTE, to, PING,
	       MPI_COMM_WORLD);

      MPI_Recv(buffer, length_of_message, MPI_BYTE, to, PONG,
	       MPI_COMM_WORLD, &status);
    }

    finish = MPI_Wtime();

    time = finish - start;
    timeused=(float)(time / (2 * number_of_messages));
  }
  else if ( rank == to) 
  {
    for (i = 1; i <= number_of_messages; i++)
    { 
      MPI_Recv(buffer, length_of_message, MPI_BYTE, from, PING,
	       MPI_COMM_WORLD, &status);

      MPI_Send(buffer, length_of_message, MPI_BYTE, from, PONG,
	       MPI_COMM_WORLD);
    }
  }
/*   MPI_Barrier(MPI_COMM_WORLD);  */

/*   printf("pingpong on %d:  from: %4d To: %4d ready \n",my_rank,from,to); */

  return(timeused);
} 

double alltoall (_linktest_spec *linktest_spec) {

  double timeused=0.0;
  char in[length_of_message_a2a];
  char out[length_of_message_a2a];
  int i,maxlen;   
  double start, finish, time;
  /*   maxlen=length_of_message_a2a/tasks; */
  maxlen=1;
  for(i=0;i<linktest_spec->size*maxlen;i++)  out[i] = (char) (linktest_spec->rank%256);
  /*   if(my_rank<10) { */
  /*     printf("starting MPI_Alltoall: maxlen=%d tasks=%d\n",maxlen,tasks); */
  /*   } */
  start = MPI_Wtime();
  MPI_Alltoall(out,maxlen,MPI_BYTE,in,maxlen,MPI_BYTE,MPI_COMM_WORLD);
  finish = MPI_Wtime();
  
  time = finish - start;
  timeused=(float)(time);
  /* '!!! ckeck results !!! */
  return(timeused);
}
