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

#include "mpilinktest_util.h"
#include "mpilinktest_datastructures_par.h"


int _linktest_distribute_spec_mpi( _linktest_spec *linktest_spec ) {

  int       rc = 1;

  if (linktest_spec == NULL) {
    return(_linktest_errorprint(0,_LINKTEST_ERROR_RETURN,"_linktest_distribute_spec: cannot distribute, data structure is not allocated, aborting ...\n"));
  }
  MPI_Bcast(&linktest_spec->do_alltoall,        1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->do_mix,             1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->do_serial,          1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->do_sion,            1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->do_sion_par,        1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->sleeptime,          1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->collectpnum,        1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->number_of_messages, 1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->number_of_warmup_messages, 1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->length_of_message,  1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->control_wtime,      1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->control_niter,      1, MPI_INT,    0, MPI_COMM_WORLD);
  MPI_Bcast(&linktest_spec->max_stest,          1, MPI_INT,    0, MPI_COMM_WORLD);

  return (rc);
}
