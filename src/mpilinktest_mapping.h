/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2009                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/
#ifndef LINKTEST_MAPPING_H
#define LINKTEST_MAPPING_H
#define MAX_LOCATION_LENGTH 256
int  get_mapping (char **locationptr, int **mappingptr);
int write_mapping (char *filename);
/* int write_cart_mapping (char *filename, MPI_Comm cartesian_comm, int dim,int *dims); */
#endif
