/****************************************************************************
**  LinkTest                                                               **
*****************************************************************************
**  Copyright (c) 2008-2011                                                **
**  Forschungszentrum Juelich, Juelich Supercomputing Centre               **
**                                                                         **
**  See the file COPYRIGHT in the package base directory for details       **
****************************************************************************/
#ifndef MPILINKTEST_KERNEL_H
#define MPILINKTEST_KERNEL_H

double pingpong ( int from,  int to, _linktest_spec *linktest_spec);
double _pingpong ( int rank,  int from,  int to, int length_of_message, int number_of_messages );

double alltoall ( _linktest_spec *linktest_spec);

#endif
